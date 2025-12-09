#!/usr/bin/env Rscript
# 06_parse_hmmer.R
# Parse HMMER output for a given SRA ID, extract only the referenced IMG/VR
# cluster sequences from a single monolithic FASTA (optionally .gz),
# then run MAFFT → trimAl → FastTree, while tarring intermediates to keep inode count low.
#
# Usage:
#   Rscript 06_parse_hmmer.R <region> <sra_id> [big_faa_path]
# Example:
#   Rscript 06_parse_hmmer.R marine_eq SRR000001 /home/vaibhavk/scratch/IMG_VR/all_clusters.faa
#
# Requirements in PATH: mafft, trimal, FastTree, awk, tar, gzip
# R packages: seqinr

suppressPackageStartupMessages(library(seqinr))
options(warn = 1)

# ----------------------------- #
# 1) Args and directories       #
# ----------------------------- #
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 06_parse_hmmer.R <region> <sra_id> [big_faa_path]")
}
region <- args[1]
sra_id <- args[2]

# Where CheckV outputs live (persistent)
base_input_dir  <- file.path("/home/vaibhavk/scratch/checkv_results", region)

# Use node-local scratch for heavy I/O
work_root <- Sys.getenv("SLURM_TMPDIR",
               unset = Sys.getenv("TMPDIR",
                        unset = "/dev/shm"))
sra_output_dir <- file.path(work_root, paste0("hmmer_", sra_id))
dir.create(sra_output_dir, recursive = TRUE, showWarnings = FALSE)

# Monolithic IMG/VR FASTA (persistent, but we only read it)
clusters_root <- "/home/vaibhavk/scratch/IMG_VR"
big_faa <- if (length(args) >= 3) args[3] else file.path(clusters_root, "all_clusters.faa")

# Inputs from your pipeline (persistent)
hmm_outfile  <- file.path(base_input_dir, "MGM_Output", paste0(sra_id, "_hmm.out"))
mgm_faa_file <- file.path(base_input_dir, "MGM_INTRM",  paste0(sra_id, ".faa"))

# Subdirs (all on local scratch)
for (d in c("hits_faa","hits_aln","hits_trimaln","hits_trees")) {
  dir.create(file.path(sra_output_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# TAR bundles to keep inode count low (on local scratch)
tar_created <- list(faa=FALSE, aln=FALSE, trimaln=FALSE)
tar_files   <- list(
  faa     = file.path(sra_output_dir, "hits_faa.tar"),
  aln     = file.path(sra_output_dir, "hits_aln.tar"),
  trimaln = file.path(sra_output_dir, "hits_trimaln.tar")
)
for (tf in tar_files) if (file.exists(tf)) file.remove(tf)

# Log
log_file <- file.path(sra_output_dir, "hmmer_processing.log")
log_msg <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ..., "\n",
                             file = log_file, append = TRUE, sep = " ")

log_msg("Processing started for", sra_id)
log_msg("Work root:", work_root)

# ----------------------------- #
# 2) Sanity checks              #
# ----------------------------- #
need <- c("mafft","trimal","awk","tar","gzip","FastTree")
miss <- need[Sys.which(need) == ""]
if (length(miss)) stop("Missing in PATH: ", paste(miss, collapse=", "))

if (!file.exists(hmm_outfile)) stop("Missing HMMER output: ", hmm_outfile)
if (!file.exists(mgm_faa_file)) stop("Missing MGM FAA: ", mgm_faa_file)
if (!file.exists(big_faa))     stop("Missing monolithic IMG/VR FASTA: ", big_faa)

# Threads
threads <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK","1")))
if (is.na(threads) || threads < 1) threads <- 1

# ----------------------------- #
# 3) Parse HMMER output         #
# ----------------------------- #
in_lines   <- readLines(hmm_outfile, warn = FALSE)
hit_name   <- list() # cluster_id -> vector of mgm seq ids
Evalues    <- list()
cluster_nm <- NA_character_
in_scores  <- FALSE

trimws1 <- function(x) trimws(x, which = "both")

for (line in in_lines) {
  # "Query:       1000009  [M=792]"
  if (grepl("^Query:\\s+", line)) {
    parts <- strsplit(trimws1(line), "\\s+")[[1]]
    cluster_nm <- parts[2]
    in_scores  <- FALSE
    log_msg("Now processing cluster", cluster_nm)
    next
  }
  if (grepl("^\\s*Scores for complete sequences", line)) { in_scores <- TRUE; next }
  if (in_scores && (grepl("^\\s*Domain annotation for each sequence", line) ||
                    grepl("^\\s*Alignments of top-scoring domains", line) ||
                    grepl("^\\s*Internal pipeline statistics", line))) {
    in_scores <- FALSE; next
  }
  if (in_scores) {
    tl <- trimws1(line)
    if (tl == "" || grepl("^-{2,}", tl) || grepl("^E-value\\s+score\\s+bias", tl)) next
    tokens <- strsplit(gsub("\\s+", " ", tl), " ", fixed = TRUE)[[1]]
    # E-value score bias E-value score bias exp N Sequence [Description...]
    if (length(tokens) >= 9) {
      e_val  <- tokens[1]
      seq_id <- tokens[9]                      # Sequence column
      seq_id <- strsplit(seq_id, "\\|")[[1]][1]
      if (is.null(hit_name[[cluster_nm]])) {
        hit_name[[cluster_nm]] <- seq_id
        Evalues[[cluster_nm]]  <- e_val
      } else {
        hit_name[[cluster_nm]] <- c(hit_name[[cluster_nm]], seq_id)
        Evalues[[cluster_nm]]  <- c(Evalues[[cluster_nm]],  e_val)
      }
    }
  }
}

# Deduplicate per cluster
if (length(hit_name)) for (k in names(hit_name)) hit_name[[k]] <- unique(hit_name[[k]])

cluster_names <- names(hit_name)
if (length(cluster_names) == 0) {
  log_msg("No clusters with hits found in", hmm_outfile, ". Exiting.")
  cat("Done. Processed: 0 Skipped: 0\n")
  quit(save = "no", status = 0)
}

# ----------------------------- #
# 4) Load MGM FAA and clean     #
# ----------------------------- #
log_msg("Loading MGM FAA:", mgm_faa_file)
mgm_faa <- read.fasta(mgm_faa_file, as.string = TRUE, seqtype = "AA", forceDNAtolower = FALSE)

first_token <- function(x) vapply(strsplit(x, "\\s+"), function(z) z[[1]], FUN.VALUE = character(1))
first_field <- function(x) vapply(strsplit(x, "\\|"),   function(z) z[[1]], FUN.VALUE = character(1))
mgm_names_raw   <- names(mgm_faa)
mgm_names_token <- first_token(mgm_names_raw)
mgm_faa_names   <- first_field(mgm_names_token)

# ----------------------------- #
# 5) Extract refs via AWK       #
# ----------------------------- #
ref_cache <- file.path(sra_output_dir, "ref_cache")
dir.create(ref_cache, showWarnings = FALSE, recursive = TRUE)

ids_file <- file.path(ref_cache, "cluster_ids.txt")
writeLines(cluster_names, ids_file)

awk_file <- file.path(ref_cache, "extract_refs.awk")
awk_src <- paste(
  'BEGIN{ while((getline k<IDS)>0){ gsub(/\\r$/,"",k); if(k!="") want[k]=1 } }',
  '/^>/{ h=$0; sub(/^>/,"",h); split(h,a,/\\|/); cid=a[1];',
  '      newout = (want[cid] ? DIR "/" cid ".faa" : "");',
  '      if(out!="" && newout!=out){ close(out) }',
  '      out=newout }',
  '{ if(out!="") print >> out }',
  'END{ if(out!="") close(out) }',
  sep="\n"
)
writeLines(awk_src, awk_file)

log_msg("Extracting", length(cluster_names), "clusters from big FASTA via AWK")
is_gz <- grepl("\\.gz$", big_faa, ignore.case = TRUE)
cmd <- if (is_gz) {
  sprintf("bash -lc 'gzip -cd %s | awk -v IDS=%s -v DIR=%s -f %s'",
          shQuote(big_faa), shQuote(ids_file), shQuote(ref_cache), shQuote(awk_file))
} else {
  sprintf("awk -v IDS=%s -v DIR=%s -f %s %s",
          shQuote(ids_file), shQuote(ref_cache), shQuote(awk_file), shQuote(big_faa))
}
rs <- system(cmd)
if (rs != 0) stop("AWK extraction failed with code ", rs)

# ----------------------------- #
# 6) Per-cluster processing     #
# ----------------------------- #
min_hits <- 1L
n_done   <- 0L
n_skip   <- 0L

for (cid in cluster_names) {
  hits <- hit_name[[cid]]
  if (length(hits) < min_hits) { n_skip <- n_skip + 1L; next }

  log_msg("Cluster", cid, "with", length(hits), "MGM hits")

  # Reference cluster from cache
  ref_path <- file.path(ref_cache, paste0(cid, ".faa"))
  if (!file.exists(ref_path) || file.info(ref_path)$size == 0) {
    log_msg("Missing or empty ref cache for cluster", cid, "- skipping")
    n_skip <- n_skip + 1L; next
  }
  ref_faa <- read.fasta(ref_path, as.string = TRUE, seqtype = "AA", forceDNAtolower = FALSE)

  # MGM sequences that matched
  idx <- which(mgm_faa_names %in% hits)
  if (length(idx) == 0) {
    log_msg("No MGM sequences found in FAA for cluster", cid, "- skipping")
    file.remove(ref_path)
    n_skip <- n_skip + 1L; next
  }

  combo    <- c(ref_faa, mgm_faa[idx])
  combo_nm <- c(names(ref_faa), names(mgm_faa)[idx])

  # Combined FASTA (local scratch)
  out_faa <- file.path(sra_output_dir, "hits_faa", paste0(cid, ".faa"))
  write.fasta(sequences = combo, names = combo_nm, file.out = out_faa)

  # MAFFT (threads + auto)
  aln_fa <- file.path(sra_output_dir, "hits_aln", paste0(cid, ".fasta"))
  cmd_mafft <- sprintf("mafft --thread %d --auto %s > %s",
                       threads, shQuote(out_faa), shQuote(aln_fa))
  rs <- system(cmd_mafft)
  if (rs != 0 || !file.exists(aln_fa) || file.info(aln_fa)$size == 0) {
    log_msg("MAFFT failed for cluster", cid, "- skipping")
    file.remove(out_faa, ref_path)
    n_skip <- n_skip + 1L; next
  }

  # trimAl
  tri_fa <- file.path(sra_output_dir, "hits_trimaln", paste0(cid, ".fasta"))
  cmd_trimal <- sprintf("trimal -in %s -out %s -gappyout", shQuote(aln_fa), shQuote(tri_fa))
  rs <- system(cmd_trimal)
  if (rs != 0 || !file.exists(tri_fa) || file.info(tri_fa)$size == 0) {
    log_msg("trimAl failed for cluster", cid, "- skipping")
    file.remove(out_faa, aln_fa, ref_path)
    n_skip <- n_skip + 1L; next
  }

  # FastTree (binary is guaranteed by sanity checks)
  tree_f <- file.path(sra_output_dir, "hits_trees", paste0(cid, ".trees"))
  cmd_ft <- sprintf("FastTree -lg -gamma %s > %s", shQuote(tri_fa), shQuote(tree_f))
  rs <- system(cmd_ft)
  if (rs != 0 || !file.exists(tree_f) || file.info(tree_f)$size == 0) {
    log_msg("FastTree failed for cluster", cid, "- skipping tree")
  }

  # Bundle + remove intermediates to keep inodes low
  # hits_faa
  rel <- file.path("hits_faa", paste0(cid, ".faa"))
  if (!tar_created$faa) {
    system(sprintf("tar -cf %s -C %s %s",
                   shQuote(tar_files$faa), shQuote(sra_output_dir), shQuote(rel)))
    tar_created$faa <- TRUE
  } else {
    system(sprintf("tar -rf %s -C %s %s",
                   shQuote(tar_files$faa), shQuote(sra_output_dir), shQuote(rel)))
  }
  file.remove(out_faa)

  # hits_aln
  rel <- file.path("hits_aln", paste0(cid, ".fasta"))
  if (!tar_created$aln) {
    system(sprintf("tar -cf %s -C %s %s",
                   shQuote(tar_files$aln), shQuote(sra_output_dir), shQuote(rel)))
    tar_created$aln <- TRUE
  } else {
    system(sprintf("tar -rf %s -C %s %s",
                   shQuote(tar_files$aln), shQuote(sra_output_dir), shQuote(rel)))
  }
  file.remove(aln_fa)

  # hits_trimaln
  rel <- file.path("hits_trimaln", paste0(cid, ".fasta"))
  if (!tar_created$trimaln) {
    system(sprintf("tar -cf %s -C %s %s",
                   shQuote(tar_files$trimaln), shQuote(sra_output_dir), shQuote(rel)))
    tar_created$trimaln <- TRUE
  } else {
    system(sprintf("tar -rf %s -C %s %s",
                   shQuote(tar_files$trimaln), shQuote(sra_output_dir), shQuote(rel)))
  }
  file.remove(tri_fa)

  # Drop cached ref
  file.remove(ref_path)

  n_done <- n_done + 1L
  if ((n_done + n_skip) %% 250 == 0) gc()
}

# ----------------------------- #
# 7) Compress TARs              #
# ----------------------------- #
for (tf in tar_files) {
  if (file.exists(tf)) {
    system(sprintf("gzip -f %s", shQuote(tf)))
  } else {
    log_msg("Skipping gzip, not found:", tf)
  }
}

# ----------------------------- #
# 8) Save session + done        #
# ----------------------------- #
save.image(file.path(sra_output_dir, paste0(sra_id, "_hmmer_out.RData")))
log_msg("Completed.", "Clusters processed:", n_done, "Skipped:", n_skip)
cat("Done. Processed:", n_done, "Skipped:", n_skip, "\n")
