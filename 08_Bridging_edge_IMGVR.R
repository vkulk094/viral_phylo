#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(data.table)
})
options(stringsAsFactors = FALSE)

# Usage: Rscript 08_Bridging_edge_IMGVR.R <hits_trees.tar|.tar.gz|.tar> <out_clades.csv>
args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript 08_Bridging_edge_IMGVR.R <hits_trees.tar|.tar.gz|.tar> <out_clades.csv>")
}
in_tar  <- normalizePath(args[1], mustWork=TRUE)
out_csv <- args[2]
members_file <- Sys.getenv("MEMBERS_FILE", unset = "")

# --- config / patterns ---
tree_ext_regex <- "\\.(tre|tree|trees|treefile|nwk|newick)$"
gene_pat  <- Sys.getenv("GENE_REGEX",  "(^|\\|)gene_")
imgvr_pat <- Sys.getenv("IMGVR_REGEX", "IMGVR")
uvig_re   <- "IMGVR_UViG_[0-9_]+"

# Cap MPD computation to avoid O(n^2) blowups on very large clades.
MAX_MPD_TIPS <- as.integer(Sys.getenv("MPD_MAX_TIPS", "400"))

# --- IMG/VR lookup ---
lk_path <- "/home/vaibhavk/scratch/IMG_VR/imgvr_uvig_phage_lookup.min.csv"
lk <- fread(lk_path, na.strings=c("","NA"))
if (!all(c("UViG","phage_flag") %in% names(lk))) {
  stop("Lookup missing required columns UViG, phage_flag: ", lk_path)
}
lk[, UViG := trimws(UViG)]
lk[, phage_flag := trimws(phage_flag)]
flag_map <- lk$phage_flag; names(flag_map) <- lk$UViG
rm(lk)

# --- helpers ---
extract_uvig <- function(x) {
  m <- regexpr(uvig_re, x, perl=TRUE)
  out <- rep(NA_character_, length(x))
  hit <- m > 0
  if (any(hit)) {
    ml <- attr(m, "match.length")
    out[hit] <- substr(x[hit], m[hit], m[hit] + ml[hit] - 1L)
  }
  trimws(out)
}
mrca_index_or_tip <- function(tr, tip_labels) {
  if (!length(tip_labels)) return(NA_integer_)
  if (length(tip_labels) == 1L) return(match(tip_labels[1], tr$tip.label))
  getMRCA(tr, tip_labels)
}
is_set_monophyletic <- function(tr, tip_labels) {
  if (length(tip_labels) <= 1L) return(TRUE)
  is.monophyletic(tr, tip_labels)
}
parent_vector <- function(tr) {
  N <- tr$Nnode + length(tr$tip.label)
  p <- integer(N); p[] <- NA_integer_
  ed <- tr$edge
  p[ed[,2]] <- ed[,1]
  p
}
lca_nodes <- function(parent, i, j) {
  if (is.na(i) || is.na(j)) return(NA_integer_)
  N <- length(parent)
  seen <- logical(N)
  ii <- i
  while (!is.na(ii) && ii > 0L && ii <= N && !seen[ii]) {
    seen[ii] <- TRUE
    ii <- parent[ii]
  }
  jj <- j
  while (!is.na(jj) && jj > 0L && jj <= N && !seen[jj]) {
    jj <- parent[jj]
  }
  if (is.na(jj) || jj <= 0L || jj > N) NA_integer_ else jj
}
distance_between_nodes_fast <- function(depth, parent, i, j) {
  if (is.na(i) || is.na(j)) return(NA_real_)
  lca <- lca_nodes(parent, i, j)
  if (is.na(lca)) return(NA_real_)
  depth[i] + depth[j] - 2 * depth[lca]
}

summarize_flags <- function(flags) {
  if (!length(flags)) return(list(nP=0L,nN=0L,nU=0L,flag="unknown"))
  nP <- sum(flags == "phage", na.rm=TRUE)
  nN <- sum(flags == "non_phage", na.rm=TRUE)
  nU <- sum(flags == "unknown" | is.na(flags))
  tree_flag <- if (nP>0 && nN==0) "phage"
          else if (nN>0 && nP==0) "non_phage"
          else if (nP==0 && nN==0) "unknown"
          else "mixed"
  list(nP=as.integer(nP), nN=as.integer(nN), nU=as.integer(nU), flag=tree_flag)
}

# --- output header ---
con <- file(out_csv, open="wt")
writeLines(paste(
  "tar","member","tree_index",
  "n_gene","n_imgvr","genes_monophyletic","imgvr_monophyletic","mrca_distance","notes",
  "n_imgvr_phage","n_imgvr_non_phage","n_imgvr_unknown","imgvr_tree_flag",
  "gene_pd","gene_mpd","gene_pairs","imgvr_pd","imgvr_mpd","imgvr_pairs",
  sep=","
), con)

# --- member set for this chunk ---
if (nzchar(members_file)) {
  members <- readLines(members_file, warn = FALSE)
  members <- members[nzchar(members)]
} else {
  members <- try(utils::untar(in_tar, list=TRUE), silent=TRUE)
  if (inherits(members, "try-error") || !length(members)) { close(con); quit(save="no", status=0) }
  members <- members[grepl(tree_ext_regex, members, ignore.case=TRUE)]
}
if (!length(members)) { close(con); quit(save="no", status=0) }

# --- stream one member from tar and parse ---
read_tree_from_tar <- function(tarfile, member) {
  # Use GNU tar to print the file to stdout. Use "--" to end options so member starting with '-' is safe.
  out <- tryCatch(
    system2("tar", c("-xOf", tarfile, "--", member), stdout = TRUE, stderr = NULL),
    error = function(e) character(0)
  )
  if (!length(out)) return(NULL)
  txt <- paste(out, collapse = "\n")
  tr <- try(read.tree(text = txt), silent = TRUE)
  if (inherits(tr, "try-error")) return(NULL)
  tr
}

process_member <- function(member_exact) {
  tr <- read_tree_from_tar(in_tar, member_exact)
  if (is.null(tr)) return(NULL)
  lst <- if (inherits(tr, "multiPhylo")) tr else list(tr)

  out <- vector("list", length(lst))
  for (k in seq_along(lst)) {
    tt <- lst[[k]]
    tips <- tt$tip.label
    if (!length(tips)) next

    gene_tips  <- tips[grep(gene_pat,  tips, perl=TRUE, ignore.case=TRUE)]
    imgvr_tips <- tips[grep(imgvr_pat, tips, perl=TRUE, ignore.case=TRUE)]
    n_gene  <- length(gene_tips)
    n_imgvr <- length(imgvr_tips)

    note <- character(0)
    # Predeclare new metrics.
    gene_pd <- NA_real_; gene_mpd <- NA_real_; gene_pairs <- NA_integer_
    imgvr_pd <- NA_real_; imgvr_mpd <- NA_real_; imgvr_pairs <- NA_integer_

    if (n_gene == 0L || n_imgvr == 0L) {
      note <- c(note, if (n_gene==0L) "no_gene" else NULL, if (n_imgvr==0L) "no_imgvr" else NULL)
      out[[k]] <- data.frame(
        tar=in_tar, member=member_exact, tree_index=k,
        n_gene=n_gene, n_imgvr=n_imgvr,
        genes_monophyletic=NA, imgvr_monophyletic=NA,
        mrca_distance=NA_real_, notes=paste(note, collapse=";"),
        n_imgvr_phage=NA_integer_, n_imgvr_non_phage=NA_integer_, n_imgvr_unknown=NA_integer_,
        imgvr_tree_flag=NA_character_,
        gene_pd=gene_pd, gene_mpd=gene_mpd, gene_pairs=gene_pairs,
        imgvr_pd=imgvr_pd, imgvr_mpd=imgvr_mpd, imgvr_pairs=imgvr_pairs,
        stringsAsFactors=FALSE
      )
      next
    }

    uvigs <- extract_uvig(imgvr_tips)
    flags <- flag_map[uvigs]
    flag_sum <- summarize_flags(unname(flags))

    mrca_g <- mrca_index_or_tip(tt, gene_tips)
    mrca_v <- mrca_index_or_tip(tt, imgvr_tips)
    depth  <- node.depth.edgelength(tt)
    parent <- parent_vector(tt)
    d_mv   <- distance_between_nodes_fast(depth, parent, mrca_g, mrca_v)

    mono_g <- is_set_monophyletic(tt, gene_tips)
    mono_v <- is_set_monophyletic(tt, imgvr_tips)

    # --- Within-clade distances (PD and MPD) ---
    if (n_gene >= 2L) {
      sub_g <- keep.tip(tt, gene_tips)
      if (!is.null(sub_g$edge.length)) {
        gene_pd <- sum(sub_g$edge.length)
        gene_pairs <- as.integer(n_gene * (n_gene - 1L) / 2L)
        if (n_gene <= MAX_MPD_TIPS) {
          Dg <- cophenetic.phylo(sub_g)
          gene_mpd <- mean(Dg[upper.tri(Dg)])
        } else {
          note <- c(note, sprintf("mpd_gene_capped_%d", n_gene))
        }
      }
    }
    if (n_imgvr >= 2L) {
      sub_v <- keep.tip(tt, imgvr_tips)
      if (!is.null(sub_v$edge.length)) {
        imgvr_pd <- sum(sub_v$edge.length)
        imgvr_pairs <- as.integer(n_imgvr * (n_imgvr - 1L) / 2L)
        if (n_imgvr <= MAX_MPD_TIPS) {
          Dv <- cophenetic.phylo(sub_v)
          imgvr_mpd <- mean(Dv[upper.tri(Dv)])
        } else {
          note <- c(note, sprintf("mpd_imgvr_capped_%d", n_imgvr))
        }
      }
    }

    out[[k]] <- data.frame(
      tar=in_tar, member=member_exact, tree_index=k,
      n_gene=n_gene, n_imgvr=n_imgvr,
      genes_monophyletic=mono_g, imgvr_monophyletic=mono_v,
      mrca_distance=as.numeric(d_mv), notes=paste(note, collapse=";"),
      n_imgvr_phage=flag_sum$nP, n_imgvr_non_phage=flag_sum$nN, n_imgvr_unknown=flag_sum$nU,
      imgvr_tree_flag=flag_sum$flag,
      gene_pd=gene_pd, gene_mpd=gene_mpd, gene_pairs=gene_pairs,
      imgvr_pd=imgvr_pd, imgvr_mpd=imgvr_mpd, imgvr_pairs=imgvr_pairs,
      stringsAsFactors=FALSE
    )
  }
  do.call(rbind, out)
}

# --- main loop over this chunk ---
for (m in members) {
  df <- process_member(m)
  if (!is.null(df) && nrow(df)) {
    write.table(df, con, sep=",", row.names=FALSE, col.names=FALSE, quote=TRUE)
  }
}
close(con)
