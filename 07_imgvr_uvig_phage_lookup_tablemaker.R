#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
})

# --- INPUT / OUTPUT ---
tsv        <- "/home/vaibhavk/scratch/IMG_VR/IMG_VR_2022-12-19_7.1/IMGVR_all_Sequence_information-high_confidence.tsv"
out_detail <- "/home/vaibhavk/scratch/IMG_VR/imgvr_uvig_phage_lookup.csv"
out_min    <- "/home/vaibhavk/scratch/IMG_VR/imgvr_uvig_phage_lookup.min.csv"

# --- helpers (vectorized fixes) ---
extract_rank <- function(x, rank) {
  m <- regexpr(paste0("\\b", rank, "__[^;]+"), x, perl = TRUE)
  out <- rep(NA_character_, length(x))
  hit <- m > 0
  if (any(hit)) {
    ml <- attr(m, "match.length")
    out[hit] <- substr(x[hit], m[hit], m[hit] + ml[hit] - 1L)
    out[hit] <- sub("^[^_]+__", "", out[hit])
  }
  out
}
host_domain <- function(host_str) {
  host_str <- as.character(host_str)
  out <- rep(NA_character_, length(host_str))
  m <- regexpr("\\bd__[^;]+", host_str, perl = TRUE)
  hit <- m > 0
  if (any(hit)) {
    ml <- attr(m, "match.length")
    tmp <- substr(host_str[hit], m[hit], m[hit] + ml[hit] - 1L)  # e.g. "d__Bacteria"
    tmp <- sub("^d__", "", tmp)                                  # "Bacteria"
    tmp <- sub("\\s+.*$", "", tmp)                               # trim trailing words
    out[hit] <- tmp
  }
  out
}
`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
nz <- function(x) ifelse(is.na(x) | x == "", NA_character_, x)

# --- taxonomy rules ---
phage_classes  <- c("Caudoviricetes")
phage_orders   <- c("Caudovirales","Tubulavirales")
phage_families <- c(
  "Inoviridae","Microviridae","Tectiviridae",
  "Siphoviridae","Myoviridae","Podoviridae",
  "Autographiviridae","Herelleviridae","Drexlerviridae","Demerecviridae","Ackermannviridae",
  "Corticoviridae","Plasmaviridae","Cystoviridae","Leviviridae"
)

phylum_map <- list(
  Artverviricota      = "non_phage",
  Cossaviricota       = "non_phage",
  Cressdnaviricota    = "non_phage",
  Dividoviricota      = "phage",
  Duplornaviricota    = "non_phage",
  Hofneiviricota      = "non_phage",
  Kitrinoviricota     = "non_phage",
  Lenarviricota       = "unknown",       # treated as unknown unless deeper taxonomy shows phage
  Negarnaviricota     = "non_phage",
  Nucleocytoviricota  = "non_phage",
  Peploviricota       = "non_phage",
  Phixviricota        = "phage",
  Pisuviricota        = "non_phage",
  Preplasmiviricota   = "phage",
  Saleviricota        = "non_phage",
  Taleaviricota       = "non_phage",
  Uroviricota         = "phage"
)

# --- load TSV ---
dt <- fread(tsv, sep = "\t", quote = "", na.strings = c("NA",""), showProgress = TRUE)
setnames(dt, old = names(dt), new = gsub("\\s+", "_", names(dt)))

# parse taxonomy
dt[, host_domain  := host_domain(Host_taxonomy_prediction)]
dt[, viral_phylum := extract_rank(Taxonomic_classification, "p")]
dt[, viral_class  := extract_rank(Taxonomic_classification, "c")]
dt[, viral_order  := extract_rank(Taxonomic_classification, "o")]
dt[, viral_family := extract_rank(Taxonomic_classification, "f")]
dt[, phylum_core  := sub("^p__", "", nz(viral_phylum))]

# 1) host-based
dt[, phage_flag := fifelse(host_domain %in% c("Bacteria","Archaea"), "phage",
                      fifelse(host_domain %in% c("Eukaryota"), "non_phage", NA_character_))]
dt[, flag_method := fifelse(!is.na(phage_flag), "host", NA_character_)]
dt[, reason      := fifelse(!is.na(phage_flag), paste0("host_domain=", host_domain), NA_character_)]

# 2) phylum-based
need_phy <- which(is.na(dt$phage_flag) & !is.na(dt$phylum_core))
if (length(need_phy)) {
  lab <- vapply(dt$phylum_core[need_phy], function(p) phylum_map[[p]] %||% NA_character_, character(1))
  ix_assign <- need_phy[!is.na(lab) & lab != "mixed"]
  if (length(ix_assign)) {
    dt$phage_flag[ix_assign]  <- lab[!is.na(lab) & lab != "mixed"]
    dt$flag_method[ix_assign] <- "phylum"
    dt$reason[ix_assign]      <- paste0("phylum=", dt$phylum_core[ix_assign])
  }
  dt[, mixed_phylum := fifelse(phylum_core %in% names(phylum_map)[unlist(phylum_map)=="mixed"], TRUE, FALSE)]
} else {
  dt[, mixed_phylum := FALSE]
}

# 3) taxonomy fallback
need_tax <- which(is.na(dt$phage_flag))
if (length(need_tax)) {
  ix <- need_tax
  tax_phage <- (!is.na(dt$viral_class[ix])  & dt$viral_class[ix]  %chin% phage_classes)  |
               (!is.na(dt$viral_order[ix])  & dt$viral_order[ix]  %chin% phage_orders)   |
               (!is.na(dt$viral_family[ix]) & dt$viral_family[ix] %chin% phage_families)
  hit <- ix[tax_phage]
  if (length(hit)) {
    pick_reason <- function(i) {
      if (!is.na(dt$viral_family[i]) && dt$viral_family[i] %chin% phage_families) return(paste0("tax_hit_family=", dt$viral_family[i]))
      if (!is.na(dt$viral_order[i])  && dt$viral_order[i]  %chin% phage_orders)   return(paste0("tax_hit_order=",  dt$viral_order[i]))
      if (!is.na(dt$viral_class[i])  && dt$viral_class[i]  %chin% phage_classes)  return(paste0("tax_hit_class=",  dt$viral_class[i]))
      return("tax_hit")
    }
    dt$phage_flag[hit]  <- "phage"
    dt$flag_method[hit] <- "taxonomy"
    dt$reason[hit]      <- vapply(hit, pick_reason, character(1))
  }
}

# 4) default unknown
dt[is.na(phage_flag) & mixed_phylum == TRUE,  `:=`(phage_flag = "unknown", flag_method = "none",
                                                   reason = "phylum_mixed=Lenarviricota")]
dt[is.na(phage_flag), `:=`(phage_flag = "unknown", flag_method = "none",
                           reason = ifelse(is.na(reason) | reason=="", "no_host_or_tax_phage_signal", reason))]

# export detailed
keep <- c("UVIG","phage_flag","flag_method","reason",
          "host_domain","Host_taxonomy_prediction",
          "viral_phylum","viral_class","viral_order","viral_family",
          "Confidence","Estimated_completeness")
out_dt <- dt[, ..keep]
fwrite(out_dt, out_detail, quote = TRUE)

# export minimal
out_min_dt <- out_dt[, .(UVIG, phage_flag)]
fwrite(out_min_dt, out_min, quote = TRUE)

# summary
cat("Wrote detailed:", out_detail, "\n")
cat("Wrote minimal :", out_min, "\n")
print(out_dt[, .N, by = .(phage_flag, flag_method)][order(-N)])
