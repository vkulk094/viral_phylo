
#!/usr/bin/env Rscript

# --------------------- Package auto-install ---------------------
# Set AUTO_INSTALL=0 to disable. Override CRAN mirror with CRAN_REPO env var.
enable_auto_install <- tolower(Sys.getenv("AUTO_INSTALL", unset = "1")) %in% c("1","true","yes")
default_cran <- Sys.getenv("CRAN_REPO", unset = "https://cloud.r-project.org")

ensure_repos <- function() {
  r <- getOption("repos")
  if (is.null(r) || length(r) == 0L || identical(r, c(CRAN = "@CRAN@"))) {
    options(repos = c(CRAN = default_cran))
  } else if (is.na(r["CRAN"]) || r["CRAN"] == "" || r["CRAN"] == "@CRAN@") {
    r["CRAN"] <- default_cran
    options(repos = r)
  }
  invisible(getOption("repos"))
}

ensure_packages <- function(pkgs, repos = NULL) {
  if (!enable_auto_install) {
    message("AUTO_INSTALL disabled. Expecting packages to be present: ", paste(pkgs, collapse = ", "))
    return(invisible(FALSE))
  }
  if (is.null(repos)) repos <- ensure_repos()
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    message("Installing missing packages from ", paste(repos, collapse = ", "), ": ", paste(missing, collapse = ", "))
    tryCatch(
      install.packages(missing, repos = repos, dependencies = c("Depends","Imports","LinkingTo")),
      error = function(e) {
        stop("Failed to install packages: ", paste(missing, collapse = ", "), " | ", conditionMessage(e))
      }
    )
  }
  invisible(TRUE)
}

# Auto-install primary dependencies before loading
primary_pkgs <- c("dplyr","tidyr","readr","ggplot2","ggpubr","dunn.test","multcompView","tibble")
ensure_repos()
ensure_packages(primary_pkgs)


suppressMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggpubr)
  library(dunn.test)
  library(multcompView)
  library(tibble)
})

# ------------------------ Reproducibility knobs ------------------------
GLOBAL_SEED <- as.integer(Sys.getenv("FIG_SEED", unset = "20250827"))
set.seed(GLOBAL_SEED)
JITTER_BUDGET_REGION <- as.integer(Sys.getenv("JITTER_BUDGET_REGION", unset = "4000"))
JITTER_BUDGET_NINE   <- as.integer(Sys.getenv("JITTER_BUDGET_NINE",   unset = "500"))
JITTER_SEED_REGION   <- as.integer(Sys.getenv("JITTER_SEED_REGION",   unset = "123"))
JITTER_SEED_NINE     <- as.integer(Sys.getenv("JITTER_SEED_NINE",     unset = "321"))
message("Figure seed: ", GLOBAL_SEED,
        " | jitter(region): budget=", JITTER_BUDGET_REGION, " seed=", JITTER_SEED_REGION,
        " | jitter(nine): budget=", JITTER_BUDGET_NINE, " seed=", JITTER_SEED_NINE)

# ------------------------ Config ------------------------
base_dir <- "/home/vaibhavk/scratch/hmmer_align"

out_long  <- file.path(base_dir, "combined_patristic_long.csv")
out_stats <- file.path(base_dir, "combined_patristic_stats.csv")
out_pdf   <- file.path(base_dir, "region_metric_significance_boxplots.pdf")

# --- Cache for per-plot statistics and letters
cache_rds <- file.path(base_dir, "boxplot_cache.rds")
plot_cache <- new.env(parent = emptyenv())
plot_cache$region <- list()
plot_cache$nine   <- list()
plot_cache$meta <- list(
  script_version      = "v5.1",
  global_seed         = GLOBAL_SEED,
  jitter_budget_region= JITTER_BUDGET_REGION,
  jitter_budget_nine  = JITTER_BUDGET_NINE,
  jitter_seed_region  = JITTER_SEED_REGION,
  jitter_seed_nine    = JITTER_SEED_NINE
)

# Optional stricter filter for including gene_mpd and imgvr_mpd. 0 disables.
min_gene_tips <- as.integer(Sys.getenv("MIN_GENE_TIPS", unset = "2"))
min_imgvr_tips <- as.integer(Sys.getenv("MIN_IMGVR_TIPS", unset = "2"))
message("MIN_GENE_TIPS = ", min_gene_tips, " (0 means disabled)")
message("MIN_IMGVR_TIPS = ", min_imgvr_tips, " (0 means disabled)")
message("Filtering to monophyletic only: genes_monophyletic == TRUE AND imgvr_monophyletic == TRUE")

# Regions and biome mapping
regions <- c("marine_eq","marine_south","marine_north",
             "terrestrial_eq","terrestrial_north","terrestrial_south")
biome_of <- function(region) if (startsWith(region, "marine_")) "marine" else "terrestrial"

# Fixed region orders for biomes: North -> Eq -> South
marine_levels      <- c("marine_north", "marine_eq", "marine_south")
terrestrial_levels <- c("terrestrial_north", "terrestrial_eq", "terrestrial_south")

# ------------------- Defensive helpers -------------------
required_headers <- c(
  "n_gene","n_imgvr","genes_monophyletic","imgvr_monophyletic",
  "mrca_distance","notes","n_imgvr_phage","n_imgvr_non_phage",
  "n_imgvr_unknown","imgvr_tree_flag","gene_pd","gene_mpd",
  "gene_pairs","imgvr_pd","imgvr_mpd","imgvr_pairs"
)

check_csv_headers <- function(fp) {
  hdr <- tryCatch({
    names(readr::read_csv(fp, n_max = 0, show_col_types = FALSE, progress = FALSE))
  }, error = function(e) NULL)
  if (is.null(hdr)) {
    warning("Unreadable CSV (cannot read header): ", fp)
    return(FALSE)
  }
  miss <- setdiff(required_headers, hdr)
  if (length(miss)) {
    warning("Skipping file due to missing columns (", paste(miss, collapse = ","), "): ", fp)
    return(FALSE)
  }
  TRUE
}

# Ensure ggpubr functions exist (install hints instead of cryptic errors)
ensure_ggpubr <- function() {
  if (!("ggpubr" %in% loadedNamespaces())) return(invisible(TRUE))
  req <- c("annotate_figure","text_grob")
  for (fn in req) {
    if (!exists(fn, where = asNamespace("ggpubr"), inherits = FALSE)) {
      stop("ggpubr::", fn, " not found. Please reinstall ggpubr.")
    }
  }
  invisible(TRUE)
}

# ------------------- Discover files ---------------------
all_csv <- list.files(base_dir,
                      pattern = "^hits_trees_clade_distance_patristic\\.csv$",
                      recursive = TRUE, full.names = TRUE)

keep <- vapply(all_csv, function(p) {
  parts <- strsplit(p, "/", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(FALSE)
  region <- parts[length(parts)-2]
  region %in% regions
}, logical(1))

skipped <- all_csv[!keep]
if (length(skipped)) {
  message("Note: ", length(skipped), " file(s) skipped due to unexpected directory layout or region not in {",
          paste(regions, collapse=","), "}.")
  # Optional: preview a few paths for debugging
  preview_n <- min(5L, length(skipped))
  if (preview_n > 0) {
    message("Examples of skipped paths:")
    for (i in seq_len(preview_n)) message("  - ", skipped[i])
  }
}

files <- all_csv[keep]
if (length(files) == 0L) {
  message("No matching CSVs found under ", base_dir, ". Exiting.")
  quit("no")
}
message("Found ", length(files), " CSV files to scan.")
ensure_ggpubr()

# --------------------- Reader ---------------------------
read_one <- function(fp) {
  parts <- strsplit(fp, "/", fixed = TRUE)[[1]]
  region <- parts[length(parts)-2]
  sra    <- parts[length(parts)-1]
  biome  <- biome_of(region)

  # Read only listed columns
  if (!check_csv_headers(fp)) return(NULL)
  df <- tryCatch({
    readr::read_csv(
      fp,
      col_types = cols(
        tar                 = col_skip(),
        member              = col_skip(),
        tree_index          = col_skip(),
        n_gene              = col_double(),
        n_imgvr             = col_double(),
        genes_monophyletic  = col_logical(),
        imgvr_monophyletic  = col_logical(),
        mrca_distance       = col_double(),
        notes               = col_character(),
        n_imgvr_phage       = col_double(),
        n_imgvr_non_phage   = col_double(),
        n_imgvr_unknown     = col_double(),
        imgvr_tree_flag     = col_character(),
        gene_pd             = col_double(),
        gene_mpd            = col_double(),
        gene_pairs          = col_double(),
        imgvr_pd            = col_double(),
        imgvr_mpd           = col_double(),
        imgvr_pairs         = col_double()
      ),
      progress = FALSE
    )
  }, error = function(e) {
    warning("Failed to read: ", fp, " - ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Keep only trees that are monophyletic for both
  df <- df %>%
    mutate(
      genes_monophyletic  = coalesce(genes_monophyletic,  FALSE),
      imgvr_monophyletic  = coalesce(imgvr_monophyletic,  FALSE)
    ) %>%
    filter(genes_monophyletic & imgvr_monophyletic)

  if (nrow(df) == 0) return(NULL)

  # Derive cluster from imgvr_tree_flag
  df <- df %>%
    mutate(
      flag = tolower(imgvr_tree_flag),
      flag = case_when(
        flag %in% c("non-phage","nonphage","non phage") ~ "non_phage",
        TRUE ~ flag
      ),
      cluster = case_when(
        flag %in% c("phage","non_phage","mixed") ~ flag,
        TRUE ~ "unknown"
      )
    )

  # Inclusion rules for MPDs
  df <- df %>%
    mutate(
      # metric-specific validity
      ok_gene_mpd  = !is.na(gene_mpd)  & (coalesce(gene_pairs,0)  >= 1 | coalesce(n_gene,0)  >= 2),
      ok_imgvr_mpd = !is.na(imgvr_mpd) & (coalesce(imgvr_pairs,0) >= 1 | coalesce(n_imgvr,0) >= 2),
      ok_bridging  = is.finite(mrca_distance)
    )

  # Apply minimum tip cutoffs (defaults set via env vars above)
  if (min_gene_tips > 0) {
    df <- df %>% mutate(ok_gene_mpd = ok_gene_mpd & coalesce(n_gene, 0) >= min_gene_tips)
  }
  if (min_imgvr_tips > 0) {
    df <- df %>% mutate(ok_imgvr_mpd = ok_imgvr_mpd & coalesce(n_imgvr, 0) >= min_imgvr_tips)
  }

  # Require that **all three** metrics are valid for a tree to be included at all
  df <- df %>% mutate(valid_all = ok_bridging & ok_gene_mpd & ok_imgvr_mpd) %>% filter(valid_all)

  # Long rows with three metrics (intersection-only)
  long <- bind_rows(
    df %>% transmute(biome = biome, region = region, sra = sra, cluster = cluster,
                     metric = "mrca_distance", value = mrca_distance),
    df %>% transmute(biome = biome, region = region, sra = sra, cluster = cluster,
                     metric = "gene_mpd", value = gene_mpd),
    df %>% transmute(biome = biome, region = region, sra = sra, cluster = cluster,
                     metric = "imgvr_mpd", value = imgvr_mpd)
  )

  long
}

# -------------------- Ingest all ------------------------
chunks <- lapply(files, read_one)
long_all <- bind_rows(chunks)

if (nrow(long_all) == 0L) {
  message("No data survived filters. Exiting.")
  quit("no")
}

# Coerce factors
long_all <- long_all %>%
  mutate(
    metric  = factor(metric, levels = c("mrca_distance","gene_mpd","imgvr_mpd"),
                     labels = c("bridging","mean_gene","mean_imgvr")),
    biome   = factor(biome, levels = c("marine","terrestrial")),
    cluster = factor(cluster, levels = c("phage","non_phage","mixed","unknown"))
  )

# Diagnostics
message("Row counts by cluster after monophyly filter:")
print(table(long_all$cluster, useNA = "ifany"))

# ------------------ Write CSV outputs -------------------
readr::write_csv(long_all, out_long)
message("Wrote: ", out_long, "  rows: ", nrow(long_all))

summarize_num <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(tibble(n = 0, mean = NA_real_, sd = NA_real_, min = NA_real_,
                                    q25 = NA_real_, median = NA_real_, q75 = NA_real_, max = NA_real_))
  tibble(
    n      = length(x),
    mean   = mean(x),
    sd     = ifelse(length(x) > 1, sd(x), NA_real_),
    min    = min(x),
    q25    = quantile(x, 0.25, names = FALSE),
    median = median(x),
    q75    = quantile(x, 0.75, names = FALSE),
    max    = max(x)
  )
}

stats <- long_all %>%
  group_by(biome, region, cluster, metric) %>%
  summarize(summarize_num(value), .groups = "drop")

readr::write_csv(stats, out_stats)
message("Wrote: ", out_stats)

# ----------------- Dunn helper for letters ----------------
get_multcomp_letters <- function(dunn_out, desired_levels) {
  comps <- dunn_out$comparisons
  if (is.null(comps) || length(comps) == 0) {
    return(data.frame(combo = desired_levels, groups = "a", stringsAsFactors = FALSE))
  }
  g1 <- sapply(comps, function(x) strsplit(x, " - ")[[1]][1])
  g2 <- sapply(comps, function(x) strsplit(x, " - ")[[1]][2])
  pvals <- dunn_out$altP  # BH-adjusted p-values

  all_groups <- unique(c(g1, g2))
  pMat <- matrix(1, nrow = length(all_groups), ncol = length(all_groups),
                 dimnames = list(all_groups, all_groups))
  for (i in seq_along(pvals)) {
    a <- g1[i]; b <- g2[i]
    pMat[a, b] <- pvals[i]
    pMat[b, a] <- pvals[i]
  }

  present <- intersect(desired_levels, rownames(pMat))
  pMat <- pMat[present, present, drop = FALSE]
  letters_vec <- multcompView::multcompLetters(pMat, threshold = 0.01, Letters = letters)$Letters
  data.frame(combo = names(letters_vec), groups = unname(letters_vec), stringsAsFactors = FALSE)
}

# ---------------------- Plot helpers ---------------------
# Colorblind-friendly (Okabe–Ito) palette for fills
# https://jfly.uni-koeln.de/color/
# NOTE: colors are reused for both fill (box) and color (jitter) to keep a clean legend-free aesthetic.
metric_fill_vals <- c(
  "bridging"   = "#0072B2", # blue
  "mean_gene"  = "#E69F00", # orange
  "mean_imgvr" = "#CC79A7"  # reddish purple
)

# Title text with KW p-value and n
page_title <- function(region = NULL, biome = NULL, cluster = NULL, n = NULL, kw_p = NULL, mode = "region") {
  if (mode == "region") {
    paste0("Region: ", region,
           "  |  Cluster: ", ifelse(cluster == "phage","Phage","Non-Phage"),
           "  |  (n=", n, ")",
           if (!is.null(kw_p)) paste0("  |  KW p=", signif(kw_p, 3)) else "")
  } else {
    paste0("Biome: ", biome,
           "  |  Cluster: ", ifelse(cluster == "phage","Phage","Non-Phage"),
           "  |  9-combo KW/Dunn across 3 regions × 3 metrics",
           if (!is.null(kw_p)) paste0("  |  KW p=", signif(kw_p, 3)) else "")
  }
}

# Per-region x-labels with per-metric n (short display names for single-region plots)
make_metric_labels <- function(dfm) {
  counts <- dfm %>% count(metric)
  lvls <- levels(dfm$metric)
  # Desired display names for single‑region plots
  disp_map <- c(
    "mean_gene"  = "Seqs",
    "bridging"   = "Bridge",
    "mean_imgvr" = "Refs"
  )
  labs <- setNames(lvls, lvls)
  for (m in lvls) {
    n_here <- counts$n[counts$metric == m]
    n_here <- ifelse(length(n_here) == 0, 0L, n_here)
    # Use mapped display name if present, otherwise fall back to the raw level
    pretty <- if (!is.null(disp_map[[m]])) disp_map[[m]] else m
    labs[m] <- paste0(pretty, " (n=", n_here, ")")
  }
  labs
}

# Precompute boxplot stats for a value column
# Returns: group, ymin, lower, middle, upper, ymax, n
box_stats <- function(df, group_col, value_col, positive_only = FALSE) {
  d <- df
  if (positive_only) d <- d %>% filter(.data[[value_col]] > 0)
  d %>%
    group_by(across(all_of(group_col))) %>%
    summarize(
      n      = n(),
      ymin   = min(.data[[value_col]], na.rm = TRUE),
      lower  = quantile(.data[[value_col]], 0.25, names = FALSE, na.rm = TRUE),
      middle = quantile(.data[[value_col]], 0.50, names = FALSE, na.rm = TRUE),
      upper  = quantile(.data[[value_col]], 0.75, names = FALSE, na.rm = TRUE),
      ymax   = max(.data[[value_col]], na.rm = TRUE),
      .groups = "drop"
    )
}

# ---------------- Jitter helpers (down-sampled raw points) ----------------
# Returns a tibble with columns: group_col, value
# - max_total controls the total number of points across all groups
# - positive_only filters to value > 0 (useful for log scale)
# - seed ensures reproducibility of the subsample
jitter_points <- function(df, group_col, value_col, max_total = 4000, positive_only = FALSE, seed = 123) {
  d <- df
  if (positive_only) d <- dplyr::filter(d, .data[[value_col]] > 0)
  if (!nrow(d)) {
    return(
      tibble::tibble(
        !!group_col := factor(character(), levels = if (!is.null(df[[group_col]]) && is.factor(df[[group_col]])) levels(df[[group_col]]) else character()),
        !!value_col := numeric()
      )
    )
  }

  # Ensure grouping column exists and is a factor with present levels
  if (is.null(d[[group_col]])) {
    stop("jitter_points(): grouping column '", group_col, "' not found in data.")
  }
  if (!is.factor(d[[group_col]])) d[[group_col]] <- factor(d[[group_col]])
  d[[group_col]] <- droplevels(d[[group_col]])

  set.seed(seed)
  # fair budget per group
  groups <- levels(d[[group_col]])
  groups <- groups[groups %in% unique(as.character(d[[group_col]]))]
  k <- length(groups)
  budget <- max(1L, floor(max_total / max(1L, k)))

  # Split explicitly so slice_sample() gets a constant n per split
  groups_split <- split(d, d[[group_col]], drop = TRUE)
  sampled_list <- lapply(groups_split, function(g) {
    n_take <- min(budget, nrow(g))
    if (n_take <= 0) return(g[FALSE, , drop = FALSE])
    dplyr::slice_sample(g, n = n_take, replace = FALSE)
  })
  dplyr::bind_rows(sampled_list)
}

# ---------------- Per-region pages (linear + log + letters) ----------------
plot_region_one_cluster <- function(df_region_cluster, region_name, cluster_name) {
  if (nrow(df_region_cluster) == 0) return(NULL)

  # KW across the three metrics (FULL data)
  df_rc <- df_region_cluster
  # Enforce requested metric order for single‑region plots: Gene, Bridging, IMGVR
  df_rc$metric <- factor(df_rc$metric, levels = c("mean_gene", "bridging", "mean_imgvr"))
  kw_p <- tryCatch(kruskal.test(value ~ metric, data = df_rc)$p.value, error = function(e) NA_real_)

  # Dunn letters if KW significant and at least 2 groups with data
  lvls <- levels(df_rc$metric)
  letters_df <- if (!is.na(kw_p) && kw_p < 0.01 && length(unique(df_rc$metric)) > 1) {
    dres <- dunn.test(x = df_rc$value, g = df_rc$metric, method = "bh", list = FALSE, altp = TRUE)
    get_multcomp_letters(dres, lvls)
  } else {
    data.frame(combo = lvls, groups = "a", stringsAsFactors = FALSE)
  }

  # Precomputed box stats for plotting (fast)
  bs_lin <- box_stats(df_rc, "metric", "value", positive_only = FALSE)
  bs_log <- box_stats(df_rc, "metric", "value", positive_only = TRUE)

  x_labels <- make_metric_labels(df_rc)
  total_n  <- nrow(df_rc)

  # Down-sampled jitter points for overlay (keeps PDFs light)
  jit_lin <- jitter_points(df_rc, "metric", "value",
                           max_total = JITTER_BUDGET_REGION,
                           positive_only = FALSE, seed = JITTER_SEED_REGION)
  jit_log <- jitter_points(df_rc, "metric", "value",
                           max_total = JITTER_BUDGET_REGION,
                           positive_only = TRUE,  seed = JITTER_SEED_REGION)

  # Cache data needed to reconstruct per-region plots later
  plot_cache$region[[paste(region_name, cluster_name, sep = "|")]] <- list(
    region    = region_name,
    cluster   = cluster_name,
    kw_p      = kw_p,
    total_n   = total_n,
    box_lin   = bs_lin,
    box_log   = bs_log,
    letters   = data.frame(metric = letters_df$combo, group = letters_df$groups, stringsAsFactors = FALSE),
    x_labels  = x_labels,
    jitter_lin = jit_lin,
    jitter_log = jit_log
  )

  # Letter positions
  letters_lin <- letters_df %>% rename(metric = combo) %>%
    left_join(bs_lin %>% transmute(metric, y = ymax*1.05), by = "metric")
  letters_log <- letters_df %>% rename(metric = combo) %>%
    left_join(bs_log %>% transmute(metric, y = ymax*1.05), by = "metric")

  p_lin <- ggplot(bs_lin,
                  aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity") +
    geom_point(data = jit_lin, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.15, height = 0, seed = JITTER_SEED_REGION),
               inherit.aes = FALSE, alpha = 0.35, size = 0.6, stroke = 0) +
    scale_fill_manual(values = metric_fill_vals, guide = "none") +
    scale_color_manual(values = metric_fill_vals, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6)
    ) +
    labs(title = "Linear scale", x = "Metric", y = "Patristic distance (substitutions/site)") +
    geom_text(data = letters_lin, aes(x = metric, y = y, label = groups),
              inherit.aes = FALSE, vjust = 0, size = 5)

  p_log <- ggplot(bs_log,
                  aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity") +
    geom_point(data = jit_log, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.15, height = 0, seed = JITTER_SEED_REGION),
               inherit.aes = FALSE, alpha = 0.35, size = 0.6, stroke = 0) +
    scale_y_log10() +
    scale_fill_manual(values = metric_fill_vals, guide = "none") +
    scale_color_manual(values = metric_fill_vals, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6)
    ) +
    labs(title = "Log10 scale", x = "Metric", y = "Patristic distance (substitutions/site)") +
    geom_text(data = letters_log, aes(x = metric, y = y, label = groups),
              inherit.aes = FALSE, vjust = 0, size = 5)

  lin_annot <- annotate_figure(
    p_lin,
    top = text_grob(page_title(region = region_name,
                               cluster = cluster_name,
                               n = total_n, kw_p = kw_p,
                               mode = "region"),
                    face = "bold", size = 14)
  )
  log_annot <- annotate_figure(
    p_log,
    top = text_grob(page_title(region = region_name,
                               cluster = cluster_name,
                               n = total_n, kw_p = kw_p,
                               mode = "region"),
                    face = "bold", size = 14)
  )
  list(lin_annot, log_annot)
}

# -------------- 9-combo pages per biome × cluster (letters across 9) ----------
plot_biome_nine_combo <- function(df_biome_cluster, biome_name, cluster_name, reg_levels) {
  if (nrow(df_biome_cluster) == 0) return(NULL)

  # Enforce region and metric order explicitly: North (left) — Equator (middle) — South (right)
  df_b <- df_biome_cluster %>%
    mutate(
      region = factor(region, levels = reg_levels, ordered = TRUE),
      metric = factor(metric, levels = c("mean_gene","bridging","mean_imgvr"), ordered = TRUE)
    ) %>%
    mutate(region_metric = paste(region, metric, sep = "|"))

  nine_levels <- as.vector(outer(reg_levels, levels(df_b$metric), paste, sep = "|"))
  df_b$region_metric <- factor(df_b$region_metric, levels = nine_levels, ordered = TRUE)

  # KW across 9 combos (FULL data)
  kw_p <- tryCatch(kruskal.test(value ~ region_metric, data = df_b)$p.value, error = function(e) NA_real_)

  letters_df <- if (!is.na(kw_p) && kw_p < 0.01 && length(unique(df_b$region_metric)) > 1) {
    dres <- dunn.test(x = df_b$value, g = df_b$region_metric, method = "bh", list = FALSE, altp = TRUE)
    get_multcomp_letters(dres, nine_levels)
  } else {
    data.frame(combo = nine_levels, groups = "a", stringsAsFactors = FALSE)
  }

  # Split combo into region and metric for caching
  letters_sep <- letters_df %>%
    tidyr::separate(combo, into = c("region","metric"), sep = "\\|", remove = TRUE)
  letters_sep <- letters_sep %>%
    mutate(
      region = factor(region, levels = reg_levels, ordered = TRUE),
      metric = factor(metric, levels = c("mean_gene","bridging","mean_imgvr"), ordered = TRUE)
    )

  # Box stats per facet for plotting
  bs <- df_b %>%
    group_by(region, metric) %>%
    summarize(
      n      = n(),
      ymin   = min(value, na.rm = TRUE),
      lower  = quantile(value, 0.25, names = FALSE, na.rm = TRUE),
      middle = quantile(value, 0.50, names = FALSE, na.rm = TRUE),
      upper  = quantile(value, 0.75, names = FALSE, na.rm = TRUE),
      ymax   = max(value, na.rm = TRUE),
      .groups = "drop"
    )
  bs <- bs %>%
    mutate(
      region = factor(region, levels = reg_levels, ordered = TRUE),
      metric = factor(metric, levels = c("mean_gene","bridging","mean_imgvr"), ordered = TRUE)
    )

  # Jitter per facet (keep light)
  # We'll sample a total of ~4500 points spread across 9 combos (~500 each)
  jit <- df_b %>%
    dplyr::group_by(region, metric) %>%
    dplyr::group_split(.keep = TRUE) %>%
    lapply(function(g) {
      jp <- jitter_points(g, group_col = "metric", value_col = "value",
                          max_total = JITTER_BUDGET_NINE, positive_only = FALSE, seed = JITTER_SEED_NINE)
      # Normalize metric to character to avoid incompatible ordered factor types across splits
      if (nrow(jp)) {
        jp$metric <- as.character(jp$metric)
      } else {
        jp$metric <- character()
      }
      jp
    }) %>%
    dplyr::bind_rows()
  if (nrow(jit)) {
    jit <- jit %>%
      mutate(
        region = factor(region, levels = reg_levels, ordered = TRUE),
        metric = factor(metric, levels = c("mean_gene","bridging","mean_imgvr"), ordered = TRUE)
      )
  }

  # Cache data needed to reconstruct 9-combo biome×cluster plots later
  plot_cache$nine[[paste(biome_name, cluster_name, sep = "|")]] <- list(
    biome      = biome_name,
    cluster    = cluster_name,
    reg_levels = reg_levels,
    kw_p       = kw_p,
    box        = bs,
    letters    = letters_sep %>%
                   dplyr::mutate(
                     region = factor(region, levels = reg_levels, ordered = TRUE),
                     metric = factor(metric, levels = c("mean_gene","bridging","mean_imgvr"), ordered = TRUE)
                   ) %>%
                   dplyr::select(region, metric, groups),
    jitter     = jit
  )

  # Add region_labels for nice facet titles
  region_labels <- if (biome_name == "terrestrial") {
    c(
      "terrestrial_north" = "Terrestrial North",
      "terrestrial_eq"    = "Terrestrial Equator",
      "terrestrial_south" = "Terrestrial South"
    )
  } else {
    c(
      "marine_north" = "Marine North",
      "marine_eq"    = "Marine Equator",
      "marine_south" = "Marine South"
    )
  }

  y_max <- max(bs$ymax, na.rm = TRUE) * 1.05
  letters_df <- letters_df %>%
    tidyr::separate(combo, into = c("region","metric"), sep = "\\|", remove = FALSE) %>%
    mutate(
      region = factor(region, levels = reg_levels, ordered = TRUE),
      metric = factor(metric, levels = c("mean_gene","bridging","mean_imgvr"), ordered = TRUE),
      y = y_max
    )

  p <- ggplot(bs,
              aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity") +
    geom_point(data = jit, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.15, height = 0, seed = JITTER_SEED_NINE),
               inherit.aes = FALSE, alpha = 0.30, size = 0.5, stroke = 0) +
    facet_wrap(~ region, scales = "fixed", drop = FALSE, labeller = labeller(region = region_labels)) +
    scale_fill_manual(values = metric_fill_vals, guide = "none") +
    scale_color_manual(values = metric_fill_vals, guide = "none") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6)
    ) +
    labs(x = "Metric", y = "Patristic distance (substitutions/site)") +
    geom_text(data = letters_df, aes(x = metric, y = y, label = groups),
              inherit.aes = FALSE, vjust = 0, size = 5)

  annotate_figure(p, top = text_grob(page_title(biome = biome_name, cluster = cluster_name, kw_p = kw_p, mode = "biome"),
                                     face = "bold", size = 14))
}

# --------------------- Re-plot from cache helpers ----------------------
plot_from_region_cache <- function(ent) {
  # ent: one element from plot_cache$region[[...]]
  bs_lin <- ent$box_lin; bs_log <- ent$box_log
  x_labels <- ent$x_labels
  letters <- ent$letters
  letters_lin <- letters %>% left_join(bs_lin %>% transmute(metric, y = ymax*1.05), by = "metric")
  letters_log <- letters %>% left_join(bs_log %>% transmute(metric, y = ymax*1.05), by = "metric")

  p_lin <- ggplot(bs_lin,
                  aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity") +
    geom_point(data = ent$jitter_lin, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.15, height = 0, seed = JITTER_SEED_REGION),
               inherit.aes = FALSE, alpha = 0.35, size = 0.6, stroke = 0) +
    scale_fill_manual(values = metric_fill_vals, guide = "none") +
    scale_color_manual(values = metric_fill_vals, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6)
    ) +
    labs(title = page_title(region = ent$region, cluster = ent$cluster, n = ent$total_n, kw_p = ent$kw_p, mode = "region"),
         x = "Metric", y = "Patristic distance (substitutions/site)") +
    geom_text(data = letters_lin, aes(x = metric, y = y, label = group),
              inherit.aes = FALSE, vjust = 0, size = 5)

  p_log <- ggplot(bs_log,
                  aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity") +
    geom_point(data = ent$jitter_log, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.15, height = 0, seed = JITTER_SEED_REGION),
               inherit.aes = FALSE, alpha = 0.35, size = 0.6, stroke = 0) +
    scale_y_log10() +
    scale_fill_manual(values = metric_fill_vals, guide = "none") +
    scale_color_manual(values = metric_fill_vals, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6)
    ) +
    labs(title = page_title(region = ent$region, cluster = ent$cluster, n = ent$total_n, kw_p = ent$kw_p, mode = "region"),
         x = "Metric", y = "Patristic distance (substitutions/site)") +
    geom_text(data = letters_log, aes(x = metric, y = y, label = group),
              inherit.aes = FALSE, vjust = 0, size = 5)

  list(p_lin, p_log)
}

plot_from_nine_cache <- function(ent) {
  bs <- ent$box; letters <- ent$letters; jit <- ent$jitter
  biome_name <- ent$biome; cluster_name <- ent$cluster
  reg_levels <- ent$reg_levels

  region_labels <- if (biome_name == "terrestrial") {
    c("terrestrial_north" = "Terrestrial North",
      "terrestrial_eq"    = "Terrestrial Equator",
      "terrestrial_south" = "Terrestrial South")
  } else {
    c("marine_north" = "Marine North",
      "marine_eq"    = "Marine Equator",
      "marine_south" = "Marine South")
  }

  y_max <- max(bs$ymax, na.rm = TRUE) * 1.05
  letters_df <- letters %>% mutate(y = y_max)

  p <- ggplot(bs,
              aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity") +
    geom_point(data = jit, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.15, height = 0, seed = JITTER_SEED_NINE),
               inherit.aes = FALSE, alpha = 0.30, size = 0.5, stroke = 0) +
    facet_wrap(~ region, scales = "fixed", labeller = labeller(region = region_labels)) +
    scale_fill_manual(values = metric_fill_vals, guide = "none") +
    scale_color_manual(values = metric_fill_vals, guide = "none") +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "gray95", color = NA),
      plot.title = element_text(face = "bold"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.6)
    ) +
    labs(x = "Metric", y = "Patristic distance (substitutions/site)") +
    geom_text(data = letters_df, aes(x = metric, y = y, label = groups),
              inherit.aes = FALSE, vjust = 0, size = 5)

  annotate_figure(p, top = text_grob(page_title(biome = biome_name, cluster = cluster_name, kw_p = ent$kw_p, mode = "biome"),
                                     face = "bold", size = 14))
}

# ------------------------- PDF --------------------------
pdf(out_pdf, width = 11, height = 8)

# A) Per-region pages: each region × {phage, non_phage}
message("Starting per-region pages...")
for (r in regions) {
  for (cl in c("phage","non_phage")) {
    message("Plotting region=", r, "  cluster=", cl)
    df_rc <- long_all %>% filter(region == r, cluster == cl)
    pg <- plot_region_one_cluster(df_rc, region_name = r, cluster_name = cl)
    if (!is.null(pg)) for (p in pg) print(p)
  }
}

# B) 9-combo pages per biome × cluster
message("Starting 9-combo biome × cluster pages...")
for (bm in c("marine","terrestrial")) {
  reg_levels <- if (bm == "marine") marine_levels else terrestrial_levels
  for (cl in c("phage","non_phage")) {
    df_bc <- long_all %>% filter(biome == bm, cluster == cl)
    message("Plotting biome=", bm, "  cluster=", cl, "  n=", nrow(df_bc))
    pg <- plot_biome_nine_combo(df_bc, biome_name = bm, cluster_name = cl, reg_levels = reg_levels)
    if (!is.null(pg)) print(pg)
  }
}

dev.off()
message("Wrote: ", out_pdf)

# Save cached computations so plots can be rebuilt without recomputing stats
saveRDS(as.list(plot_cache), cache_rds)
message("Cached plot data to: ", cache_rds)

# ---- Human-readable CSV exports built from cache ----
out_region_boxes   <- file.path(base_dir, "boxstats_region.csv")
out_region_letters <- file.path(base_dir, "letters_region.csv")
out_nine_boxes     <- file.path(base_dir, "boxstats_nine.csv")
out_nine_letters   <- file.path(base_dir, "letters_nine.csv")

# Region-level box stats (linear + log10)
if (length(plot_cache$region)) {
  region_box_rows <- lapply(plot_cache$region, function(ent) {
    lin <- ent$box_lin %>% mutate(scale = "linear")
    lg  <- ent$box_log %>% mutate(scale = "log10")
    bind_rows(lin, lg) %>% mutate(region = ent$region, cluster = ent$cluster)
  })
  region_boxes <- bind_rows(region_box_rows) %>%
    select(region, cluster, scale, metric, n, ymin, lower, middle, upper, ymax)
  readr::write_csv(region_boxes, out_region_boxes)
  message("Wrote: ", out_region_boxes, "  rows: ", nrow(region_boxes))

  region_letters_rows <- lapply(plot_cache$region, function(ent) {
    ent$letters %>% transmute(region = ent$region, cluster = ent$cluster, metric = metric, group = group)
  })
  region_letters <- bind_rows(region_letters_rows)
  readr::write_csv(region_letters, out_region_letters)
  message("Wrote: ", out_region_letters, "  rows: ", nrow(region_letters))
}

if (length(plot_cache$nine)) {
  nine_box_rows <- lapply(plot_cache$nine, function(ent) {
    ent$box %>%
      dplyr::mutate(
        biome   = ent$biome,
        cluster = ent$cluster,
        # Normalize to character to avoid ordered-factor level mismatches across biomes
        region  = as.character(region),
        metric  = as.character(metric)
      )
  })
  nine_boxes <- dplyr::bind_rows(nine_box_rows) %>%
    dplyr::select(biome, region, cluster, metric, n, ymin, lower, middle, upper, ymax)
  readr::write_csv(nine_boxes, out_nine_boxes)
  message("Wrote: ", out_nine_boxes, "  rows: ", nrow(nine_boxes))

  nine_letters_rows <- lapply(plot_cache$nine, function(ent) {
    ent$letters %>%
      dplyr::mutate(
        biome   = ent$biome,
        cluster = ent$cluster,
        region  = as.character(region),
        metric  = as.character(metric)
      )
  })
  nine_letters <- dplyr::bind_rows(nine_letters_rows) %>%
    dplyr::select(biome, region, cluster, metric, groups)
  readr::write_csv(nine_letters, out_nine_letters)
  message("Wrote: ", out_nine_letters, "  rows: ", nrow(nine_letters))
}

message("Tip: You can rebuild figures later using plot_from_region_cache()/plot_from_nine_cache() with entries from readRDS(cache_rds)$region / $nine.")