#!/usr/bin/env Rscript

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
})

# ---- Inputs & Outputs (read-only cache; write-only figures) ----
BASE_DIR   <- Sys.getenv("BASE_DIR",   unset = setwd(getwd()))
CACHE_FILE <- Sys.getenv("CACHE_FILE", unset = "boxplot_cache.rds")
CACHE_PATH <- file.path(BASE_DIR, CACHE_FILE)

OUT_DIR    <- file.path(BASE_DIR, "rebuild_figures")
SVG_DIR    <- file.path(OUT_DIR, "svgs")
PDF_REGION_PATH <- file.path(OUT_DIR, "figures_region.pdf")
PDF_GLOBAL_PATH <- file.path(OUT_DIR, "figures_global.pdf")
OUT_SVG    <- tolower(Sys.getenv("OUT_SVG", unset = "1")) %in% c("1","true","yes")

# ---- Locked-in styling ----
THEME_BASE_SIZE <- 13
AXIS_X_LABEL    <- "Metric"
AXIS_Y_LABEL    <- "Patristic distance (substitutions/site)"
TITLE_PREFIX    <- ""
LETTERS_Y_MULT <- 1.05

# Metric order and display labels (presentation)
METRIC_ORDER  <- c("mean_gene","bridging","mean_imgvr")
METRIC_LABELS <- c("Reconstructed","Bridge","IMG/VR")
metric_factor <- function(x) factor(x, levels = METRIC_ORDER, ordered = TRUE)
metric_label_map <- setNames(METRIC_LABELS, METRIC_ORDER)

# Region orders (for facet order)
REGION_ORDER_MARINE <- c("marine_north","marine_eq","marine_south")
REGION_ORDER_TERR   <- c("terrestrial_north","terrestrial_eq","terrestrial_south")

# Clusters: include both by default
FILTER_CLUSTER <- ""

# Scales (log for region plots only)
USE_LOG_Y_REGION <- TRUE
USE_LOG_Y_BIOME9 <- FALSE

# Jitter defaults
JITTER_ALPHA <- 0.30
JITTER_SIZE  <- 0.5

# Minimal, colorblind-safe palette (Okabe–Ito)
metric_colors <- c(
  "bridging"   = "#0072B2",
  "mean_gene"  = "#E69F00",
  "mean_imgvr" = "#CC79A7"
)

# ---- Load cache (read-only) ----
if (!file.exists(CACHE_PATH)) stop("Cache not found: ", CACHE_PATH)
cache <- readRDS(CACHE_PATH)

# ---- Helpers ----
to_char_box <- function(df) {
  nm <- names(df)
  if ("region" %in% nm) df$region <- as.character(df$region)
  if ("metric" %in% nm) df$metric <- as.character(df$metric)
  df
}
to_char_letters <- function(df) {
  nm <- names(df)
  if ("region" %in% nm) df$region <- as.character(df$region)
  if ("metric" %in% nm) df$metric <- as.character(df$metric)
  df
}
metric_pretty_labels <- function(lvls) {
  labs <- setNames(lvls, lvls)
  for (m in lvls) if (!is.null(metric_label_map[[m]])) labs[m] <- metric_label_map[[m]]
  labs
}
apply_region_levels <- function(x, biome) {
  if (biome == "marine") factor(x, levels = REGION_ORDER_MARINE, ordered = TRUE)
  else                   factor(x, levels = REGION_ORDER_TERR,   ordered = TRUE)
}

# ---- Plotters ----
plot_region <- function(ent) {
  bs_lin <- to_char_box(ent$box_lin)
  bs_log <- to_char_box(ent$box_log)
  bs_lin$metric <- metric_factor(bs_lin$metric)
  bs_log$metric <- metric_factor(bs_log$metric)
  labs_x <- metric_pretty_labels(levels(bs_lin$metric))
  # Ensure no '(n=...)' appears in x-axis tick labels
  labs_x <- setNames(gsub("\\s*\\(n\\s*=.*?\\)\\s*", "", unname(labs_x)), names(labs_x))

  # Significance letters (use cached table; place at top of each box scale)
  letters_df <- to_char_letters(ent$letters)
  letters_lin <- letters_df %>% dplyr::left_join(bs_lin %>% dplyr::transmute(metric, y = ymax*LETTERS_Y_MULT), by = "metric")
  letters_log <- letters_df %>% dplyr::left_join(bs_log %>% dplyr::transmute(metric, y = ymax*LETTERS_Y_MULT), by = "metric")

  # Derive a single per-region tree count from box stats (same across metrics)
  trees_n <- suppressWarnings(as.integer(stats::median(bs_lin$n, na.rm = TRUE)))
  if (is.na(trees_n)) trees_n <- 0L

  ttl <- paste0(if (nzchar(TITLE_PREFIX)) paste0(TITLE_PREFIX, " ") else "",
                "Region: ", ent$region, " | Cluster: ", ent$cluster, " | n=", trees_n)

  p_lin <- ggplot(bs_lin,
                  aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity", alpha = 0.5, color = "black") +
    geom_point(data = ent$jitter_lin, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.2, height = 0),
               inherit.aes = FALSE, alpha = 0.6, size = 0.8, stroke = 0) +
    scale_fill_manual(values = metric_colors, guide = "none") +
    scale_color_manual(values = metric_colors, guide = "none") +
    scale_x_discrete(labels = labs_x) +
    geom_text(data = letters_lin, aes(x = metric, y = y, label = group),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    theme_bw(base_size = THEME_BASE_SIZE) +
    labs(title = ttl, x = AXIS_X_LABEL, y = AXIS_Y_LABEL)

  p_log <- ggplot(bs_log,
                  aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity", alpha = 0.5, color = "black") +
    geom_point(data = ent$jitter_log, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.2, height = 0),
               inherit.aes = FALSE, alpha = 0.6, size = 0.8, stroke = 0) +
    scale_fill_manual(values = metric_colors, guide = "none") +
    scale_color_manual(values = metric_colors, guide = "none") +
    scale_x_discrete(labels = labs_x) +
    geom_text(data = letters_log, aes(x = metric, y = y, label = group),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    theme_bw(base_size = THEME_BASE_SIZE) +
    labs(title = ttl, x = AXIS_X_LABEL, y = AXIS_Y_LABEL)

  if (USE_LOG_Y_REGION) p_log <- p_log + scale_y_log10()

  list(p_lin, p_log)
}

plot_biome9 <- function(ent, use_log = FALSE) {
  bs <- to_char_box(ent$box)
  jit <- ent$jitter
  biome_name <- ent$biome; cluster_name <- ent$cluster

  bs$region <- apply_region_levels(bs$region, biome_name)
  bs$metric <- metric_factor(bs$metric)

  # X‑axis labels for metrics (pretty + sanitized)
  labs_x <- metric_pretty_labels(levels(bs$metric))
  labs_x <- setNames(gsub("\n*\\s*\\(n\\s*=.*?\\)\\s*", "", unname(labs_x)), names(labs_x))

  # Compute per-region n (expect identical across metrics within a region; if not, use max)
  n_split <- split(as.integer(bs$n), bs$region)
  n_by_region <- vapply(n_split, function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0L) 0L else if (length(unique(v)) == 1L) unique(v) else max(v, na.rm = TRUE)
  }, integer(1))
  trees_n_biome <- sum(n_by_region, na.rm = TRUE)

  # Significance letters for 3×3; cache can use column name 'groups' or 'group'
  letters_raw <- to_char_letters(ent$letters)
  letters_raw$label <- if ("groups" %in% names(letters_raw)) letters_raw$groups else if ("group" %in% names(letters_raw)) letters_raw$group else ""

  y_max <- max(bs$ymax, na.rm = TRUE) * LETTERS_Y_MULT
  letters_9 <- letters_raw %>% dplyr::mutate(y = y_max,
                                            region = apply_region_levels(region, biome_name),
                                            metric = metric_factor(metric))

  # Base pretty names per biome
  base_labels <- if (biome_name == "terrestrial") {
    c("terrestrial_north" = "Terrestrial North",
      "terrestrial_eq"    = "Terrestrial Equator",
      "terrestrial_south" = "Terrestrial South")
  } else {
    c("marine_north" = "Marine North",
      "marine_eq"    = "Marine Equator",
      "marine_south" = "Marine South")
  }
  # Append (n=...) per region; fill missing with 0 if needed
  all_keys <- names(base_labels)
  n_by_region <- n_by_region[all_keys]
  n_by_region[is.na(n_by_region)] <- 0L
  region_labels <- setNames(paste0(base_labels, " (n=", n_by_region, ")"), all_keys)

  ttl <- paste0(if (nzchar(TITLE_PREFIX)) paste0(TITLE_PREFIX, " ") else "",
                "Biome: ", biome_name, " | Cluster: ", ifelse(cluster_name == "phage","Phage","Non-Phage"),
                " | n=", trees_n_biome)

  p <- ggplot(bs,
              aes(x = metric, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = metric)) +
    geom_boxplot(stat = "identity", alpha = 0.5, color = "black") +
    geom_point(data = jit, aes(x = metric, y = value, color = metric),
               position = position_jitter(width = 0.2, height = 0),
               inherit.aes = FALSE, alpha = 0.6, size = 0.8, stroke = 0) +
    facet_wrap(~ region, scales = "fixed", labeller = labeller(region = region_labels)) +
    scale_x_discrete(labels = labs_x) +
    geom_text(data = letters_9, aes(x = metric, y = y, label = label),
              inherit.aes = FALSE, vjust = 0, size = 5) +
    scale_fill_manual(values = metric_colors, guide = "none") +
    scale_color_manual(values = metric_colors, guide = "none") +
    theme_bw(base_size = THEME_BASE_SIZE) +
    labs(title = ttl, x = AXIS_X_LABEL, y = AXIS_Y_LABEL)

  if (use_log) p <- p + scale_y_log10()

  full_title <- paste0("Global Virome Analysis (KW test with BH) - Biome: ", biome_name,
                       " | Cluster: ", ifelse(cluster_name == "phage","Phage","Non-Phage"),
                       " | total n=", trees_n_biome)

  p + ggtitle(full_title) + theme(plot.title = element_text(hjust = 0.5, face = "plain"))
}

# ---- Build figures (vector only) ----
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (OUT_SVG) dir.create(SVG_DIR, recursive = TRUE, showWarnings = FALSE)

region_entries <- cache$region
nine_entries   <- cache$nine
if (nzchar(FILTER_CLUSTER)) {
  region_entries <- Filter(function(ent) ent$cluster == FILTER_CLUSTER, region_entries)
  nine_entries   <- Filter(function(ent) ent$cluster == FILTER_CLUSTER, nine_entries)
}

# PDF (vector) - Per‑region only
pdf(PDF_REGION_PATH, width = 11, height = 8)

# Per‑region
for (k in names(region_entries)) {
  ent <- region_entries[[k]]
  pg <- plot_region(ent)
  for (pp in pg) print(pp)
  if (OUT_SVG) {
    svglite::svglite(file.path(SVG_DIR, paste0("region_", ent$region, "_", ent$cluster, "_linear.svg")), width = 11, height = 8)
    print(pg[[1]]); dev.off()
    svglite::svglite(file.path(SVG_DIR, paste0("region_", ent$region, "_", ent$cluster, "_log.svg")), width = 11, height = 8)
    print(pg[[2]]); dev.off()
  }
}

dev.off()

# PDF (vector) - Global 3×3 pages
pdf(PDF_GLOBAL_PATH, width = 11, height = 8)

# Biome × cluster (ordered: Phage first, then Non‑Phage), and both linear & log pages
biome_order <- c("marine", "terrestrial")
cluster_order <- c("phage", "non_phage")
for (bm in biome_order) {
  for (cl in cluster_order) {
    entries <- Filter(function(e) e$biome == bm && e$cluster == cl, nine_entries)
    if (!length(entries)) next
    ent <- entries[[1]]

    # Linear
    p9_lin <- plot_biome9(ent, use_log = FALSE)
    print(p9_lin)
    if (OUT_SVG) {
      svglite::svglite(file.path(SVG_DIR, paste0("biome_", ent$biome, "_", ent$cluster, "_nine_linear.svg")), width = 11, height = 8)
      print(p9_lin); dev.off()
    }

    # Log
    p9_log <- plot_biome9(ent, use_log = TRUE)
    print(p9_log)
    if (OUT_SVG) {
      svglite::svglite(file.path(SVG_DIR, paste0("biome_", ent$biome, "_", ent$cluster, "_nine_log.svg")), width = 11, height = 8)
      print(p9_log); dev.off()
    }
  }
}

dev.off()
message("Figures written:\n  PDF (per‑region): ", PDF_REGION_PATH,
        "\n  PDF (global 3×3): ", PDF_GLOBAL_PATH,
        if (OUT_SVG) paste0("\n  SVGs: ", SVG_DIR) else "")
