# =============================================================================
# Wire Myography Analysis
# =============================================================================
# Author : Agnieszka Karas, PhD
# Contact: agaakaras@gmail.com | linkedin.com/in/agnieszka-karas
# =============================================================================

library(dplyr)
library(ggplot2)
library(writexl)
library(scales)

# =============================================================================
#Load file: insert full path of the txt file for analysis and desired output directory
# =============================================================================

INPUT_FILE  <- "/../sample_data.txt"
OUTPUT_DIR  <- "/../output"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
df <- read.table(
  INPUT_FILE, sep = "\t", dec = ",", skip = 9, header = FALSE,
  col.names = c("t", 1:8, "comment"),
  fill = TRUE, quote = "", comment.char = "", na.strings = ""
)
cat(sprintf("Loaded %d rows | %.0f s\n", nrow(df), max(df$t) - min(df$t)))

df$comment <- trimws(as.character(df$comment))
df$comment[is.na(df$comment)] <- ""

preview_comments <- function(df) {
  df[grepl("^#\\*", df$comment), c("t", "comment")]
}
cat("\nProtocol comments:\n")
print(preview_comments(df))
# =============================================================================
# Settings — edit here Phe mode and comment labels
# =============================================================================

# Phe protocol: 'dose_response' or 'single_dose'
PHE_MODE <- "dose_response"

# Cumulative dose labels — must match LabChart annotations exactly
PHE_LABELS <- c("Phe0,01", "Phe0,03", "Phe0,1", "Phe0,3", "Phe1", "Phe3")
PHE_CONC   <- c(0.01, 0.03, 0.1, 0.3, 1.0, 3.0)   # uM

ACH_LABELS <- c("Ach0,001", "Ach0,01", "Ach0,1", "Ach1", "Ach10")
ACH_CONC   <- c(0.001, 0.01, 0.1, 1.0, 10.0)       # uM

SNP_LABELS <- c("SNP0,001", "SNP0,01", "SNP0,1", "SNP1")
SNP_CONC   <- c(0.001, 0.01, 0.1, 1.0)              # uM

# Fixed protocol labels — edit if your file uses different annotations
L_KCL60       <- "KCl60"    # KCl 60 mM addition
L_KCL60_END   <- "P2"       # rinse after KCl60
L_PHE_END     <- "PP"       # plateau marker after Phe
L_PHE_SINGLE  <- "Phe3uM"  # single-dose Phe label (single_dose mode only)
L_SUBPHE      <- "subPhe"   # submaximal Phe before ACh
L_SUBPHE_END  <- "Ach0,001" # first ACh dose (= end of subPhe window)
L_SUBPHE2     <- "2subPhe"  # submaximal Phe before SNP
L_SUBPHE2_END <- "SNP0,001" # first SNP dose (= end of 2subPhe window)
L_ACH_END     <- "P3"       # rinse after ACh
L_SNP_END     <- "K"        # end of recording

CHANNEL_COLS <- paste0("X", 1:8)   

# =============================================================================
# Functions
# =============================================================================

find_row <- function(df, label) {
  stripped <- trimws(sub("^#\\*\\s*", "", df$comment))
  target   <- trimws(sub("^#\\*\\s*", "", label))
  idx      <- which(stripped == target)
  if (length(idx) == 0) {
    avail <- df$comment[grepl("^#\\*", df$comment)]
    stop(sprintf("Label '%s' not found.\nAvailable: %s",
                 label, paste(avail, collapse = ", ")))
  }
  idx[1]
}

compute_contraction <- function(df, start_label, end_label) {
  s        <- find_row(df, start_label)
  e        <- find_row(df, end_label)
  baseline <- as.numeric(df[s - 1, CHANNEL_COLS])
  max_vals <- apply(df[s:e, CHANNEL_COLS], 2, max, na.rm = TRUE)
  setNames(as.numeric(max_vals - baseline), CHANNEL_COLS)
}

compute_cumulative_contraction <- function(df, dose_labels, end_label) {
  idxs     <- sapply(dose_labels, find_row, df = df)
  baseline <- as.numeric(df[idxs[1] - 1, CHANNEL_COLS])
  result   <- lapply(seq_along(dose_labels), function(i) {
    s   <- idxs[i]
    e   <- if (i < length(dose_labels)) idxs[i + 1] else find_row(df, end_label)
    setNames(as.numeric(apply(df[s:e, CHANNEL_COLS], 2, max, na.rm = TRUE) - baseline),
             CHANNEL_COLS)
  })
  df_out <- as.data.frame(do.call(rbind, result))
  rownames(df_out) <- dose_labels
  df_out
}

compute_cumulative_relaxation <- function(df, dose_labels, end_label) {
  idxs     <- sapply(dose_labels, find_row, df = df)
  baseline <- as.numeric(df[idxs[1] - 1, CHANNEL_COLS])
  result   <- lapply(seq_along(dose_labels), function(i) {
    s   <- idxs[i]
    e   <- if (i < length(dose_labels)) idxs[i + 1] else find_row(df, end_label)
    setNames(as.numeric(baseline - apply(df[s:e, CHANNEL_COLS], 2, min, na.rm = TRUE)),
             CHANNEL_COLS)
  })
  df_out <- as.data.frame(do.call(rbind, result))
  rownames(df_out) <- dose_labels
  df_out
}

# =============================================================================
# Main analysis
# =============================================================================

analyse_experiment <- function(df, phe_mode, phe_labels, ach_labels, snp_labels,
                                l = list(
                                  KCl60       = L_KCL60,
                                  KCl60_end   = L_KCL60_END,
                                  phe_end     = L_PHE_END,
                                  phe_single  = L_PHE_SINGLE,
                                  subPhe      = L_SUBPHE,
                                  subPhe_end  = L_SUBPHE_END,
                                  subPhe2     = L_SUBPHE2,
                                  subPhe2_end = L_SUBPHE2_END,
                                  ach_end     = L_ACH_END,
                                  snp_end     = L_SNP_END
                                )) {

  make_row <- function(metric, vals) {
    row        <- as.data.frame(t(round(as.numeric(unname(vals)), 4)))
    names(row) <- paste0("ch", 1:8)
    cbind(data.frame(metric = metric, stringsAsFactors = FALSE), row)
  }

  rows_mn  <- list()
  rows_pct <- list()

  # KCl 60 mM
  kcl60 <- compute_contraction(df, l$KCl60, l$KCl60_end)
  rows_mn[["KCl60_mN"]] <- make_row("KCl60_mN", kcl60)

  # Phenylephrine
  if (phe_mode == "dose_response") {
    phe_abs <- compute_cumulative_contraction(df, phe_labels, l$phe_end)
    phe_pct <- sweep(phe_abs, 2, kcl60, "/") * 100
    for (label in phe_labels) {
      rows_mn [[paste0("Phe_", label, "_mN")]]    <- make_row(paste0("Phe_", label, "_mN"),    phe_abs[label, ])
      rows_pct[[paste0("Phe_", label, "_%_KCl")]] <- make_row(paste0("Phe_", label, "_%_KCl"), phe_pct[label, ])
    }
  } else {
    phe_max <- compute_contraction(df, l$phe_single, l$phe_end)
    rows_mn [["Phe_3uM_mN"]]    <- make_row("Phe_3uM_mN",    phe_max)
    rows_pct[["Phe_3uM_%_KCl"]] <- make_row("Phe_3uM_%_KCl", phe_max / kcl60 * 100)
    phe_pct <- NULL
  }

  # Submaximal Phe
  phe_sub_ach <- compute_contraction(df, l$subPhe,  l$subPhe_end)
  phe_sub_snp <- compute_contraction(df, l$subPhe2, l$subPhe2_end)
  rows_mn[["subPhe_ACh_mN"]] <- make_row("subPhe_ACh_mN", phe_sub_ach)
  rows_mn[["subPhe_SNP_mN"]] <- make_row("subPhe_SNP_mN", phe_sub_snp)

  # ACh relaxation
  ach_abs <- compute_cumulative_relaxation(df, ach_labels, l$ach_end)
  ach_pct <- sweep(ach_abs, 2, phe_sub_ach, "/") * 100
  for (label in ach_labels) {
    rows_mn [[paste0("ACh_", label, "_mN")]] <- make_row(paste0("ACh_", label, "_mN"), ach_abs[label, ])
    rows_pct[[paste0("ACh_", label, "_%")]]  <- make_row(paste0("ACh_", label, "_%"),  ach_pct[label, ])
  }

  # SNP relaxation
  snp_abs <- compute_cumulative_relaxation(df, snp_labels, l$snp_end)
  snp_pct <- sweep(snp_abs, 2, phe_sub_snp, "/") * 100
  for (label in snp_labels) {
    rows_mn [[paste0("SNP_", label, "_mN")]] <- make_row(paste0("SNP_", label, "_mN"), snp_abs[label, ])
    rows_pct[[paste0("SNP_", label, "_%")]]  <- make_row(paste0("SNP_", label, "_%"),  snp_pct[label, ])
  }

  # mN rows first, then % rows
  result_df <- do.call(rbind, c(rows_mn, rows_pct))
  rownames(result_df) <- NULL

  attr(result_df, "kcl60")   <- kcl60
  attr(result_df, "ach_pct") <- ach_pct
  attr(result_df, "snp_pct") <- snp_pct
  attr(result_df, "phe_pct") <- phe_pct

  result_df
}

# =============================================================================
# Plotting
# =============================================================================

CH_COLORS <- hue_pal()(8)

ylim_auto <- function(vals, margin = 0.05) {
  mn <- min(vals, na.rm = TRUE); mx <- max(vals, na.rm = TRUE)
  span <- if (mx != mn) mx - mn else max(abs(mx), 1)
  c(mn - span * margin, mx + span * margin)
}

to_long <- function(df_wide, conc_vec) {
  df_wide$dose <- conc_vec
  tidyr::pivot_longer(df_wide, cols = all_of(CHANNEL_COLS),
                      names_to = "channel", values_to = "value")
}

plot_summary <- function(results, phe_mode, phe_conc, ach_conc, snp_conc) {

  if (!requireNamespace("tidyr",     quietly = TRUE)) install.packages("tidyr")
  if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
  library(tidyr); library(patchwork)

  kcl60   <- attr(results, "kcl60")
  ach_pct <- attr(results, "ach_pct")
  snp_pct <- attr(results, "snp_pct")
  phe_pct <- attr(results, "phe_pct")

  theme_mg <- theme_classic(base_size = 11) +
    theme(plot.title      = element_text(size = 11, face = "bold"),
          legend.key.size = unit(0.4, "cm"),
          legend.text     = element_text(size = 8))

  # Phe
  if (!is.null(phe_pct)) {
    phe_long <- to_long(phe_pct, phe_conc)
    yl <- ylim_auto(phe_long$value)
    p_phe <- ggplot(phe_long, aes(dose, value, colour = channel, group = channel)) +
      geom_line() +
      geom_point(shape = 21, fill = "white", size = 2.5, stroke = 1.2) +
      scale_x_log10(labels = label_scientific()) +
      scale_colour_manual(values = setNames(CH_COLORS, paste0("X", 1:8)),
                          labels = paste0("Ch ", 1:8)) +
      coord_cartesian(ylim = yl) +
      labs(x = "Phenylephrine (uM)", y = "Contraction (% KCl 60 mM)",
           title = "Phenylephrine dose-response", colour = NULL) +
      theme_mg
  } else {
    bar_df <- data.frame(
      channel = paste0("Ch ", 1:8),
      value   = as.numeric(results[results$metric == "Phe_3uM_%_KCl",
                                    paste0("ch", 1:8)])
    )
    yl <- ylim_auto(bar_df$value)
    p_phe <- ggplot(bar_df, aes(channel, value, fill = channel)) +
      geom_col(width = 0.65, colour = "white") +
      geom_text(aes(label = round(value, 1)), vjust = -0.4, size = 3) +
      geom_hline(yintercept = mean(bar_df$value), linetype = "dashed") +
      scale_fill_manual(values = setNames(CH_COLORS, paste0("Ch ", 1:8))) +
      coord_cartesian(ylim = yl) +
      labs(x = NULL, y = "Contraction (% KCl 60 mM)",
           title = "Phenylephrine 3 uM — max contraction") +
      theme_mg + theme(legend.position = "none",
                       axis.text.x = element_text(angle = 30, hjust = 1))
  }

  # ACh — inverted Y
  ach_long <- to_long(ach_pct, ach_conc)
  yl_ach   <- ylim_auto(ach_long$value)
  p_ach <- ggplot(ach_long, aes(dose, value, colour = channel, group = channel)) +
    geom_line() +
    geom_point(shape = 21, fill = "white", size = 2.5, stroke = 1.2) +
    scale_x_log10(labels = label_scientific()) +
    scale_colour_manual(values = setNames(CH_COLORS, paste0("X", 1:8)),
                        labels = paste0("Ch ", 1:8)) +
    scale_y_reverse(limits = rev(yl_ach)) +
    labs(x = "Acetylcholine (uM)", y = "Relaxation (% Phe pre-contraction)",
         title = "Endothelium-dependent relaxation (ACh)", colour = NULL) +
    theme_mg

  # SNP — inverted Y
  snp_long <- to_long(snp_pct, snp_conc)
  yl_snp   <- ylim_auto(snp_long$value)
  p_snp <- ggplot(snp_long, aes(dose, value, colour = channel, group = channel)) +
    geom_line() +
    geom_point(shape = 21, fill = "white", size = 2.5, stroke = 1.2) +
    scale_x_log10(labels = label_scientific()) +
    scale_colour_manual(values = setNames(CH_COLORS, paste0("X", 1:8)),
                        labels = paste0("Ch ", 1:8)) +
    scale_y_reverse(limits = rev(yl_snp)) +
    labs(x = "SNP (uM)", y = "Relaxation (% Phe pre-contraction)",
         title = "Endothelium-independent relaxation (SNP)", colour = NULL) +
    theme_mg

  # KCl barplot
  kcl_df <- data.frame(channel = paste0("Ch ", 1:8),
                        value   = as.numeric(kcl60))
  p_kcl <- ggplot(kcl_df, aes(channel, value, fill = channel)) +
    geom_col(width = 0.65, colour = "white") +
    geom_text(aes(label = round(value, 2)), vjust = -0.4, size = 3) +
    scale_fill_manual(values = setNames(CH_COLORS, paste0("Ch ", 1:8))) +
    labs(x = NULL, y = "Contraction (mN)", title = "KCl 60 mM — viability check") +
    theme_mg + theme(legend.position = "none",
                     axis.text.x = element_text(angle = 30, hjust = 1))

  combined <- (p_phe | p_ach) / (p_snp | p_kcl) +
    plot_annotation(
      title = paste0("Wire Myography — ", basename(INPUT_FILE), "  [", phe_mode, "]"),
      theme = theme(plot.title = element_text(size = 12, face = "bold"))
    )

  print(combined)
  out_png <- file.path(OUTPUT_DIR, paste0(tools::file_path_sans_ext(basename(INPUT_FILE)), ".png"))
  ggsave(out_png, combined, width = 12, height = 9, dpi = 150)
  cat(sprintf("Plot saved: %s\n", out_png))
}

# =============================================================================
# Run analysis, preview and download results
# =============================================================================

cat("Running analysis...\n")
results <- analyse_experiment(df, PHE_MODE, PHE_LABELS, ACH_LABELS, SNP_LABELS)
View(results)

plot_summary(results, PHE_MODE, PHE_CONC, ACH_CONC, SNP_CONC)
out_xlsx <- file.path(OUTPUT_DIR, paste0(tools::file_path_sans_ext(basename(INPUT_FILE)), ".xlsx"))
write_xlsx(list(results = results), path = out_xlsx)
cat(sprintf("Excel saved: %s\n", out_xlsx))
