###############################################################################
# SUPPLEMENTARY FIGURES:
# Relative loss of each estimator compared with SKRST
#
# Layout:
#   x-axis  : n
#   rows    : p
#   columns : alpha0
#   color   : nu
#
# Y-axis:
#   Relative loss = Loss(estimator) / Loss(SKRST)
#
# Interpretation:
#   1   = same as SKRST
#   > 1 = estimator is worse than SKRST
#
# Output:
#   One PDF per estimator for Frobenius loss
#   One PDF per estimator for Stein loss
###############################################################################
# =========================
# 1. Packages
# =========================
required_packages <- c("readr", "dplyr", "ggplot2")
install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    install.packages(new_packages, dependencies = TRUE)
  }
}
install_if_missing(required_packages)
library(readr)
library(dplyr)
library(ggplot2)
# =========================
# 2. Settings
# =========================
input_file <- "simulation_results_server_raw.csv"
output_dir <- "supplementary_relative_loss_figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
# Optional: use log scale for Stein figures
use_log_scale_frobenius <- FALSE
use_log_scale_stein     <- TRUE
# =========================
# 3. Load raw data
# =========================
raw <- read_csv(input_file, show_col_types = FALSE)
# =========================
# 4. Estimator map
# =========================
estimator_map <- data.frame(
  short_name = c("Sample", "LW", "RestrictedMLE", "Robust", "Glasso", "Oracle"),
  label = c("Sample Covariance", "Ledoit-Wolf", "Restricted MLE",
            "Robust Shrinkage", "Glasso", "Oracle"),
  frob_col = c("FL_Sample", "FL_LW", "FL_RestrictedMLE", "FL_Robust", "FL_Glasso", "FL_Oracle"),
  stein_col = c("SL_Sample", "SL_LW", "SL_RestrictedMLE", "SL_Robust", "SL_Glasso", "SL_Oracle"),
  stringsAsFactors = FALSE
)
# =========================
# 5. Helper to summarise one relative-loss variable
# =========================
summarise_rel <- function(df, rel_col) {
  df %>%
    group_by(n, p, alpha0, nu) %>%
    summarise(
      mean_rel = mean(.data[[rel_col]], na.rm = TRUE),
      se_rel   = sd(.data[[rel_col]], na.rm = TRUE) / sqrt(sum(!is.na(.data[[rel_col]]))),
      .groups = "drop"
    ) %>%
    mutate(
      n = factor(n, levels = c(100, 200, 500, 1000)),
      p = factor(p, levels = c(5, 10, 20)),
      alpha0_lab = factor(
        alpha0,
        levels = c(0.0, 0.5, 1.0, 1.5),
        labels = c("0", "0.5", "1.0", "1.5")
      ),
      nu = factor(nu, levels = c(5, 10, 20, 50))
    )
}
# =========================
# 6. Function to plot one estimator vs SKRST
# =========================
make_rel_plot <- function(df, estimator_label, loss_label = c("Frobenius", "Stein"),
                          log_scale = FALSE) {
  loss_label <- match.arg(loss_label)
  p <- ggplot(
    df,
    aes(x = n, y = mean_rel, color = nu, group = nu)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.5, color = "black") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_rel - se_rel, 1e-8),
        ymax = mean_rel + se_rel
      ),
      width = 0.08,
      linewidth = 0.45
    ) +
    facet_grid(
      p ~ alpha0_lab,
      labeller = labeller(
        p = function(x) paste0("p = ", x),
        alpha0_lab = function(x) paste0("\u03B1", "\u2080", " = ", x)
      )
    ) +
    labs(
      title = paste0("Relative ", loss_label, " Loss: ", estimator_label, " vs SKRST"),
      subtitle = paste0("Values above 1 indicate worse ", tolower(loss_label), " loss than SKRST"),
      x = "Sample size n",
      y = paste0(loss_label, " loss / SKRST"),
      color = expression(nu)
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_text(size = 13),
      axis.title.y = element_text(size = 13),
      strip.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  if (log_scale) {
    p <- p + scale_y_log10()
  }
  p
}
# =========================
# 7. Build and save all Frobenius figures
# =========================
for (i in seq_len(nrow(estimator_map))) {
  est_short <- estimator_map$short_name[i]
  est_label <- estimator_map$label[i]
  est_col   <- estimator_map$frob_col[i]
  rel_col_name <- paste0("rel_", est_short, "_Frob")
  raw_tmp <- raw %>%
    mutate(!!rel_col_name := .data[[est_col]] / FL_SKRST)
  df_plot <- summarise_rel(raw_tmp, rel_col_name)
  fig <- make_rel_plot(
    df = df_plot,
    estimator_label = est_label,
    loss_label = "Frobenius",
    log_scale = use_log_scale_frobenius
  )
  out_file <- file.path(
    output_dir,
    paste0("fig_frobenius_", est_short, "_vs_skrst.pdf")
  )
  ggsave(
    filename = out_file,
    plot = fig,
    width = 12,
    height = 8,
    device = "pdf"
  )
  print(fig)
}
# =========================
# 8. Build and save all Stein figures
# =========================
for (i in seq_len(nrow(estimator_map))) {
  est_short <- estimator_map$short_name[i]
  est_label <- estimator_map$label[i]
  est_col   <- estimator_map$stein_col[i]
  rel_col_name <- paste0("rel_", est_short, "_Stein")
  raw_tmp <- raw %>%
    mutate(!!rel_col_name := .data[[est_col]] / SL_SKRST)
  df_plot <- summarise_rel(raw_tmp, rel_col_name)
  fig <- make_rel_plot(
    df = df_plot,
    estimator_label = est_label,
    loss_label = "Stein",
    log_scale = use_log_scale_stein
  )
  out_file <- file.path(
    output_dir,
    paste0("fig_stein_", est_short, "_vs_skrst.pdf")
  )
  ggsave(
    filename = out_file,
    plot = fig,
    width = 12,
    height = 8,
    device = "pdf"
  )
  print(fig)
}
# =========================
# 9. Final message
# =========================
cat("\nAll supplementary figures saved in:\n", normalizePath(output_dir), "\n")