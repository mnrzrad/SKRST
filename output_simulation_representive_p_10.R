###############################################################################
# CLEAN SCRIPT: TABLES AND FIGURES FOR SKRST SIMULATION STUDY
# Source file: skrst_classical_raw.csv
# Output: LaTeX tables + PDF figures only
###############################################################################

###############################################################################
# 1. PACKAGES
###############################################################################
required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "xtable")

install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages) > 0) {
    install.packages(new_packages, dependencies = TRUE)
  }
}
install_if_missing(required_packages)

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(xtable)

###############################################################################
# 2. SETTINGS
###############################################################################
input_file <- "simulation_results_server_raw.csv"
output_dir <- "mc_tables_figures"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

###############################################################################
# 3. LOAD DATA
###############################################################################
mc <- read_csv(input_file, show_col_types = FALSE)

###############################################################################
# 4. CLEAN STEIN COLUMNS AND DERIVED QUANTITIES
###############################################################################
stein_cols <- c(
  "SL_Sample", "SL_LW", "SL_RestrictedMLE",
  "SL_Robust", "SL_Glasso", "SL_SKRST", "SL_Oracle"
)

mc <- mc %>%
  mutate(across(all_of(stein_cols), ~ ifelse(is.finite(.x), .x, NA_real_))) %>%
  mutate(
    gain_SKRST_vs_LW         = 100 * (FL_LW - FL_SKRST) / FL_LW,
    gain_SKRST_vs_Restricted = 100 * (FL_RestrictedMLE - FL_SKRST) / FL_RestrictedMLE,
    gain_SKRST_vs_Robust     = 100 * (FL_Robust - FL_SKRST) / FL_Robust,
    gain_SKRST_vs_Glasso     = 100 * (FL_Glasso - FL_SKRST) / FL_Glasso,
    gain_SKRST_vs_Sample     = 100 * (FL_Sample - FL_SKRST) / FL_Sample,
    gap_to_oracle            = FL_SKRST - FL_Oracle
  )

###############################################################################
# 5. HELPERS
###############################################################################
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

safe_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) NA_real_ else sd(x) / sqrt(length(x))
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) NA_real_ else sd(x)
}


fmt_mean_se <- function(mean, se, digits_mean = 3, digits_se = 3) {
  out <- ifelse(
    is.na(mean) | !is.finite(mean),
    "NA",
    ifelse(
      is.na(se) | !is.finite(se),
      sprintf(paste0("%.", digits_mean, "f"), mean),
      sprintf(paste0("%.", digits_mean, "f (%.", digits_se, "f)"), mean, se)
    )
  )
  as.character(out)
}

save_tex_table <- function(df, file_name, caption, label) {
  out_file <- file.path(output_dir, file_name)
  
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(x) as.character(unlist(x)))
  
  print(
    xtable(df, caption = caption, label = label),
    file = out_file,
    include.rownames = FALSE,
    booktabs = TRUE,
    floating = TRUE,
    sanitize.text.function = identity,
    comment = FALSE
  )
}

summarise_mc <- function(data, group_vars, value_vars) {
  data %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      across(
        all_of(value_vars),
        list(
          mean = ~ safe_mean(.x),
          sd   = ~ safe_sd(.x),
          se   = ~ safe_se(.x)
        ),
        .names = "{.col}__{.fn}"
      ),
      .groups = "drop"
    )
}

make_mean_se_table <- function(data, filter_expr, col_var, value_vars, labels,
                               digits_mean = 3, digits_se = 3) {
  sub <- subset(data, eval(filter_expr, data, parent.frame()))
  if (nrow(sub) == 0) stop("No rows matched the requested filter.")
  
  out <- data.frame(
    Estimator = labels,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  col_levels <- sort(unique(sub[[col_var]]))
  
  for (cv in col_levels) {
    vals <- sapply(value_vars, function(v) {
      z <- sub[[v]][sub[[col_var]] == cv]
      m <- safe_mean(z)
      s <- safe_se(z)
      fmt_mean_se(m, s, digits_mean = digits_mean, digits_se = digits_se)
    })
    out[[as.character(cv)]] <- vals
  }
  
  out
}

###############################################################################
# 6. ESTIMATOR LABELS
###############################################################################
fl_vars <- c(
  "FL_Sample", "FL_LW", "FL_RestrictedMLE", "FL_Robust", "FL_Glasso", "FL_SKRST", "FL_Oracle"
)

sl_vars <- c(
  "SL_Sample", "SL_LW", "SL_RestrictedMLE", "SL_Robust", "SL_Glasso", "SL_SKRST", "SL_Oracle"
)

labels <- c(
  "Sample Covariance",
  "Ledoit-Wolf",
  "Restricted MLE",
  "Robust Shrinkage",
  "Glasso",
  "SKRST",
  "Oracle"
)

###############################################################################
# 7. MAIN TEXT TABLES — FROBENIUS
###############################################################################
tab_frob_n <- make_mean_se_table(
  data = mc,
  filter_expr = quote(p == 10 & alpha0 == 1.0 & nu == 10),
  col_var = "n",
  value_vars = fl_vars,
  labels = labels
)

save_tex_table(
  tab_frob_n,
  "table_mc_main_benchmark.tex",
  "Mean Frobenius loss with Monte Carlo standard errors in parentheses for the benchmark setting $p = 10$, $\\alpha_0 = 1.0$, $\\nu = 10$, and varying sample size $n$.",
  "tab:mc_main_benchmark"
)

tab_frob_nu <- make_mean_se_table(
  data = mc,
  filter_expr = quote(p == 10 & n == 200 & alpha0 == 1.0),
  col_var = "nu",
  value_vars = fl_vars,
  labels = labels
)

save_tex_table(
  tab_frob_nu,
  "table_mc_tail_weight.tex",
  "Mean Frobenius loss with Monte Carlo standard errors in parentheses for varying degrees of freedom $\\nu$ with $p = 10$, $n = 200$, and $\\alpha_0 = 1.0$.",
  "tab:mc_tail_weight"
)

tab_frob_p <- make_mean_se_table(
  data = mc,
  filter_expr = quote(n == 200 & alpha0 == 1.0 & nu == 10),
  col_var = "p",
  value_vars = fl_vars,
  labels = labels
)

save_tex_table(
  tab_frob_p,
  "table_mc_dimensionality.tex",
  "Mean Frobenius loss with Monte Carlo standard errors in parentheses for varying dimension $p$ with $n = 200$, $\\alpha_0 = 1.0$, and $\\nu = 10$.",
  "tab:mc_dimensionality"
)

###############################################################################
# 8. MAIN TEXT TABLES — STEIN
###############################################################################
tab_stein_n <- make_mean_se_table(
  data = mc,
  filter_expr = quote(p == 10 & alpha0 == 1.0 & nu == 10),
  col_var = "n",
  value_vars = sl_vars,
  labels = labels
)

save_tex_table(
  tab_stein_n,
  "table_mc_main_benchmark_stein.tex",
  "Mean Stein loss with Monte Carlo standard errors in parentheses for the benchmark setting $p = 10$, $\\alpha_0 = 1.0$, $\\nu = 10$, and varying sample size $n$.",
  "tab:mc_main_benchmark_stein"
)

tab_stein_nu <- make_mean_se_table(
  data = mc,
  filter_expr = quote(p == 10 & n == 200 & alpha0 == 1.0),
  col_var = "nu",
  value_vars = sl_vars,
  labels = labels
)

save_tex_table(
  tab_stein_nu,
  "table_mc_tail_weight_stein.tex",
  "Mean Stein loss with Monte Carlo standard errors in parentheses for varying degrees of freedom $\\nu$ with $p = 10$, $n = 200$, and $\\alpha_0 = 1.0$.",
  "tab:mc_tail_weight_stein"
)

tab_stein_p <- make_mean_se_table(
  data = mc,
  filter_expr = quote(n == 200 & alpha0 == 1.0 & nu == 10),
  col_var = "p",
  value_vars = sl_vars,
  labels = labels
)

save_tex_table(
  tab_stein_p,
  "table_mc_dimensionality_stein.tex",
  "Mean Stein loss with Monte Carlo standard errors in parentheses for varying dimension $p$ with $n = 200$, $\\alpha_0 = 1.0$, and $\\nu = 10$.",
  "tab:mc_dimensionality_stein"
)

###############################################################################
# 9. MAIN TEXT FIGURES — FROBENIUS AND STEIN
###############################################################################
frob_sum <- summarise_mc(
  data = mc %>% filter(p == 10, nu == 10),
  group_vars = c("n", "alpha0"),
  value_vars = fl_vars
)

frob_df <- frob_sum %>%
  pivot_longer(
    cols = -(n:alpha0),
    names_to = c("Estimator", ".value"),
    names_pattern = "(FL_[^_]+(?:MLE|Robust|Glasso|SKRST|Oracle|Sample|LW)?)__(mean|se)"
  ) %>%
  mutate(
    Estimator = factor(Estimator, levels = fl_vars, labels = labels),
    n = factor(n, levels = c(100, 200, 500, 1000)),
    alpha0 = factor(alpha0, levels = c(0.0, 0.5, 1.0, 1.5),
                    labels = c("0", "0.5", "1.0", "1.5"))
  )

fig_frob <- ggplot(frob_df, aes(x = n, y = mean, color = Estimator, group = Estimator)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.08, linewidth = 0.45) +
  facet_wrap(~ alpha0, nrow = 1, labeller = label_bquote(alpha[0] == .(alpha0))) +
  labs(
    title = "Frobenius Loss Across Sample Sizes",
    x = "Sample size n",
    y = "Frobenius loss",
    color = "Estimator"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "fig_frobenius_vs_n.pdf"),
  fig_frob,
  width = 11,
  height = 4.8,
  device = "pdf"
)

stein_sum <- summarise_mc(
  data = mc %>% filter(p == 10, nu == 10),
  group_vars = c("n", "alpha0"),
  value_vars = sl_vars
)

stein_df <- stein_sum %>%
  pivot_longer(
    cols = -(n:alpha0),
    names_to = c("Estimator", ".value"),
    names_pattern = "(SL_[^_]+(?:MLE|Robust|Glasso|SKRST|Oracle|Sample|LW)?)__(mean|se)"
  ) %>%
  mutate(
    Estimator = factor(Estimator, levels = sl_vars, labels = labels),
    n = factor(n, levels = c(100, 200, 500, 1000)),
    alpha0 = factor(alpha0, levels = c(0.0, 0.5, 1.0, 1.5),
                    labels = c("0", "0.5", "1.0", "1.5"))
  )

fig_stein <- ggplot(stein_df, aes(x = n, y = mean, color = Estimator, group = Estimator)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(mean - se, 1e-8), ymax = mean + se), width = 0.08, linewidth = 0.45) +
  facet_wrap(~ alpha0, nrow = 1, labeller = label_bquote(alpha[0] == .(alpha0))) +
  labs(
    title = "Stein Loss Across Sample Sizes",
    x = "Sample size n",
    y = "Stein loss",
    color = "Estimator"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "fig_stein_vs_n.pdf"),
  fig_stein,
  width = 11,
  height = 4.8,
  device = "pdf"
)

###############################################################################
# 10. MAIN/SUPP FIGURE — TUNING PARAMETERS
###############################################################################
fig_tune_df <- mc %>%
  group_by(alpha0, nu) %>%
  summarise(
    lambda_LW = safe_mean(lambda_LW),
    lambda_Robust = safe_mean(lambda_Robust),
    lambda_SKRST = safe_mean(lambda_SKRST),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(lambda_LW, lambda_Robust, lambda_SKRST),
    names_to = "Parameter",
    values_to = "Value"
  ) %>%
  mutate(
    alpha0 = factor(alpha0, levels = c(0.0, 0.5, 1.0, 1.5),
                    labels = c("0", "0.5", "1.0", "1.5")),
    nu = factor(nu, levels = c(5, 10, 20, 50))
  )

fig_tune <- ggplot(fig_tune_df, aes(x = nu, y = Value, color = Parameter, group = Parameter)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  facet_wrap(~ alpha0, nrow = 1, labeller = label_bquote(alpha[0] == .(alpha0))) +
  labs(
    title = "Average Shrinkage Intensities Across Scenarios",
    x = expression(nu),
    y = "Average tuning value",
    color = "Parameter"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom"
  )

ggsave(
  file.path(output_dir, "fig_tuning_parameters.pdf"),
  fig_tune,
  width = 11,
  height = 4.8,
  device = "pdf"
)

###############################################################################
# 11. APPENDIX TABLES — FROBENIUS BY (p, alpha0), rows = estimator/nu, cols = n
###############################################################################
appendix_fl_raw <- summarise_mc(
  data = mc,
  group_vars = c("p", "n", "alpha0", "nu"),
  value_vars = fl_vars
)

make_appendix_fl_table <- function(p_val, alpha_val) {
  appendix_fl_raw %>%
    filter(p == p_val, alpha0 == alpha_val) %>%
    dplyr::select(p, n, alpha0, nu, matches("^(FL_).*(__(mean|se))$")) %>%
    pivot_longer(
      cols = -(p:n:alpha0:nu),
      names_to = c("Estimator", ".value"),
      names_pattern = "(FL_[^_]+(?:MLE|Robust|Glasso|SKRST|Oracle|Sample|LW)?)__(mean|se)"
    ) %>%
    mutate(
      Estimator = factor(Estimator, levels = fl_vars, labels = labels),
      value = fmt_mean_se(mean, se),
      nu = paste0("nu=", nu)
    ) %>%
    dplyr::select(Estimator, nu, n, value) %>%
    pivot_wider(names_from = n, values_from = value, values_fn = \(x) x[1]) %>%
    arrange(Estimator, nu) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

###############################################################################
# 12. APPENDIX TABLES — STEIN BY (p, alpha0), rows = estimator/nu, cols = n
###############################################################################
appendix_sl_raw <- summarise_mc(
  data = mc,
  group_vars = c("p", "n", "alpha0", "nu"),
  value_vars = sl_vars
)

make_appendix_sl_table <- function(p_val, alpha_val) {
  appendix_sl_raw %>%
    filter(p == p_val, alpha0 == alpha_val) %>%
    dplyr::select(p, n, alpha0, nu, matches("^(SL_).*(__(mean|se))$")) %>%
    pivot_longer(
      cols = -(p:n:alpha0:nu),
      names_to = c("Estimator", ".value"),
      names_pattern = "(SL_[^_]+(?:MLE|Robust|Glasso|SKRST|Oracle|Sample|LW)?)__(mean|se)"
    ) %>%
    mutate(
      Estimator = factor(Estimator, levels = sl_vars, labels = labels),
      value = fmt_mean_se(mean, se),
      nu = paste0("nu=", nu)
    ) %>%
    dplyr::select(Estimator, nu, n, value) %>%
    pivot_wider(names_from = n, values_from = value, values_fn = \(x) x[1]) %>%
    arrange(Estimator, nu) %>%
    as.data.frame(stringsAsFactors = FALSE)
}

p_vals <- sort(unique(mc$p))
alpha_vals <- sort(unique(mc$alpha0))

for (pp in p_vals) {
  for (aa in alpha_vals) {
    tab_fl <- make_appendix_fl_table(pp, aa)
    save_tex_table(
      tab_fl,
      paste0("appendix_table_frob_p", pp, "_a", gsub("\\.", "_", aa), ".tex"),
      paste0("Mean Frobenius loss with Monte Carlo standard errors in parentheses for p = ", pp,
             " and $\\alpha_0 = ", aa, "$, across values of $n$ and $\\nu$."),
      paste0("tab:appendix_frob_p", pp, "_a", gsub("\\.", "_", aa))
    )
    
    tab_sl <- make_appendix_sl_table(pp, aa)
    save_tex_table(
      tab_sl,
      paste0("appendix_table_stein_p", pp, "_a", gsub("\\.", "_", aa), ".tex"),
      paste0("Mean Stein loss with Monte Carlo standard errors in parentheses for p = ", pp,
             " and $\\alpha_0 = ", aa, "$, across values of $n$ and $\\nu$."),
      paste0("tab:appendix_stein_p", pp, "_a", gsub("\\.", "_", aa))
    )
  }
}

###############################################################################
# 13. FINAL MESSAGE
###############################################################################
cat("\nAll LaTeX tables and PDF figures saved in:\n", normalizePath(output_dir), "\n")