#!/usr/bin/env Rscript
# =============================================================================
# SKRST Simulation Study — Classical Regime Only
# Saves BOTH:
#   1) replication-level raw results
#   2) scenario-level summary with mean / sd / se
#
# USAGE:
#   Rscript SKRST_simulation.R [n_cores] [n_rep]
#
# EXAMPLES:
#   Rscript SKRST_simulation.R
#   Rscript SKRST_simulation.R 8 1000
# =============================================================================

# =============================================================================
# SECTION 0 — COMMAND-LINE CONFIGURATION
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

default_cores <- 12L
default_rep   <- 1000L

n_cores <- if (length(args) >= 1) as.integer(args[1]) else default_cores
n_rep   <- if (length(args) >= 2) as.integer(args[2]) else default_rep

n_cores <- max(1L, min(n_cores, parallel::detectCores()))
n_rep   <- max(1L, n_rep)

cat(sprintf("Config: n_cores=%d | n_rep=%d\n\n", n_cores, n_rep))

# =============================================================================
# SECTION 1 — PACKAGE LOADING
# =============================================================================

required_pkgs <- c("MASS", "parallel")
optional_pkgs <- c("sn", "glasso")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Required package '%s' not found. Install it and retry.", pkg))
  }
  library(pkg, character.only = TRUE)
}

has_sn <- requireNamespace("sn", quietly = TRUE)
if (has_sn) library(sn)

has_glasso <- requireNamespace("glasso", quietly = TRUE)
if (has_glasso) library(glasso)

# =============================================================================
# SECTION 2 — UTILITIES
# =============================================================================

make_pd <- function(S, eps = 1e-8) {
  S <- (S + t(S)) / 2
  eig <- eigen(S, symmetric = TRUE)
  eig$values[eig$values < eps] <- eps
  V <- eig$vectors %*% diag(eig$values, nrow = length(eig$values)) %*% t(eig$vectors)
  (V + t(V)) / 2
}

safe_solve <- function(S) {
  tryCatch(solve(make_pd(S)), error = function(e) MASS::ginv(make_pd(S)))
}

frobenius_loss <- function(S_hat, S_true) {
  sum((S_hat - S_true)^2)
}

stein_loss <- function(S_hat, S_true) {
  tryCatch({
    M <- make_pd(S_true) %*% safe_solve(S_hat)
    ev <- Re(eigen(M, only.values = TRUE)$values)
    ev <- ev[ev > 0]
    sum(ev) - sum(log(ev)) - nrow(S_true)
  }, error = function(e) NA_real_)
}

safe_mean <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) NA_real_ else sd(x)
}

safe_se <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) <= 1) NA_real_ else sd(x) / sqrt(length(x))
}

make_toeplitz <- function(p, rho) {
  rho ^ abs(outer(seq_len(p), seq_len(p), "-"))
}

# Proposition 2.4:
# Sigma = c_nu * Omega - s_nu^2 * delta delta^T
true_sigma_mvst <- function(Omega, alpha_vec, nu) {
  stopifnot(nu > 2)
  c_nu  <- nu / (nu - 2)
  s_nu  <- sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2)
  denom <- sqrt(1 + as.numeric(t(alpha_vec) %*% Omega %*% alpha_vec))
  delta <- as.vector(Omega %*% alpha_vec) / denom
  make_pd(c_nu * Omega - s_nu^2 * outer(delta, delta))
}

# =============================================================================
# SECTION 3 — DATA GENERATION
# =============================================================================

# Full multivariate skew-t generator
# Tries sn::rmst() first, falls back to stochastic representation.
generate_mvst <- function(n, Omega, alpha_vec, nu) {
  p <- nrow(Omega)

  if (has_sn) {
    try_fit <- tryCatch(
      sn::rmst(n, xi = rep(0, p), Omega = Omega, alpha = alpha_vec, nu = nu),
      error = function(e) NULL
    )
    if (!is.null(try_fit)) return(try_fit)
  }

  denom <- sqrt(1 + as.numeric(t(alpha_vec) %*% Omega %*% alpha_vec))
  delta <- as.vector(Omega %*% alpha_vec) / denom
  Gamma <- make_pd(Omega - outer(delta, delta))
  GSq   <- t(chol(Gamma))

  W  <- rchisq(n, df = nu)
  Z0 <- rnorm(n)
  Y  <- matrix(rnorm(n * p), n, p)
  V  <- sqrt(nu / W)

  U <- outer(abs(Z0), delta) + tcrossprod(Y, GSq)
  V * U
}

# =============================================================================
# SECTION 4 — RESTRICTION FRAMEWORK
# =============================================================================
# Band restriction: entries with |i-j| > band_width are fixed to true Sigma values

make_restriction <- function(Sigma_true, band_width = 1) {
  p   <- nrow(Sigma_true)
  idx <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  far <- abs(idx[, 1] - idx[, 2]) > band_width

  list(
    pairs  = idx[far, , drop = FALSE],
    values = Sigma_true[idx[far, , drop = FALSE]]
  )
}

project_R <- function(S, restr) {
  S_R <- S
  nr  <- nrow(restr$pairs)

  if (nr > 0) {
    for (k in seq_len(nr)) {
      i <- restr$pairs[k, 1]
      j <- restr$pairs[k, 2]
      S_R[i, j] <- restr$values[k]
      S_R[j, i] <- restr$values[k]
    }
  }

  make_pd(S_R)
}

free_mask <- function(p, restr) {
  m  <- matrix(TRUE, p, p)
  nr <- nrow(restr$pairs)

  if (nr > 0) {
    for (k in seq_len(nr)) {
      i <- restr$pairs[k, 1]
      j <- restr$pairs[k, 2]
      m[i, j] <- FALSE
      m[j, i] <- FALSE
    }
  }

  m
}

# =============================================================================
# SECTION 5 — VARIANCE TERM ESTIMATOR
# =============================================================================
# Uses 1/n convention consistently inside LW and SKRST.

estimate_kappa <- function(X, restr) {
  n  <- nrow(X)
  p  <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n
  free_v <- as.vector(free_mask(p, restr))

  acc <- 0
  for (i in seq_len(n)) {
    xi  <- Xc[i, ]
    d   <- as.vector(outer(xi, xi)) - as.vector(S)
    acc <- acc + sum(d[free_v]^2)
  }

  acc / n^2
}

# =============================================================================
# SECTION 6 — ESTIMATORS
# =============================================================================

sample_cov_estimator <- function(X) {
  make_pd(cov(X))
}

# Ledoit-Wolf shrinkage toward scaled identity
# Internal 1/n convention, final output on 1/(n-1) scale.
lw_cov <- function(X) {
  n  <- nrow(X)
  p  <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n

  mu   <- mean(diag(S))
  Tmat <- mu * diag(p)

  pi_hat <- sum(sapply(seq_len(n), function(i) {
    xi <- Xc[i, ]
    sum((outer(xi, xi) - S)^2)
  })) / n^2

  gamma_hat <- sum((S - Tmat)^2)
  alpha_lw  <- if (gamma_hat < 1e-12) 1 else min(1, pi_hat / gamma_hat)

  lw <- (1 - alpha_lw) * S + alpha_lw * Tmat

  list(
    cov    = make_pd(lw * (n / (n - 1))),
    lambda = alpha_lw,
    S      = S,
    Tmat   = Tmat
  )
}

# Robust covariance + shrinkage
robust_shrink_cov <- function(X) {
  n <- nrow(X)
  p <- ncol(X)

  rob <- tryCatch(MASS::cov.rob(X, method = "mcd"), error = function(e) NULL)
  if (is.null(rob) || any(!is.finite(rob$cov))) {
    rob <- tryCatch(MASS::cov.rob(X, method = "mve"), error = function(e) NULL)
  }

  if (is.null(rob) || any(!is.finite(rob$cov))) {
    lw <- lw_cov(X)
    return(list(cov = lw$cov, lambda = NA_real_))
  }

  S <- make_pd(rob$cov)
  center <- rob$center

  mu   <- mean(diag(S))
  Tmat <- mu * diag(p)
  Xc   <- sweep(X, 2, center, "-")

  pi_hat <- sum(sapply(seq_len(n), function(i) {
    xi <- Xc[i, ]
    sum((outer(xi, xi) - S)^2)
  })) / n^2

  gamma_hat <- sum((S - Tmat)^2)
  alpha_rb  <- if (gamma_hat < 1e-12) 1 else min(1, pi_hat / gamma_hat)

  rob_shr <- (1 - alpha_rb) * S + alpha_rb * Tmat

  list(
    cov = make_pd(rob_shr),
    lambda = alpha_rb
  )
}

# Graphical lasso
glasso_cov <- function(X, rho = NULL) {
  S <- make_pd(cov(X))
  p <- ncol(X)

  if (is.null(rho)) {
    rho <- 0.05 * mean(diag(S))
  }

  if (!has_glasso) {
    return(list(cov = S, rho = rho))
  }

  fit <- tryCatch(glasso::glasso(S, rho = rho), error = function(e) NULL)

  if (is.null(fit) || is.null(fit$w) || any(!is.finite(fit$w))) {
    return(list(cov = S, rho = rho))
  }

  list(
    cov = make_pd(fit$w),
    rho = rho
  )
}

restricted_cov <- function(X, restr) {
  project_R(cov(X), restr)
}

# SKRST
skrst_cov <- function(X, restr, Sigma0 = NULL) {
  n  <- nrow(X)
  p  <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n

  if (is.null(Sigma0)) {
    sigma2 <- sum(diag(S)) / p
    Sigma0 <- sigma2 * diag(p)
  }

  Sigma0_R   <- project_R(Sigma0, restr)
  SigmaR_hat <- project_R(S, restr)

  kappa_hat <- estimate_kappa(X, restr)
  denom     <- sum((SigmaR_hat - Sigma0_R)^2)
  num       <- denom - kappa_hat

  lambda_hat <- if (denom > 1e-14) max(0, min(1, num / denom)) else 0

  Sigma_hat_1n <- lambda_hat * SigmaR_hat + (1 - lambda_hat) * Sigma0_R
  Sigma_hat    <- make_pd(Sigma_hat_1n * (n / (n - 1)))

  list(
    cov       = Sigma_hat,
    lambda    = lambda_hat,
    SigmaR    = SigmaR_hat,
    Sigma0_R  = Sigma0_R,
    kappa_hat = kappa_hat
  )
}

# Oracle SKRST
oracle_skrst_cov <- function(X, Sigma_true, restr, Sigma0 = NULL) {
  n  <- nrow(X)
  p  <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n

  if (is.null(Sigma0)) {
    sigma2 <- sum(diag(S)) / p
    Sigma0 <- sigma2 * diag(p)
  }

  Sigma0_R     <- project_R(Sigma0, restr)
  SigmaR_hat   <- project_R(S, restr)
  Sigma_true_1n <- Sigma_true * ((n - 1) / n)

  A   <- sum((Sigma_true_1n - Sigma0_R)^2)
  B   <- max(0, estimate_kappa(X, restr))
  lam <- if ((A + B) > 1e-14) max(0, min(1, A / (A + B))) else 1

  make_pd((lam * SigmaR_hat + (1 - lam) * Sigma0_R) * (n / (n - 1)))
}

# =============================================================================
# SECTION 7 — SINGLE REPLICATION
# =============================================================================

one_rep <- function(n, p, alpha0, nu, toeplitz_rho, band_width, rep_id = NA_integer_) {
  Omega      <- make_toeplitz(p, toeplitz_rho)
  alpha_vec  <- rep(alpha0, p)
  Sigma_true <- true_sigma_mvst(Omega, alpha_vec, nu)
  restr      <- make_restriction(Sigma_true, band_width)

  X <- generate_mvst(n, Omega, alpha_vec, nu)

  lw_fit <- lw_cov(X)
  rb_fit <- robust_shrink_cov(X)
  gl_fit <- glasso_cov(X)

  sk_fit <- tryCatch(
    skrst_cov(X, restr),
    error = function(e) list(cov = cov(X), lambda = NA_real_)
  )

  oc_fit <- tryCatch(
    oracle_skrst_cov(X, Sigma_true, restr),
    error = function(e) cov(X)
  )

  ests <- list(
    Sample        = sample_cov_estimator(X),
    LW            = lw_fit$cov,
    RestrictedMLE = restricted_cov(X, restr),
    Robust        = rb_fit$cov,
    Glasso        = gl_fit$cov,
    SKRST         = sk_fit$cov,
    Oracle        = oc_fit
  )

  FL <- sapply(ests, frobenius_loss, S_true = Sigma_true)
  SL <- sapply(ests, stein_loss,     S_true = Sigma_true)

  data.frame(
    rep = rep_id,
    n = n,
    p = p,
    alpha0 = alpha0,
    nu = nu,
    toeplitz_rho = toeplitz_rho,
    band_width = band_width,

    FL_Sample        = unname(FL["Sample"]),
    FL_LW            = unname(FL["LW"]),
    FL_RestrictedMLE = unname(FL["RestrictedMLE"]),
    FL_Robust        = unname(FL["Robust"]),
    FL_Glasso        = unname(FL["Glasso"]),
    FL_SKRST         = unname(FL["SKRST"]),
    FL_Oracle        = unname(FL["Oracle"]),

    SL_Sample        = unname(SL["Sample"]),
    SL_LW            = unname(SL["LW"]),
    SL_RestrictedMLE = unname(SL["RestrictedMLE"]),
    SL_Robust        = unname(SL["Robust"]),
    SL_Glasso        = unname(SL["Glasso"]),
    SL_SKRST         = unname(SL["SKRST"]),
    SL_Oracle        = unname(SL["Oracle"]),

    lambda_LW = lw_fit$lambda,
    lambda_Robust = rb_fit$lambda,
    lambda_SKRST = sk_fit$lambda,
    rho_Glasso = gl_fit$rho,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# SECTION 8 — MONTE CARLO RUNNER
# =============================================================================

run_scenario <- function(n, p, alpha0, nu, toeplitz_rho, band_width,
                         R, seed, n_cores) {
  if (p >= n - 2 || p < 3) return(NULL)

  raw_list <- parallel::mclapply(
    seq_len(R),
    function(r) {
      set.seed(seed + r)
      tryCatch(
        one_rep(n, p, alpha0, nu, toeplitz_rho, band_width, rep_id = r),
        error = function(e) data.frame(
          rep = r,
          n = n,
          p = p,
          alpha0 = alpha0,
          nu = nu,
          toeplitz_rho = toeplitz_rho,
          band_width = band_width,

          FL_Sample = NA_real_, FL_LW = NA_real_, FL_RestrictedMLE = NA_real_,
          FL_Robust = NA_real_, FL_Glasso = NA_real_, FL_SKRST = NA_real_, FL_Oracle = NA_real_,

          SL_Sample = NA_real_, SL_LW = NA_real_, SL_RestrictedMLE = NA_real_,
          SL_Robust = NA_real_, SL_Glasso = NA_real_, SL_SKRST = NA_real_, SL_Oracle = NA_real_,

          lambda_LW = NA_real_, lambda_Robust = NA_real_,
          lambda_SKRST = NA_real_, rho_Glasso = NA_real_,
          stringsAsFactors = FALSE
        )
      )
    },
    mc.cores = n_cores
  )

  raw_df <- do.call(rbind, raw_list)

  metric_cols <- c(
    "FL_Sample", "FL_LW", "FL_RestrictedMLE", "FL_Robust", "FL_Glasso", "FL_SKRST", "FL_Oracle",
    "SL_Sample", "SL_LW", "SL_RestrictedMLE", "SL_Robust", "SL_Glasso", "SL_SKRST", "SL_Oracle",
    "lambda_LW", "lambda_Robust", "lambda_SKRST", "rho_Glasso"
  )

  summary_list <- list(
    n = n,
    p = p,
    alpha0 = alpha0,
    nu = nu,
    toeplitz_rho = toeplitz_rho,
    band_width = band_width,
    n_rep = R
  )

  for (nm in metric_cols) {
    summary_list[[paste0(nm, "__mean")]] <- safe_mean(raw_df[[nm]])
    summary_list[[paste0(nm, "__sd")]]   <- safe_sd(raw_df[[nm]])
    summary_list[[paste0(nm, "__se")]]   <- safe_se(raw_df[[nm]])
  }

  summary_df <- as.data.frame(summary_list, stringsAsFactors = FALSE)

  list(raw = raw_df, summary = summary_df)
}

run_all <- function(grid, R, base_seed, n_cores, label) {
  cat(sprintf("\n=== %s: %d scenarios x %d reps ===\n", label, nrow(grid), R))
  t0 <- proc.time()

  raw_results <- vector("list", nrow(grid))
  summary_results <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]

    cat(sprintf("  [%d/%d] n=%d p=%d a=%.1f nu=%2d rho=%.1f bw=%d",
                i, nrow(grid), row$n, row$p, row$alpha0, row$nu,
                row$toeplitz_rho, row$band_width))

    res <- run_scenario(
      n = row$n,
      p = row$p,
      alpha0 = row$alpha0,
      nu = row$nu,
      toeplitz_rho = row$toeplitz_rho,
      band_width = row$band_width,
      R = R,
      seed = base_seed + i * 100L,
      n_cores = n_cores
    )

    if (!is.null(res)) {
      s <- res$summary

      cat(sprintf(
        " | FL_rel: Sn=%.3f LW=%.3f RMLE=%.3f Robust=%.3f Glasso=%.3f Oracle=%.3f | lam_SK=%.3f\n",
        s$FL_Sample__mean        / s$FL_SKRST__mean,
        s$FL_LW__mean            / s$FL_SKRST__mean,
        s$FL_RestrictedMLE__mean / s$FL_SKRST__mean,
        s$FL_Robust__mean        / s$FL_SKRST__mean,
        s$FL_Glasso__mean        / s$FL_SKRST__mean,
        s$FL_Oracle__mean        / s$FL_SKRST__mean,
        s$lambda_SKRST__mean
      ))

      raw_results[[i]] <- res$raw
      summary_results[[i]] <- res$summary
    } else {
      cat(" [SKIPPED]\n")
    }
  }

  elapsed <- proc.time() - t0
  cat(sprintf("  Done in %.1f minutes.\n", elapsed["elapsed"] / 60))

  list(
    raw = do.call(rbind, raw_results[!sapply(raw_results, is.null)]),
    summary = do.call(rbind, summary_results[!sapply(summary_results, is.null)])
  )
}

# =============================================================================
# SECTION 9 — SCENARIO GRID
# =============================================================================

alpha_vals <- c(0, 0.5, 1.0, 1.5)
nu_vals    <- c(5, 10, 20, 50)

grid_classical <- expand.grid(
  n = c(100, 200, 500, 1000),
  p = c(5, 10, 20),
  alpha0 = alpha_vals,
  nu = nu_vals,
  toeplitz_rho = 0.7,
  band_width = 1L,
  stringsAsFactors = FALSE
)

# =============================================================================
# SECTION 10 — RUN AND SAVE
# =============================================================================

base_seed <- 2024L
out_dir   <- "."

res_classical <- run_all(
  grid_classical,
  R = n_rep,
  base_seed = base_seed,
  n_cores = n_cores,
  label = "Classical regime"
)

f_raw     <- file.path(out_dir, "skrst_classical_raw.csv")
f_summary <- file.path(out_dir, "skrst_classical_summary.csv")
f_legacy  <- file.path(out_dir, "skrst_classical.csv")

write.csv(res_classical$raw,     f_raw,     row.names = FALSE)
write.csv(res_classical$summary, f_summary, row.names = FALSE)
write.csv(res_classical$summary, f_legacy,  row.names = FALSE)

cat(sprintf("  Saved raw reps: %s\n", f_raw))
cat(sprintf("  Saved summary : %s\n", f_summary))
cat(sprintf("  Saved legacy  : %s\n", f_legacy))

# =============================================================================
# SECTION 11 — SUMMARY PRINTER
# =============================================================================

print_summary <- function(res, label) {
  cat(sprintf("\n\n========== %s ==========\n", label))

  cat("\n--- Relative Frobenius Loss (ratio / FL_SKRST mean; >1 = SKRST better) ---\n")
  rel <- res[, c("n", "p", "alpha0", "nu")]
  rel$Sn     <- round(res$FL_Sample__mean        / res$FL_SKRST__mean, 3)
  rel$LW     <- round(res$FL_LW__mean            / res$FL_SKRST__mean, 3)
  rel$RMLE   <- round(res$FL_RestrictedMLE__mean / res$FL_SKRST__mean, 3)
  rel$Robust <- round(res$FL_Robust__mean        / res$FL_SKRST__mean, 3)
  rel$Glasso <- round(res$FL_Glasso__mean        / res$FL_SKRST__mean, 3)
  rel$Oracle <- round(res$FL_Oracle__mean        / res$FL_SKRST__mean, 3)
  print(rel, row.names = FALSE)

  cat("\n--- Relative Stein Loss (ratio / SL_SKRST mean; >1 = SKRST better) ---\n")
  rel_sl <- res[, c("n", "p", "alpha0", "nu")]
  rel_sl$Sn     <- round(res$SL_Sample__mean        / res$SL_SKRST__mean, 3)
  rel_sl$LW     <- round(res$SL_LW__mean            / res$SL_SKRST__mean, 3)
  rel_sl$RMLE   <- round(res$SL_RestrictedMLE__mean / res$SL_SKRST__mean, 3)
  rel_sl$Robust <- round(res$SL_Robust__mean        / res$SL_SKRST__mean, 3)
  rel_sl$Glasso <- round(res$SL_Glasso__mean        / res$SL_SKRST__mean, 3)
  rel_sl$Oracle <- round(res$SL_Oracle__mean        / res$SL_SKRST__mean, 3)
  print(rel_sl, row.names = FALSE)

  cat("\n--- Mean Frobenius Loss ---\n")
  fl_cols <- c(
    "n", "p", "alpha0", "nu",
    "FL_Sample__mean", "FL_LW__mean", "FL_RestrictedMLE__mean",
    "FL_Robust__mean", "FL_Glasso__mean", "FL_SKRST__mean", "FL_Oracle__mean"
  )
  print(round(res[, fl_cols], 4), row.names = FALSE)

  cat("\n--- SE Frobenius Loss ---\n")
  fl_se_cols <- c(
    "n", "p", "alpha0", "nu",
    "FL_Sample__se", "FL_LW__se", "FL_RestrictedMLE__se",
    "FL_Robust__se", "FL_Glasso__se", "FL_SKRST__se", "FL_Oracle__se"
  )
  print(round(res[, fl_se_cols], 4), row.names = FALSE)

  cat("\n--- Mean Stein Loss ---\n")
  sl_cols <- c(
    "n", "p", "alpha0", "nu",
    "SL_Sample__mean", "SL_LW__mean", "SL_RestrictedMLE__mean",
    "SL_Robust__mean", "SL_Glasso__mean", "SL_SKRST__mean", "SL_Oracle__mean"
  )
  print(round(res[, sl_cols], 4), row.names = FALSE)

  cat("\n--- SE Stein Loss ---\n")
  sl_se_cols <- c(
    "n", "p", "alpha0", "nu",
    "SL_Sample__se", "SL_LW__se", "SL_RestrictedMLE__se",
    "SL_Robust__se", "SL_Glasso__se", "SL_SKRST__se", "SL_Oracle__se"
  )
  print(round(res[, sl_se_cols], 4), row.names = FALSE)

  cat("\n--- Mean tuning parameters ---\n")
  lam_cols <- c(
    "n", "p", "alpha0", "nu",
    "lambda_SKRST__mean", "lambda_LW__mean", "lambda_Robust__mean", "rho_Glasso__mean"
  )
  print(round(res[, lam_cols], 3), row.names = FALSE)
}

print_summary(res_classical$summary, "CLASSICAL REGIME")

for (nv in nu_vals) {
  sub <- res_classical$summary[res_classical$summary$nu == nv, ]
  if (nrow(sub) == 0) next

  cat(sprintf("\n--- Publication Table: nu=%d, Frobenius Loss ---\n", nv))
  fl_cols <- c(
    "n", "p", "alpha0",
    "FL_Sample__mean", "FL_LW__mean", "FL_RestrictedMLE__mean",
    "FL_Robust__mean", "FL_Glasso__mean", "FL_SKRST__mean", "FL_Oracle__mean"
  )
  print(round(sub[, fl_cols], 4), row.names = FALSE)

  cat(sprintf("\n--- Publication Table: nu=%d, Stein Loss ---\n", nv))
  sl_cols <- c(
    "n", "p", "alpha0",
    "SL_Sample__mean", "SL_LW__mean", "SL_RestrictedMLE__mean",
    "SL_Robust__mean", "SL_Glasso__mean", "SL_SKRST__mean", "SL_Oracle__mean"
  )
  print(round(sub[, sl_cols], 4), row.names = FALSE)
}

cat("\nAll done.\n")
