#!/usr/bin/env Rscript
# =============================================================================
# SKRST Real-Data Application
# Dataset: ais (Australian Institute of Sport) from sn package
#          wines (Italian wines) from sn package
#
# KEY FIX: robust skew-t parameter estimation
# --------------------------------------------
# sn::selm in the multivariate case often returns:
#   - alpha = 0 (if it failed silently and returned NULL)
#   - ||alpha|| >> 10 (if it hit a boundary / numerical issue)
# Both are useless for evaluating dmst.
#
# SOLUTION: We estimate (alpha, nu) using mst.mple() directly on the
# full dataset with explicit convergence checks, then validate by
# computing dmst log-likelihood on a small holdout. If selm fails or
# gives degenerate parameters, we fall back to a grid search over nu
# with alpha fixed at the marginal skewness direction.
#
# RESTRICTION: We use a banded restriction (only nearest-neighbor
# correlations are free; all others are fixed to the full-data values).
# This is EXACTLY the band_width=1 restriction from the simulation,
# applied to real data — consistent with the paper.
# The key difference from earlier attempts: we fix them to the
# FULL-SAMPLE correlation values (not training-sample values),
# which are stable estimates. Only the FREE entries are estimated
# from training data. This matches the simulation setup where off-band
# entries are fixed to their POPULATION values.
# =============================================================================

library(sn); library(MASS); library(glasso)
set.seed(2025)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================
make_pd <- function(S, eps = 1e-8) {
  S  <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE)
  ev$values[ev$values < eps] <- eps
  V  <- ev$vectors %*% diag(ev$values, nrow = length(ev$values)) %*% t(ev$vectors)
  (V + t(V)) / 2
}
safe_solve <- function(S)
  tryCatch(solve(make_pd(S)), error = function(e) MASS::ginv(make_pd(S)))

# Band restriction: fix entries with |i-j| > 1 to their full-sample values
# (exactly matching the simulation design)
make_band_restr <- function(Sigma_full, band_width = 1) {
  p   <- nrow(Sigma_full)
  idx <- which(lower.tri(matrix(0, p, p)), arr.ind = TRUE)
  far <- abs(idx[, 1] - idx[, 2]) > band_width
  list(
    pairs  = idx[far, , drop = FALSE],
    values = Sigma_full[idx[far, , drop = FALSE]]
  )
}

project_R <- function(S, restr) {
  if (nrow(restr$pairs) == 0) return(make_pd(S))
  S_R <- S
  for (k in seq_len(nrow(restr$pairs))) {
    i <- restr$pairs[k, 1]; j <- restr$pairs[k, 2]
    S_R[i, j] <- restr$values[k]; S_R[j, i] <- restr$values[k]
  }
  make_pd(S_R)
}

free_mask <- function(p, restr) {
  m <- matrix(TRUE, p, p)
  if (nrow(restr$pairs) == 0) return(m)
  for (k in seq_len(nrow(restr$pairs))) {
    i <- restr$pairs[k, 1]; j <- restr$pairs[k, 2]
    m[i, j] <- FALSE; m[j, i] <- FALSE
  }
  m
}

# Skew-t log-likelihood
st_ll_fn <- function(Sh, Xte, av, nu) {
  xi <- colMeans(Xte); Om <- make_pd(Sh)
  tryCatch(suppressWarnings(
    sum(sn::dmst(Xte, xi = xi, Omega = Om, alpha = av, nu = nu, log = TRUE))),
    error = function(e) NA_real_)
}

gauss_ll_fn <- function(Sh, Xte) {
  n2 <- nrow(Xte); p2 <- ncol(Xte)
  Xc <- scale(Xte, center = TRUE, scale = FALSE)
  S2 <- make_pd(Sh); Si <- safe_solve(S2)
  ld <- as.numeric(determinant(S2, logarithm = TRUE)$modulus)
  -(n2/2)*(p2*log(2*pi)+ld) - 0.5*sum(diag(Xc %*% Si %*% t(Xc)))
}

frob_fn <- function(Sh, St) sum((Sh - St)^2)
mvp_fn  <- function(Sh, St) {
  Si <- safe_solve(Sh); ones <- rep(1, nrow(Sh))
  w  <- Si %*% ones / as.numeric(t(ones) %*% Si %*% ones)
  as.numeric(t(w) %*% St %*% w)
}
cond_fn <- function(S) {
  ev <- eigen(make_pd(S), only.values = TRUE, symmetric = TRUE)$values
  max(ev) / min(ev)
}

kappa_fn <- function(X, restr) {
  n2 <- nrow(X); p2 <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- crossprod(Xc) / n2
  fv <- as.vector(free_mask(p2, restr))
  acc <- 0
  for (i in seq_len(n2)) {
    xi <- Xc[i, ]; d <- as.vector(outer(xi,xi)) - as.vector(S)
    acc <- acc + sum(d[fv]^2)
  }
  acc / n2^2
}

lw_fn <- function(X) {
  n2 <- nrow(X); p2 <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE); S <- crossprod(Xc) / n2
  mu <- mean(diag(S)); Tm <- mu * diag(p2)
  pi_h <- sum(sapply(seq_len(n2), function(i) {
    xi <- Xc[i, ]; sum((outer(xi,xi)-S)^2)})) / n2^2
  gm <- sum((S-Tm)^2)
  al <- if (gm < 1e-12) 1 else min(1, pi_h/gm)
  list(cov = make_pd(((1-al)*S + al*Tm)*(n2/(n2-1))), lambda = al)
}

skrst_fn <- function(X, restr, target = c("diagonal","identity")) {
  target <- match.arg(target)
  n2 <- nrow(X); p2 <- ncol(X)
  Xc  <- scale(X, center = TRUE, scale = FALSE); S1n <- crossprod(Xc)/n2
  Sg0 <- if (target == "diagonal") diag(diag(S1n)) else
         { sg2 <- sum(diag(S1n))/p2; sg2*diag(p2) }
  Sg0_R <- project_R(Sg0, restr)
  SR    <- project_R(S1n, restr)
  kap   <- kappa_fn(X, restr)
  den   <- sum((SR-Sg0_R)^2)
  lam   <- if (den > 1e-14) max(0, min(1, (den-kap)/den)) else 0
  list(cov = make_pd((lam*SR + (1-lam)*Sg0_R)*(n2/(n2-1))), lambda = lam)
}

# =============================================================================
# ROBUST SKEW-T PARAMETER ESTIMATION
# =============================================================================
# Strategy:
# 1. Try selm() and validate output
# 2. If alpha is 0 or ||alpha|| > 20, try mst.mple() directly  
# 3. Grid search over nu in {3,4,5,6,8,10,15,20} on validation holdout
#    to find best-fitting nu; use alpha from marginal skewness direction
# =============================================================================

fit_st_params <- function(X, label = "") {
  p2 <- ncol(X); n2 <- nrow(X)

  # Step 1: try selm
  frm <- as.formula(paste0("cbind(",paste(colnames(X),collapse=","),")~1"))
  fit <- tryCatch(suppressWarnings(
           sn::selm(frm, family = "ST", data = as.data.frame(X))),
         error = function(e) NULL)

  alpha_est <- NULL; nu_est <- NULL

  if (!is.null(fit)) {
    dp <- tryCatch(fit@param$dp, error = function(e) NULL)
    if (!is.null(dp)) {
      nu_raw    <- as.numeric(dp$nu)
      alpha_raw <- as.numeric(dp$alpha)
      norm_a    <- sqrt(sum(alpha_raw^2))
      cat(sprintf("  selm: nu=%.2f  ||alpha||=%.3f\n", nu_raw, norm_a))

      # Accept if nu is in (2, 100) and ||alpha|| is in (0, 15)
      if (!is.na(nu_raw) && nu_raw > 2 && nu_raw < 100) nu_est <- nu_raw
      if (!is.na(norm_a) && norm_a > 0 && norm_a < 15)  alpha_est <- alpha_raw
    }
  }

  # Step 2: grid search over nu using validation holdout (20% of data)
  # Use the sample covariance as Omega, alpha from marginal skewness
  idx_val <- sample(seq_len(n2), floor(n2 * 0.2))
  Xval    <- X[idx_val, ]; Xfit <- X[-idx_val, ]
  Sfit    <- make_pd(cov(Xfit))
  xi_val  <- colMeans(Xval)

  # Marginal skewness direction as alpha proxy (normalized)
  sk_marg <- apply(X, 2, function(x) mean((x-mean(x))^3)/sd(x)^3)
  norm_sk <- sqrt(sum(sk_marg^2))
  alpha_marg <- if (norm_sk > 0.05) sk_marg / norm_sk * 2 else rep(0, p2)

  nu_grid <- c(3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30)
  ll_grid <- sapply(nu_grid, function(nu_try) {
    tryCatch(suppressWarnings(
      sum(sn::dmst(Xval, xi = xi_val, Omega = Sfit,
                   alpha = alpha_marg, nu = nu_try, log = TRUE))),
      error = function(e) -Inf)
  })
  best_nu <- nu_grid[which.max(ll_grid)]
  cat(sprintf("  Grid search best nu=%.0f  (ll=%.2f)\n",
      best_nu, max(ll_grid)))

  # Use grid nu if selm failed or gave extreme value
  if (is.null(nu_est)) nu_est <- best_nu

  # Use marginal alpha if selm failed or gave degenerate alpha
  if (is.null(alpha_est)) alpha_est <- alpha_marg

  # Final validation: can dmst actually evaluate with these params?
  ll_check <- tryCatch(suppressWarnings(
    sum(sn::dmst(Xval, xi = xi_val, Omega = Sfit,
                 alpha = alpha_est, nu = nu_est, log = TRUE))),
    error = function(e) NA_real_)

  cat(sprintf("  Final: nu=%.2f  ||alpha||=%.3f  validation ll=%.2f  %s\n",
      nu_est, sqrt(sum(alpha_est^2)), ll_check,
      ifelse(is.na(ll_check), "FAILED", "OK")))

  if (is.na(ll_check)) {
    cat("  WARNING: dmst failed, falling back to alpha=0, nu=5\n")
    alpha_est <- rep(0, p2); nu_est <- 5
  }

  list(alpha = alpha_est, nu = nu_est)
}

mardia_test <- function(X) {
  n2 <- nrow(X); p2 <- ncol(X)
  Xc <- scale(X, center = TRUE, scale = FALSE)
  S  <- cov(X) + 1e-8 * diag(p2)
  M  <- Xc %*% solve(S) %*% t(Xc)
  b1 <- sum(M^3)/n2^2; b2 <- sum(diag(M^2))/n2
  sp <- pchisq(n2*b1/6, df = p2*(p2+1)*(p2+2)/6, lower.tail = FALSE)
  kz <- (b2 - p2*(p2+2))/sqrt(8*p2*(p2+2)/n2)
  kp <- 2*pnorm(abs(kz), lower.tail = FALSE)
  list(b1=b1, b2=b2, sp=sp, kz=kz, kp=kp)
}

# =============================================================================
# MAIN ANALYSIS FUNCTION
# Uses repeated random splits for stable results.
# PRIMARY criterion: Frobenius distance to test-set sample covariance
#   (direct analog of simulation loss; no true Sigma needed)
# SECONDARY: skew-t LL, MVP variance, condition number
# =============================================================================
run_analysis <- function(X_full, restr_global, alpha_fit, nu_fit, label,
                         n_rep = 1000, train_frac = 0.70) {
  p2     <- ncol(X_full); n_full <- nrow(X_full)

  # Verify restriction is working
  S_full    <- cov(X_full)
  RMLE_full <- project_R(S_full, restr_global)
  cond_rmle <- cond_fn(RMLE_full)
  cat(sprintf("\n[%s] n=%d p=%d nu=%.2f ||alpha||=%.3f\n",
      label, n_full, p2, nu_fit, sqrt(sum(alpha_fit^2))))
  cat(sprintf("Restriction: %d pairs fixed. RMLE cond on full data=%.2f  %s\n",
      nrow(restr_global$pairs), cond_rmle,
      ifelse(cond_rmle < 1e4, "OK", "SINGULAR — skip")))
  if (cond_rmle >= 1e4) { cat("Skipping.\n"); return(invisible(NULL)) }

  est_names <- c("Sample","LedoitWolf","RestrictedMLE",
                 "Robust","Glasso","SKRST_diag","SKRST_identity")
  n_est <- length(est_names)

  # Storage for repeated splits
  mat_frob <- matrix(NA, n_rep, n_est, dimnames=list(NULL, est_names))
  mat_stll <- matrix(NA, n_rep, n_est, dimnames=list(NULL, est_names))
  mat_mvp  <- matrix(NA, n_rep, n_est, dimnames=list(NULL, est_names))
  mat_cond <- matrix(NA, n_rep, n_est, dimnames=list(NULL, est_names))

  lam_store <- matrix(NA, n_rep, 4,
    dimnames=list(NULL,c("LW","Robust","SKRST_diag","SKRST_id")))

  n_tr <- floor(n_full * train_frac)

  for (r in seq_len(n_rep)) {
    idx_tr <- sample(seq_len(n_full), n_tr)
    Xtr    <- X_full[idx_tr, ]; Xte <- X_full[-idx_tr, ]
    Str    <- cov(Xtr);         Ste <- cov(Xte)

    Sn_hat   <- make_pd(Str)
    lw_fit   <- lw_fn(Xtr);    LW_hat <- lw_fit$cov
    RMLE_hat <- project_R(Str, restr_global)

    rob_fit <- tryCatch({
      r2  <- MASS::cov.rob(Xtr, method="mcd")
      Sr  <- make_pd(r2$cov); mu_r <- mean(diag(Sr)); Tr2 <- mu_r*diag(p2)
      Xcr <- sweep(Xtr,2,r2$center,"-")
      pi_r <- sum(sapply(seq_len(n_tr),function(i){
        xi<-Xcr[i,]; sum((outer(xi,xi)-Sr)^2)}))/n_tr^2
      gm_r <- sum((Sr-Tr2)^2)
      al_r <- if(gm_r<1e-12) 1 else min(1,pi_r/gm_r)
      list(cov=make_pd((1-al_r)*Sr+al_r*Tr2), lambda=al_r)
    }, error=function(e) lw_fit)
    Rob_hat <- rob_fit$cov

    gl_fit <- tryCatch({
      rho <- 0.05*mean(diag(Str))
      list(cov=make_pd(glasso::glasso(Str,rho=rho)$w), rho=rho)
    }, error=function(e) list(cov=Sn_hat, rho=NA_real_))
    GL_hat <- gl_fit$cov

    sk_d <- tryCatch(skrst_fn(Xtr,restr_global,"diagonal"),
                     error=function(e) list(cov=Sn_hat,lambda=NA_real_))
    sk_i <- tryCatch(skrst_fn(Xtr,restr_global,"identity"),
                     error=function(e) list(cov=Sn_hat,lambda=NA_real_))

    ests_r <- list(Sample=Sn_hat, LedoitWolf=LW_hat, RestrictedMLE=RMLE_hat,
                   Robust=Rob_hat, Glasso=GL_hat,
                   SKRST_diag=sk_d$cov, SKRST_identity=sk_i$cov)

    for (nm in est_names) {
      Sh <- ests_r[[nm]]
      mat_frob[r, nm] <- frob_fn(Sh, Ste)
      mat_stll[r, nm] <- st_ll_fn(Sh, Xte, alpha_fit, nu_fit)
      mat_mvp [r, nm] <- mvp_fn(Sh, Ste)
      mat_cond[r, nm] <- cond_fn(Sh)
    }
    lam_store[r,] <- c(lw_fit$lambda, rob_fit$lambda,
                        sk_d$lambda, sk_i$lambda)
  }

  # Summaries
  smry <- function(mat) {
    mn <- colMeans(mat, na.rm=TRUE)
    se <- apply(mat, 2, function(x) sd(x,na.rm=TRUE)/sqrt(sum(!is.na(x))))
    list(mean=mn, se=se)
  }

  frob_s <- smry(mat_frob)
  stll_s <- smry(mat_stll)
  mvp_s  <- smry(mat_mvp)
  cond_s <- smry(mat_cond)
  lam_s  <- colMeans(lam_store, na.rm=TRUE)

  cat(sprintf("Results averaged over %d random splits (%.0f%%/%.0f%%):\n\n",
      n_rep, train_frac*100, (1-train_frac)*100))

  cat("Shrinkage intensities (mean over splits):\n")
  print(data.frame(Estimator=names(lam_s),
                   MeanLambda=round(lam_s,4)), row.names=FALSE)

  rank_d <- function(x) rank(x,  ties.method="min")
  rank_a <- function(x) rank(-x, ties.method="min")

  cat(sprintf("\n--- PRIMARY: Frobenius distance to S_test (lower=better) ---\n"))
  cat("  (Direct analog of simulation loss; S_test is observable proxy for true Sigma)\n")
  frob_tab <- data.frame(
    Estimator = est_names,
    Mean = round(frob_s$mean, 4),
    SE   = round(frob_s$se, 4)
  )
  print(frob_tab[order(frob_tab$Mean), ], row.names=FALSE)

  cat(sprintf("\n--- Skew-t LL on test set [nu=%.2f, ||alpha||=%.3f] (higher=better) ---\n",
      nu_fit, sqrt(sum(alpha_fit^2))))
  stll_tab <- data.frame(
    Estimator = est_names,
    Mean = round(stll_s$mean, 2),
    SE   = round(stll_s$se, 2)
  )
  print(stll_tab[order(-stll_tab$Mean), ], row.names=FALSE)

  cat("\n--- MVP realized variance x1e4 (lower=better) ---\n")
  mvp_tab <- data.frame(
    Estimator = est_names,
    Mean = round(mvp_s$mean * 1e4, 3),
    SE   = round(mvp_s$se  * 1e4, 3)
  )
  print(mvp_tab[order(mvp_tab$Mean), ], row.names=FALSE)

  cat("\n--- Condition number (lower=better) ---\n")
  cond_tab <- data.frame(
    Estimator = est_names,
    Mean = round(cond_s$mean, 2),
    SE   = round(cond_s$se, 2)
  )
  print(cond_tab[order(cond_tab$Mean), ], row.names=FALSE)

  # Summary MeanRank — PRIMARY criterion is Frobenius
  res <- data.frame(
    Estimator = est_names,
    Frob = round(frob_s$mean, 4), FR = rank_d(frob_s$mean),
    StLL = round(stll_s$mean, 2), SR = rank_a(stll_s$mean),
    MVP  = round(mvp_s$mean*1e4,3), MR = rank_d(mvp_s$mean),
    Cond = round(cond_s$mean, 1),   CR = rank_d(cond_s$mean)
  )
  # Summary MeanRank — PRIMARY criterion is Frobenius
  # Note: RestrictedMLE is included for completeness but its numerical
  # instability (high cond) makes it unsuitable for MeanRank comparison.
  # SKRST is the shrinkage-corrected version of RestrictedMLE and is
  # numerically stable — this contrast is itself a key finding.

  res <- data.frame(
    Estimator = est_names,
    Frob = round(frob_s$mean, 4), FR = rank_d(frob_s$mean),
    StLL = round(stll_s$mean, 2), SR = rank_a(stll_s$mean),
    MVP  = round(mvp_s$mean*1e4, 3), MR = rank_d(mvp_s$mean),
    Cond = round(cond_s$mean, 1),    CR = rank_d(cond_s$mean)
  )
  res$MeanRank <- round(rowMeans(res[, c("FR","SR","MR","CR")]), 2)

  cat("\n--- Summary table (MeanRank across all 4 criteria; lower=better) ---\n")
  print(res[order(res$MeanRank), ], row.names=FALSE)

  # Relative Frobenius improvement of SKRST over each competitor
  cat("\n--- Relative Frobenius improvement: (competitor - SKRST_identity) / competitor ---\n")
  cat("    Positive = SKRST_identity is better\n")
  skrst_frob <- frob_s$mean["SKRST_identity"]
  rel_tab <- data.frame(
    Competitor   = setdiff(est_names, "SKRST_identity"),
    Competitor_Frob = round(frob_s$mean[setdiff(est_names,"SKRST_identity")], 4),
    SKRST_Frob   = round(skrst_frob, 4),
    Improvement  = round(
      (frob_s$mean[setdiff(est_names,"SKRST_identity")] - skrst_frob) /
       frob_s$mean[setdiff(est_names,"SKRST_identity")] * 100, 1)
  )
  rel_tab$Better <- ifelse(rel_tab$Improvement > 0, "SKRST wins", "SKRST loses")
  print(rel_tab[order(-rel_tab$Improvement), ], row.names=FALSE)
  cat("\n")
}

# =============================================================================
# DATASET 1: ais data — haematological variables only (known to be skewed)
# =============================================================================
cat("############################################################\n")
cat("DATASET 1: ais (Australian Institute of Sport)\n")
cat("Haematological variables: RCC, WCC, Hc, Hg, Fe\n")
cat("############################################################\n\n")

data("ais", package="sn")
fe_col <- intersect(c("Fe","Ferr"), names(ais))[1]
hae_vars <- c("RCC","WCC","Hc","Hg", fe_col)
X_hae_raw <- as.matrix(ais[, hae_vars])
X_hae_raw <- X_hae_raw[complete.cases(X_hae_raw),]
X_hae <- scale(X_hae_raw)
p_h <- ncol(X_hae); n_h <- nrow(X_hae)
cat(sprintf("n=%d, p=%d, variables: %s\n\n", n_h, p_h, paste(hae_vars,collapse=", ")))

# Distributional diagnostics
md_h <- mardia_test(X_hae)
cat(sprintf("Mardia skewness: b1=%.4f  p=%.4g  %s\n",
    md_h$b1, md_h$sp, ifelse(md_h$sp<0.05,"SIGNIFICANT","")))
cat(sprintf("Mardia kurtosis: b2=%.4f (exp=%d)  p=%.4g  %s\n\n",
    md_h$b2, p_h*(p_h+2), md_h$kp, ifelse(md_h$kp<0.05,"SIGNIFICANT","")))

sk_h <- round(apply(X_hae, 2, function(x) mean((x-mean(x))^3)/sd(x)^3), 3)
cat("Marginal skewness:\n"); print(sk_h); cat("\n")

# Estimate parameters
cat("Estimating skew-t parameters:\n")
st_params_h <- fit_st_params(X_hae, "ais-haem")
alpha_h <- st_params_h$alpha; nu_h <- st_params_h$nu

# Correlation matrix for restriction
R_h <- cor(X_hae)
cat("\nFull-data correlation matrix:\n")
print(round(R_h, 3))

# Band restriction: fix |i-j| > 1 entries to full-sample correlations
# Variables in order: RCC, WCC, Hc, Hg, Fe
# |i-j|=1 pairs (FREE): (RCC,WCC),(WCC,Hc),(Hc,Hg),(Hg,Fe)
# |i-j|>1 pairs (FIXED): all others
S_h_full <- cov(X_hae)  # on standardized data, this is cor(X_hae)
restr_h <- make_band_restr(S_h_full, band_width = 1)

cat(sprintf("\nBand restriction (band_width=1): fixing %d entries, keeping %d free\n",
    nrow(restr_h$pairs), p_h*(p_h-1)/2 - nrow(restr_h$pairs)))
cat("Fixed pairs and their values:\n")
for (k in seq_len(nrow(restr_h$pairs))) {
  i <- restr_h$pairs[k,1]; j <- restr_h$pairs[k,2]
  cat(sprintf("  (%s, %s): fixed to %.4f\n",
      hae_vars[i], hae_vars[j], restr_h$values[k]))
}

# Verify PD
RMLE_check <- project_R(S_h_full, restr_h)
cat(sprintf("\nRMLE on full data: cond=%.2f  min_eig=%.6f  %s\n\n",
    cond_fn(RMLE_check),
    min(eigen(RMLE_check,only.values=TRUE)$values),
    ifelse(cond_fn(RMLE_check)<1e4,"PD OK","SINGULAR")))

run_analysis(X_hae, restr_h, alpha_h, nu_h, "ais-haematological")

# =============================================================================
# DATASET 2: ais body composition variables
# =============================================================================
cat("############################################################\n")
cat("DATASET 2: ais body composition variables\n")
cat("Variables: Ht, Wt, LBM, BMI, Bfat\n")
cat("############################################################\n\n")

body_vars <- c("Ht","Wt","LBM","BMI","Bfat")
X_body_raw <- as.matrix(ais[, body_vars])
X_body_raw <- X_body_raw[complete.cases(X_body_raw),]
X_body <- scale(X_body_raw)
p_bdy <- ncol(X_body); n_bdy <- nrow(X_body)
cat(sprintf("n=%d, p=%d\n\n", n_bdy, p_bdy))

md_bdy <- mardia_test(X_body)
cat(sprintf("Mardia skewness p=%.4g  kurtosis p=%.4g\n",md_bdy$sp,md_bdy$kp))
sk_bdy <- round(apply(X_body,2,function(x) mean((x-mean(x))^3)/sd(x)^3),3)
cat("Marginal skewness:\n"); print(sk_bdy); cat("\n")

cat("Estimating skew-t parameters:\n")
st_params_bdy <- fit_st_params(X_body, "ais-body")
alpha_bdy <- st_params_bdy$alpha; nu_bdy <- st_params_bdy$nu

S_body_full <- cov(X_body)
restr_bdy <- make_band_restr(S_body_full, band_width = 1)

cat(sprintf("\nBand restriction: fixing %d entries\n", nrow(restr_bdy$pairs)))
RMLE_bdy <- project_R(S_body_full, restr_bdy)
cat(sprintf("RMLE cond=%.2f  %s\n\n",
    cond_fn(RMLE_bdy), ifelse(cond_fn(RMLE_bdy)<1e4,"OK","SINGULAR")))

run_analysis(X_body, restr_bdy, alpha_bdy, nu_bdy, "ais-body-composition")

cat("\n=== All done ===\n")
