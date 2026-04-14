#!/usr/bin/env Rscript
# =============================================================================
# Plots for Section 5 (Application) of the SKRST paper
# Requires: ggplot2, ggcorrplot (or corrplot), reshape2
# Run AFTER SKRST_real_data_analysis.R has been sourced so that
# X_hae, X_body, restr_h, restr_bdy, alpha_h, nu_h, alpha_bdy, nu_bdy
# are available in the workspace.
# =============================================================================

library(ggplot2)
library(reshape2)

# Install ggcorrplot if needed
if (!requireNamespace("ggcorrplot", quietly = TRUE))
  install.packages("ggcorrplot", repos = "https://cloud.r-project.org")
library(ggcorrplot)

theme_set(theme_bw(base_size = 11))
pal_est <- c(
  "SKRST"          = "#E41A1C",   # red
  "RestrictedMLE"  = "#FF7F00",   # orange
  "LedoitWolf"     = "#377EB8",   # blue
  "Sample"         = "#4DAF4A",   # green
  "Glasso"         = "#984EA3",   # purple
  "Robust"         = "#A65628"    # brown
)

# ── helper: run 50 splits and collect per-split losses ──────────────────────
collect_splits <- function(X_full, restr_global, alpha_fit, nu_fit,
                           n_rep = 50, train_frac = 0.70) {
  p2   <- ncol(X_full); n_full <- nrow(X_full)
  n_tr <- floor(n_full * train_frac)

  est_names <- c("SKRST","RestrictedMLE","LedoitWolf","Sample","Glasso","Robust")

  rows <- vector("list", n_rep * length(est_names))
  idx  <- 1L

  for (r in seq_len(n_rep)) {
    id_tr <- sample(seq_len(n_full), n_tr)
    Xtr   <- X_full[id_tr, ]; Xte <- X_full[-id_tr, ]
    Str   <- cov(Xtr);        Ste <- cov(Xte)

    make_pd_loc <- function(S, eps = 1e-8) {
      S <- (S+t(S))/2; ev <- eigen(S, symmetric=TRUE)
      ev$values[ev$values < eps] <- eps
      V <- ev$vectors %*% diag(ev$values,nrow=length(ev$values)) %*% t(ev$vectors)
      (V+t(V))/2
    }
    safe_solve_loc <- function(S)
      tryCatch(solve(make_pd_loc(S)), error=function(e) MASS::ginv(make_pd_loc(S)))
    project_loc <- function(S, restr) {
      S_R <- S
      for (k in seq_len(nrow(restr$pairs))) {
        i<-restr$pairs[k,1]; j<-restr$pairs[k,2]
        S_R[i,j]<-restr$values[k]; S_R[j,i]<-restr$values[k]
      }
      make_pd_loc(S_R)
    }
    kappa_loc <- function(X, restr) {
      n2<-nrow(X); p3<-ncol(X)
      Xc<-scale(X,center=TRUE,scale=FALSE); S<-crossprod(Xc)/n2
      fm<-matrix(TRUE,p3,p3)
      for(k in seq_len(nrow(restr$pairs))){i<-restr$pairs[k,1];j<-restr$pairs[k,2];fm[i,j]<-FALSE;fm[j,i]<-FALSE}
      fv<-as.vector(fm); acc<-0
      for(i in seq_len(n2)){xi<-Xc[i,];d<-as.vector(outer(xi,xi))-as.vector(S);acc<-acc+sum(d[fv]^2)}
      acc/n2^2
    }

    # Sample
    Sn <- make_pd_loc(Str)

    # LedoitWolf
    n2<-nrow(Xtr); p3<-ncol(Xtr)
    Xc<-scale(Xtr,center=TRUE,scale=FALSE); S<-crossprod(Xc)/n2
    mu<-mean(diag(S)); Tm<-mu*diag(p3)
    pi_h<-sum(sapply(seq_len(n2),function(i){xi<-Xc[i,];sum((outer(xi,xi)-S)^2)}))/n2^2
    gm<-sum((S-Tm)^2); al<-if(gm<1e-12) 1 else min(1,pi_h/gm)
    LW <- make_pd_loc(((1-al)*S+al*Tm)*(n2/(n2-1)))

    # RestrictedMLE
    RM <- project_loc(Str, restr_global)

    # Robust
    Rob <- tryCatch({
      r2<-MASS::cov.rob(Xtr,method="mcd"); Sr<-make_pd_loc(r2$cov)
      mu_r<-mean(diag(Sr)); Tr2<-mu_r*diag(p3)
      Xcr<-sweep(Xtr,2,r2$center,"-")
      pi_r<-sum(sapply(seq_len(n_tr),function(i){xi<-Xcr[i,];sum((outer(xi,xi)-Sr)^2)}))/n_tr^2
      gm_r<-sum((Sr-Tr2)^2); al_r<-if(gm_r<1e-12) 1 else min(1,pi_r/gm_r)
      make_pd_loc((1-al_r)*Sr+al_r*Tr2)
    }, error=function(e) LW)

    # Glasso
    GL <- tryCatch(
      make_pd_loc(glasso::glasso(Str,rho=0.05*mean(diag(Str)))$w),
      error=function(e) Sn)

    # SKRST (diagonal target)
    SK <- tryCatch({
      S1n<-crossprod(scale(Xtr,center=TRUE,scale=FALSE))/n_tr
      Sg0<-diag(diag(S1n)); Sg0_R<-project_loc(Sg0,restr_global)
      SR<-project_loc(S1n,restr_global); kap<-kappa_loc(Xtr,restr_global)
      den<-sum((SR-Sg0_R)^2); lam<-if(den>1e-14) max(0,min(1,(den-kap)/den)) else 0
      make_pd_loc((lam*SR+(1-lam)*Sg0_R)*(n_tr/(n_tr-1)))
    }, error=function(e) Sn)

    ests_r <- list(SKRST=SK, RestrictedMLE=RM, LedoitWolf=LW,
                   Sample=Sn, Glasso=GL, Robust=Rob)

    for (nm in est_names) {
      Sh <- ests_r[[nm]]
      frob <- sum((Sh - Ste)^2)
      si   <- safe_solve_loc(Sh); ones<-rep(1,p3)
      w    <- si%*%ones/as.numeric(t(ones)%*%si%*%ones)
      mvp  <- as.numeric(t(w)%*%Ste%*%w)
      ev   <- eigen(make_pd_loc(Sh),only.values=TRUE,symmetric=TRUE)$values
      cond <- max(ev)/min(ev)
      stll <- tryCatch(suppressWarnings(
        sum(sn::dmst(Xte,xi=colMeans(Xte),Omega=make_pd_loc(Sh),
                     alpha=alpha_fit,nu=nu_fit,log=TRUE))),
        error=function(e) NA_real_)
      rows[[idx]] <- data.frame(rep=r, Estimator=nm,
        Frob=frob, MVP=mvp*1e4, Cond=cond, StLL=stll,
        stringsAsFactors=FALSE)
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

# ── Re-run data loading (in case run standalone) ─────────────────────────────
library(sn); library(MASS); library(glasso)
set.seed(2025)
data("ais", package="sn")
fe_col    <- intersect(c("Fe","Ferr"), names(ais))[1]
hae_vars  <- c("RCC","WCC","Hc","Hg",fe_col)
body_vars <- c("Ht","Wt","LBM","BMI","Bfat")

X_hae  <- scale(as.matrix(ais[complete.cases(ais[,hae_vars]),  hae_vars]))
X_body <- scale(as.matrix(ais[complete.cases(ais[,body_vars]), body_vars]))

make_band_restr_plot <- function(S, band_width=1) {
  p<-nrow(S); idx<-which(lower.tri(matrix(0,p,p)),arr.ind=TRUE)
  far<-abs(idx[,1]-idx[,2])>band_width
  list(pairs=idx[far,,drop=FALSE], values=S[idx[far,,drop=FALSE]])
}
make_pd_g <- function(S,eps=1e-8){S<-(S+t(S))/2;ev<-eigen(S,symmetric=TRUE);ev$values[ev$values<eps]<-eps;V<-ev$vectors%*%diag(ev$values,nrow=length(ev$values))%*%t(ev$vectors);(V+t(V))/2}

restr_h   <- make_band_restr_plot(cov(X_hae))
restr_bdy <- make_band_restr_plot(cov(X_body))

# Skew-t params (use grid search result directly)
nu_h    <- 13.49; alpha_h   <- apply(X_hae,  2, function(x) mean((x-mean(x))^3)/sd(x)^3)
nu_bdy  <-  3.76; alpha_bdy <- apply(X_body, 2, function(x) mean((x-mean(x))^3)/sd(x)^3)
norm_h   <- sqrt(sum(alpha_h^2));   alpha_h   <- if(norm_h>0.05)   alpha_h/norm_h*2   else rep(0,5)
norm_bdy <- sqrt(sum(alpha_bdy^2)); alpha_bdy <- if(norm_bdy>0.05) alpha_bdy/norm_bdy*5.442 else rep(0,5)

cat("Collecting split results for Dataset 1 (haematological)...\n")
df_hae  <- collect_splits(X_hae,  restr_h,   alpha_h,   nu_h)
cat("Collecting split results for Dataset 2 (body composition)...\n")
df_body <- collect_splits(X_body, restr_bdy, alpha_bdy, nu_bdy)

df_hae$Dataset  <- "Haematological\n(RCC, WCC, Hc, Hg, Fe)"
df_body$Dataset <- "Body Composition\n(Ht, Wt, LBM, BMI, Bfat)"
df_all <- rbind(df_hae, df_body)

# Fix factor ordering: SKRST first
df_all$Estimator <- factor(df_all$Estimator,
  levels = c("SKRST","RestrictedMLE","LedoitWolf","Sample","Glasso","Robust"))

# =============================================================================
# PLOT 1: Boxplot of Frobenius distance across 50 splits
# =============================================================================
# Cap Robust outliers for body composition for visibility
df_frob <- df_all
df_frob$Frob_capped <- pmin(df_frob$Frob, 3)   # cap at 3 for display

p1 <- ggplot(df_frob, aes(x = Estimator, y = Frob_capped, fill = Estimator)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.5, width = 0.6) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_fill_manual(values = pal_est, guide = "none") +
  labs(
    title = "Frobenius Distance to Test-Set Sample Covariance",
    subtitle = "Averaged over 50 random 70/30 splits; lower is better",
    x = NULL, y = expression(paste("||", hat(Sigma)[train] - S[test], "||"[F]^2))
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        strip.text  = element_text(size = 10),
        plot.title  = element_text(size = 12, face = "bold"))

ggsave("fig_frob_boxplot.pdf", p1,
       width = 8, height = 4.5, device = cairo_pdf)
# ggsave("fig_frob_boxplot.png", p1,
       # width = 8, height = 4.5, dpi = 200)
cat("Saved fig_frob_boxplot\n")

# =============================================================================
# PLOT 2: Mean Frobenius with error bars (SE), both datasets side by side
# =============================================================================
smry_frob <- aggregate(Frob ~ Estimator + Dataset, df_all,
  FUN = function(x) c(mean=mean(x), se=sd(x)/sqrt(length(x))))
smry_frob <- do.call(data.frame, smry_frob)
names(smry_frob)[3:4] <- c("Mean","SE")

# Remove Robust from body (outlier distorts scale)
smry_frob_plot <- smry_frob[!(smry_frob$Estimator=="Robust" &
                               grepl("Body",smry_frob$Dataset)),]

p2 <- ggplot(smry_frob_plot,
             aes(x = Estimator, y = Mean, ymin = Mean-SE, ymax = Mean+SE,
                 colour = Estimator, shape = Estimator)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.25) +
  facet_wrap(~Dataset, scales = "free_y") +
  scale_colour_manual(values = pal_est, guide = "none") +
  scale_shape_manual(values = c(16,17,15,18,8,4), guide = "none") +
  labs(
    title = "Mean Frobenius Distance to Test-Set Covariance ± SE",
    # subtitle = "50 random splits; Robust excluded from Body Composition panel (scale)",
    x = NULL,
    y = expression(paste("Mean ||", hat(Sigma)[train] - S[test], "||"[F]^2))
  ) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 9),
        strip.text  = element_text(size = 10))

ggsave("fig_frob_errbar.pdf", p2,
       width = 8, height = 4.5, device = cairo_pdf)
# ggsave("/mnt/user-data/outputs/fig_frob_errbar.png", p2,
       # width = 8, height = 4.5, dpi = 200)
cat("Saved fig_frob_errbar\n")

# =============================================================================
# PLOT 3: Correlation matrix heatmaps with band structure overlay
# (shows which entries are free vs. restricted)
# =============================================================================
plot_corrmat <- function(R, restr, var_names, title_str) {
  p <- nrow(R)
  df_r <- melt(R); names(df_r) <- c("Var1","Var2","value")
  df_r$Var1 <- factor(df_r$Var1, levels=var_names)
  df_r$Var2 <- factor(df_r$Var2, levels=rev(var_names))

  # Mark restricted entries
  restricted_pairs <- restr$pairs
  df_r$restricted <- FALSE
  for (k in seq_len(nrow(restricted_pairs))) {
    i <- restricted_pairs[k,1]; j <- restricted_pairs[k,2]
    df_r$restricted[df_r$Var1==var_names[i] & df_r$Var2==rev(var_names)[p+1-j]] <- TRUE
    df_r$restricted[df_r$Var1==var_names[j] & df_r$Var2==rev(var_names)[p+1-i]] <- TRUE
  }

  ggplot(df_r, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile(colour="white", size=0.5) +
    geom_tile(data=df_r[df_r$restricted,],
              aes(x=Var1, y=Var2), fill=NA,
              colour="black", size=1.2, linetype="dashed") +
    geom_text(aes(label=sprintf("%.2f", value)), size=3) +
    scale_fill_gradient2(low="#2166AC", mid="white", high="#D6604D",
                         midpoint=0, limits=c(-1,1),
                         name="Correlation") +
    labs(title=title_str,
         # subtitle="Dashed border = restricted entry (fixed to full-sample value)",
         x=NULL, y=NULL) +
    theme(axis.text=element_text(size=10),
          plot.title=element_text(size=11, face="bold"),
          legend.position="right")
}

R_hae  <- cor(X_hae);  colnames(R_hae) <- rownames(R_hae) <- hae_vars
R_body <- cor(X_body); colnames(R_body) <- rownames(R_body) <- body_vars

p3a <- plot_corrmat(R_hae,  restr_h,   hae_vars,
                    "Haematological Variables: Sample Correlation Matrix")
p3b <- plot_corrmat(R_body, restr_bdy, body_vars,
                    "Body Composition Variables: Sample Correlation Matrix")

ggsave("fig_corrmat_haem.pdf",  p3a, width=6, height=5, device=cairo_pdf)
# ggsave("/mnt/user-data/outputs/fig_corrmat_haem.png",  p3a, width=6, height=5, dpi=200)
ggsave("fig_corrmat_body.pdf",  p3b, width=6, height=5, device=cairo_pdf)
# ggsave("/mnt/user-data/outputs/fig_corrmat_body.png",  p3b, width=6, height=5, dpi=200)
cat("Saved correlation matrix figures\n")

# =============================================================================
# PLOT 4: Marginal density plots showing skewness (key motivation)
# =============================================================================
df_hae_long  <- melt(as.data.frame(X_hae),  variable.name="Variable", value.name="Value")
df_body_long <- melt(as.data.frame(X_body), variable.name="Variable", value.name="Value")
df_hae_long$Dataset  <- "Haematological"
df_body_long$Dataset <- "Body Composition"
df_dens <- rbind(df_hae_long, df_body_long)

p4 <- ggplot(df_dens, aes(x=Value)) +
  geom_histogram(aes(y=after_stat(density)), bins=20,
                 fill="steelblue", colour="white", alpha=0.7) +
  geom_density(colour="#E41A1C", linewidth=0.8) +
  facet_wrap(Dataset ~ Variable, scales="free", ncol=5) +
  labs(title="Marginal Distributions of AIS Variables (Standardized)",
       # subtitle="Red curve: kernel density estimate. Skewness and heavy tails visible.",
       x="Standardised value", y="Density") +
  theme(strip.text=element_text(size=8),
        axis.text=element_text(size=7),
        plot.title=element_text(size=11, face="bold"))

ggsave("fig_marginal_densities.pdf", p4,
       width=10, height=5.5, device=cairo_pdf)
# ggsave("/mnt/user-data/outputs/fig_marginal_densities.png", p4,
       # width=10, height=5.5, dpi=180)
cat("Saved fig_marginal_densities\n")

# =============================================================================
# PLOT 5: Relative Frobenius improvement of SKRST over each competitor
# =============================================================================
smry2 <- aggregate(Frob ~ Estimator + Dataset, df_all, mean)
skrst_vals <- smry2[smry2$Estimator == "SKRST", ]
smry2 <- merge(smry2, skrst_vals[,c("Dataset","Frob")], by="Dataset", suffixes=c("","_SKRST"))
smry2$RelImprovement <- (smry2$Frob - smry2$Frob_SKRST) / smry2$Frob * 100
smry2 <- smry2[smry2$Estimator != "SKRST", ]

# Cap Robust for display
smry2$RelImprovement_capped <- pmin(smry2$RelImprovement, 100)

p5 <- ggplot(smry2, aes(x=reorder(Estimator, -RelImprovement_capped),
                         y=RelImprovement_capped, fill=Estimator)) +
  geom_col(width=0.6) +
  geom_hline(yintercept=0, linetype="dashed") +
  facet_wrap(~Dataset) +
  scale_fill_manual(values=pal_est, guide="none") +
  labs(
    title="Frobenius Improvement of SKRST Over Competitors",
    # subtitle=expression(paste("(Competitor Frob - SKRST Frob) / Competitor Frob × 100%;  positive = SKRST better")),
    x=NULL, y="Relative Improvement (%)"
  ) +
  theme(axis.text.x=element_text(angle=35, hjust=1, size=9),
        strip.text=element_text(size=10),
        plot.title=element_text(size=11, face="bold"))

ggsave("fig_relative_improvement.pdf", p5,
       width=8, height=4.5, device=cairo_pdf)
# ggsave("/mnt/user-data/outputs/fig_relative_improvement.png", p5,
       # width=8, height=4.5, dpi=200)
cat("Saved fig_relative_improvement\n")

cat("\nAll plots saved to /mnt/user-data/outputs/\n")
