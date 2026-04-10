# ----------------------------- #
# 1. Package management
# ----------------------------- #

required_packages <- c("quantmod", "ggplot2", "MASS", "scatterplot3d")

install_if_missing <- function(packages) {
  installed <- rownames(installed.packages())
  missing <- setdiff(packages, installed)
  if (length(missing) > 0) {
    install.packages(missing, dependencies = TRUE)
  }
}

install_if_missing(required_packages)

suppressPackageStartupMessages({
  library(quantmod)
  library(ggplot2)
  library(MASS)
  library(scatterplot3d)
})

# ----------------------------- #
# 2. Configuration
# ----------------------------- #

FROM_DATE <- "2018-01-01"
TO_DATE <- "2023-12-31"
OUTPUT_DIR <- "skewt_candidates"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Parameters for the 3D figure
T_DF <- 5
ELLIPSOID_LEVEL <- 0.95
PLOT_ANGLE <- 42

# Manual tuning for custom y-axis label in scatterplot3d
Y_LABEL_ANGLE <- 30
Y_LABEL_OFFSET_IN <- -1

# ----------------------------- #
# 3. Helper functions
# ----------------------------- #

sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  
  s <- sd(x)
  if (s == 0) return(NA_real_)
  
  m <- mean(x)
  mean((x - m)^3) / s^3
}

bowley_skewness <- function(x) {
  x <- x[is.finite(x)]
  q <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
  denom <- q[3] - q[1]
  
  if (denom == 0) return(NA_real_)
  
  (q[3] + q[1] - 2 * q[2]) / denom
}

compute_summary_stats <- function(x) {
  x <- x[is.finite(x)]
  
  data.frame(
    N = length(x),
    Mean = mean(x),
    Median = median(x),
    SD = sd(x),
    Moment_Skewness = sample_skewness(x),
    Bowley_Skewness = bowley_skewness(x),
    Min = min(x),
    Max = max(x)
  )
}

log_returns <- function(x) {
  na.omit(diff(log(x)))
}

download_financial_data <- function(from_date, to_date) {
  getSymbols(
    "NFLX",
    src = "yahoo",
    from = from_date,
    to = to_date,
    adjust = TRUE,
    auto.assign = TRUE,
    warnings = FALSE
  )
  
  getSymbols(
    c("VIXCLS", "DCOILWTICO"),
    src = "FRED",
    from = from_date,
    to = to_date,
    auto.assign = TRUE
  )
  
  nflx_ret <- log_returns(Ad(NFLX))
  vix_ret <- log_returns(na.omit(VIXCLS))
  wti_ret <- log_returns(na.omit(DCOILWTICO))
  
  colnames(nflx_ret) <- "NFLX"
  colnames(vix_ret) <- "VIX"
  colnames(wti_ret) <- "WTI"
  
  list(
    NFLX = as.numeric(nflx_ret),
    VIX = as.numeric(vix_ret),
    WTI = as.numeric(wti_ret),
    panel = na.omit(merge(NFLX = nflx_ret, VIX = vix_ret, WTI = wti_ret))
  )
}

save_histogram <- function(x, series_name, output_dir) {
  x <- x[is.finite(x)]
  
  mu <- mean(x)
  med <- median(x)
  sig <- sd(x)
  sk <- sample_skewness(x)
  bsk <- bowley_skewness(x)
  
  binwidth <- diff(range(x)) / 40
  
  p <- ggplot(data.frame(Returns = x), aes(x = Returns)) +
    geom_histogram(
      aes(y = after_stat(density)),
      binwidth = binwidth,
      fill = "steelblue",
      alpha = 0.60,
      color = "white"
    ) +
    geom_density(linewidth = 1.2) +
    stat_function(
      fun = dnorm,
      args = list(mean = mu, sd = sig),
      linetype = "dashed",
      linewidth = 1
    ) +
    geom_vline(xintercept = mu, linewidth = 1) +
    geom_vline(xintercept = med, linetype = "dotted", linewidth = 1) +
    labs(
      title = paste(series_name, "Log-Returns"),
      subtitle = sprintf(
        "Mean = %.4f | Median = %.4f | Moment skewness = %.3f | Bowley skewness = %.3f",
        mu, med, sk, bsk
      ),
      x = "Log-Return",
      y = "Density"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16)
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(series_name, "_hist_full.pdf")),
    plot = p,
    width = 8,
    height = 6,
    device = "pdf"
  )
}

trim_central_mass <- function(X, probs = c(0.005, 0.995)) {
  keep <- rep(TRUE, nrow(X))
  
  for (j in seq_len(ncol(X))) {
    qj <- quantile(X[, j], probs = probs, na.rm = TRUE)
    keep <- keep & X[, j] >= qj[1] & X[, j] <= qj[2]
  }
  
  X[keep, , drop = FALSE]
}

get_ellipsoid_lines <- function(mu, Sigma, df = 5, level = 0.95,
                                n_phi = 26, n_theta = 60) {
  p <- length(mu)
  radius <- p * qf(level, df1 = p, df2 = df)
  
  eig <- eigen(Sigma, symmetric = TRUE)
  A <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0) * radius), nrow = p)
  
  phi <- seq(0, pi, length.out = n_phi)
  theta <- seq(0, 2 * pi, length.out = n_theta)
  
  lines_out <- list()
  
  for (ph in phi) {
    u <- cbind(
      sin(ph) * cos(theta),
      sin(ph) * sin(theta),
      cos(ph)
    )
    pts <- sweep(u %*% t(A), 2, mu, "+")
    lines_out[[length(lines_out) + 1]] <- pts
  }
  
  for (th in seq(0, 2 * pi, length.out = 14)) {
    u <- cbind(
      sin(phi) * cos(th),
      sin(phi) * sin(th),
      cos(phi)
    )
    pts <- sweep(u %*% t(A), 2, mu, "+")
    lines_out[[length(lines_out) + 1]] <- pts
  }
  
  lines_out
}

add_projected_kde_contours <- function(s3d, x, y, plane = c("xy", "xz", "yz"),
                                       const, n = 100, levels = 4,
                                       col = "steelblue", lwd = 1.1) {
  plane <- match.arg(plane)
  
  kd <- kde2d(x, y, n = n)
  
  z_levels <- pretty(range(kd$z), n = levels + 2)
  z_levels <- z_levels[z_levels > min(kd$z) & z_levels < max(kd$z)]
  if (length(z_levels) > levels) {
    z_levels <- tail(z_levels, levels)
  }
  
  contour_lines <- contourLines(kd$x, kd$y, kd$z, levels = z_levels)
  
  for (cl in contour_lines) {
    if (plane == "xy") {
      xx <- cl$x
      yy <- cl$y
      zz <- rep(const, length(xx))
    } else if (plane == "xz") {
      xx <- cl$x
      yy <- rep(const, length(cl$x))
      zz <- cl$y
    } else {
      xx <- rep(const, length(cl$x))
      yy <- cl$x
      zz <- cl$y
    }
    
    proj <- s3d$xyz.convert(xx, yy, zz)
    lines(proj$x, proj$y, col = col, lwd = lwd)
  }
}

add_custom_y_label <- function(s3d, xlim, ylim, zlim,
                               label = "VIX (standardized)",
                               angle = 30,
                               offset_in = -1,
                               cex = 1) {
  p1 <- s3d$xyz.convert(xlim[2], ylim[1], zlim[1])
  p2 <- s3d$xyz.convert(xlim[2], ylim[2], zlim[1])
  
  x1i <- grconvertX(p1$x, from = "user", to = "inches")
  y1i <- grconvertY(p1$y, from = "user", to = "inches")
  x2i <- grconvertX(p2$x, from = "user", to = "inches")
  y2i <- grconvertY(p2$y, from = "user", to = "inches")
  
  xm <- (p1$x + p2$x) / 2
  ym <- (p1$y + p2$y) / 2
  
  dxi <- x2i - x1i
  dyi <- y2i - y1i
  leni <- sqrt(dxi^2 + dyi^2)
  
  nxi <- -dyi / leni
  nyi <- dxi / leni
  
  x_off <- grconvertX(x1i + nxi * offset_in, from = "inches", to = "user") -
    grconvertX(x1i, from = "inches", to = "user")
  
  y_off <- grconvertY(y1i + nyi * offset_in, from = "inches", to = "user") -
    grconvertY(y1i, from = "inches", to = "user")
  
  text(
    x = xm + x_off,
    y = ym + y_off,
    labels = label,
    srt = angle,
    xpd = NA,
    adj = c(0.5, 0.5),
    cex = cex
  )
}

save_trivariate_3d_plot <- function(X_panel, output_dir,
                                    df = 5,
                                    ellipsoid_level = 0.95,
                                    plot_angle = 42,
                                    y_label_angle = 30,
                                    y_label_offset_in = -1) {
  X <- data.frame(coredata(X_panel))
  colnames(X) <- c("NFLX", "VIX", "WTI")
  
  X_trim <- trim_central_mass(X, probs = c(0.005, 0.995))
  X_std <- as.data.frame(scale(X_trim))
  colnames(X_std) <- c("NFLX", "VIX", "WTI")
  
  fit_t <- cov.trob(X_std, nu = df, maxit = 100)
  mu_hat <- fit_t$center
  Sigma_hat <- fit_t$cov
  
  xlim <- range(X_std$NFLX)
  ylim <- range(X_std$VIX)
  zlim <- range(X_std$WTI)
  
  x_wall <- xlim[1]
  y_wall <- ylim[2]
  z_wall <- zlim[1]
  
  pdf(
    file = file.path(output_dir, "figure_trivariate_3d_clean_contours.pdf"),
    width = 8.8,
    height = 7.2
  )
  
  par(mar = c(2.5, 2.5, 3.2, 1))
  
  s3d <- scatterplot3d(
    x = X_std$NFLX,
    y = X_std$VIX,
    z = X_std$WTI,
    pch = 16,
    cex.symbols = 0.42,
    color = adjustcolor("black", alpha.f = 0.35),
    angle = plot_angle,
    scale.y = 1,
    grid = TRUE,
    box = FALSE,
    xlab = "NFLX (standardized)",
    ylab = "",
    zlab = "WTI (standardized)",
    main = "Joint Returns with Fitted Multivariate-t Ellipsoid",
    xlim = xlim,
    ylim = ylim,
    zlim = zlim
  )
  
  add_custom_y_label(
    s3d = s3d,
    xlim = xlim,
    ylim = ylim,
    zlim = zlim,
    label = "VIX (standardized)",
    angle = y_label_angle,
    offset_in = y_label_offset_in,
    cex = 1
  )
  
  # Shadow projections
  floor_proj <- s3d$xyz.convert(X_std$NFLX, X_std$VIX, rep(z_wall, nrow(X_std)))
  points(
    floor_proj$x, floor_proj$y,
    pch = 16, cex = 0.18,
    col = adjustcolor("steelblue", alpha.f = 0.15)
  )
  
  back_proj <- s3d$xyz.convert(X_std$NFLX, rep(y_wall, nrow(X_std)), X_std$WTI)
  points(
    back_proj$x, back_proj$y,
    pch = 16, cex = 0.18,
    col = adjustcolor("firebrick", alpha.f = 0.12)
  )
  
  side_proj <- s3d$xyz.convert(rep(x_wall, nrow(X_std)), X_std$VIX, X_std$WTI)
  points(
    side_proj$x, side_proj$y,
    pch = 16, cex = 0.18,
    col = adjustcolor("darkgreen", alpha.f = 0.12)
  )
  
  # KDE contours on walls
  add_projected_kde_contours(
    s3d, X_std$NFLX, X_std$VIX,
    plane = "xy", const = z_wall,
    n = 120, levels = 4,
    col = "steelblue4", lwd = 1.2
  )
  
  add_projected_kde_contours(
    s3d, X_std$NFLX, X_std$WTI,
    plane = "xz", const = y_wall,
    n = 120, levels = 4,
    col = "firebrick4", lwd = 1.2
  )
  
  add_projected_kde_contours(
    s3d, X_std$VIX, X_std$WTI,
    plane = "yz", const = x_wall,
    n = 120, levels = 4,
    col = "darkgreen", lwd = 1.2
  )
  
  # Fitted multivariate-t ellipsoid
  ellipsoid_lines <- get_ellipsoid_lines(
    mu = mu_hat,
    Sigma = Sigma_hat,
    df = df,
    level = ellipsoid_level
  )
  
  for (pts in ellipsoid_lines) {
    proj <- s3d$xyz.convert(pts[, 1], pts[, 2], pts[, 3])
    lines(proj$x, proj$y, col = adjustcolor("black", alpha.f = 0.8), lwd = 0.9)
  }
  
  legend(
    "topleft",
    inset = c(0.01, 0.01),
    legend = c(
      "3D point cloud",
      "NFLX-VIX floor projection + contours",
      "NFLX-WTI back-wall projection + contours",
      "VIX-WTI side-wall projection + contours",
      "95% multivariate-t ellipsoid"
    ),
    col = c(
      adjustcolor("black", alpha.f = 0.55),
      "steelblue4",
      "firebrick4",
      "darkgreen",
      "black"
    ),
    pch = c(16, 16, 16, 16, NA),
    pt.cex = c(0.7, 0.7, 0.7, 0.7, NA),
    lty = c(NA, 1, 1, 1, 1),
    lwd = c(NA, 1.2, 1.2, 1.2, 0.9),
    bty = "n",
    cex = 0.82,
    y.intersp = 1.05,
    x.intersp = 0.8
  )
  
  dev.off()
}

# ----------------------------- #
# 4. Main workflow
# ----------------------------- #

data_list <- download_financial_data(FROM_DATE, TO_DATE)

series_list <- list(
  NFLX = data_list$NFLX,
  VIX = data_list$VIX,
  WTI = data_list$WTI
)

compare_table <- do.call(
  rbind,
  lapply(names(series_list), function(name) {
    out <- compute_summary_stats(series_list[[name]])
    out$Series <- name
    out
  })
)

compare_table <- compare_table[, c(
  "Series", "N", "Mean", "Median", "SD",
  "Moment_Skewness", "Bowley_Skewness", "Min", "Max"
)]

compare_table$Abs_Moment_Skewness <- abs(compare_table$Moment_Skewness)
compare_table <- compare_table[order(-compare_table$Abs_Moment_Skewness), ]

print(compare_table)
write.csv(
  compare_table,
  file = file.path(OUTPUT_DIR, "skewness_comparison.csv"),
  row.names = FALSE
)

save_histogram(series_list$NFLX, "NFLX", OUTPUT_DIR)
save_histogram(series_list$VIX, "VIX", OUTPUT_DIR)
save_histogram(series_list$WTI, "WTI", OUTPUT_DIR)

save_trivariate_3d_plot(
  X_panel = data_list$panel,
  output_dir = OUTPUT_DIR,
  df = T_DF,
  ellipsoid_level = ELLIPSOID_LEVEL,
  plot_angle = PLOT_ANGLE,
  y_label_angle = Y_LABEL_ANGLE,
  y_label_offset_in = Y_LABEL_OFFSET_IN
)
