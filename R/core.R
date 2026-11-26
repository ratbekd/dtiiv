## Estimation machinery: first step, rotation, corrected and clustered statistics,
## critical values. Port of indep_iiv.R / het_lib.R (and their Python version).
##
## First step: residualize y and x on the controls W (with an intercept), take the
## OLS slope b and the residual r, and whiten: xw = xt / sx, rw = r / sr. At angle
## theta the exogenous component and the error are
##     s = cos(theta) xw + sin(theta) rw,   u = cos(theta) rw - sin(theta) xw,
## and the slope is beta(theta) = b + tan(theta) sr / sx. theta = 0 is OLS.

.prepare <- function(y, x, W, cluster) {
  n <- length(y)
  W1 <- cbind(`(Intercept)` = 1, W)
  q0 <- qr(W1)
  if (q0$rank < ncol(W1)) W1 <- W1[, sort(q0$pivot[seq_len(q0$rank)]), drop = FALSE]
  if (ncol(W1) > 1) {
    sdv <- apply(W1[, -1, drop = FALSE], 2, stats::sd)
    W1 <- cbind(1, sweep(sweep(W1[, -1, drop = FALSE], 2, colMeans(W1[, -1, drop = FALSE])), 2, sdv, "/"))
  }
  qW <- qr(W1)
  xt <- as.vector(qr.resid(qW, x)); yt <- as.vector(qr.resid(qW, y))
  if (sum(xt^2) < 1e-12 * max(1, sum(x^2)))
    stop("the regressor has no variation left after partialling out the controls", call. = FALSE)
  b <- sum(xt * yt) / sum(xt^2); r <- yt - xt * b
  sx <- sqrt(mean(xt^2)); sr <- sqrt(mean(r^2))
  clustered <- !is.null(cluster)
  cl <- if (clustered) as.integer(factor(cluster)) else seq_len(n)
  G <- max(cl)
  list(n = n, k = ncol(W1), W1 = W1, A_W = W1 %*% solve(crossprod(W1) / n), qW = qW,
       xt = xt, yt = yt, r = r, b = b, sx = sx, sr = sr, xw = xt / sx, rw = r / sr,
       psi_b = xt * r / mean(xt^2), psi_sx = (xt^2 - sx^2) / (2 * sx), psi_sr = (r^2 - sr^2) / (2 * sr),
       clustered = clustered, cl = cl, G = G)
}

.beta_at  <- function(th, P) P$b + tan(th) * P$sr / P$sx
.theta_of <- function(beta, P) atan((beta - P$b) * P$sx / P$sr)
.QBIG <- 1e100
## reciprocal condition number of a symmetric matrix on the correlation scale
.rc <- function(S) {
  d <- sqrt(diag(S))
  if (any(!is.finite(d)) || any(d <= 0)) return(0)
  tryCatch(rcond(S / outer(d, d)), error = function(e) 0)
}
.su <- function(th, P) list(s = cos(th) * P$xw + sin(th) * P$rw, u = cos(th) * P$rw - sin(th) * P$xw)

## Uncorrected i.i.d. CUE objective (used for point estimates, as in the paper)
.Q_u <- function(th, P, moms) {
  z <- .su(th, P); g <- .gmat(z$s, z$u, moms); m <- colMeans(g)
  S <- crossprod(sweep(g, 2, m)) / P$n
  w <- if (.rc(S) < 1e-12) NULL else tryCatch(solve(S, m), error = function(e) NULL)
  if (is.null(w)) return(.QBIG)                   # singular moment variance (e.g. binary data)
  q <- P$n * sum(m * w)
  if (is.finite(q)) q else .QBIG
}

## First-step corrected contributions (n x m, not centred). With beta non-NULL the
## slope is held fixed and the angle moves with the first step (the angle channel).
.gcorr <- function(th, P, moms, beta = NULL) {
  if (!is.null(beta)) th <- .theta_of(beta, P)
  cs <- cos(th); sn <- sin(th)
  z <- .su(th, P); s <- z$s; u <- z$u
  D <- .gmat(s, u, moms, deriv = TRUE); g <- D$g; gs <- D$gs; gu <- D$gu
  d_xw <- gs * cs - gu * sn; d_rw <- gs * sn + gu * cs
  d_xt <- d_xw / P$sx; d_r <- d_rw / P$sr
  G_px <- crossprod(P$W1, -d_xt + P$b * d_r) / P$n
  G_py <- crossprod(P$W1, -d_r) / P$n
  G_b  <- colMeans(-d_r * P$xt)
  G_sx <- colMeans(-d_xw * P$xw) / P$sx
  G_sr <- colMeans(-d_rw * P$rw) / P$sr
  if (!is.null(beta)) {
    G_th <- colMeans(gs * u - gu * s); c2 <- cs^2
    G_b  <- G_b  + G_th * (-(P$sx / P$sr) * c2)
    G_sx <- G_sx + G_th * (((beta - P$b) / P$sr) * c2)
    G_sr <- G_sr + G_th * (-((beta - P$b) * P$sx / P$sr^2) * c2)
  }
  gc <- g + P$xt * (P$A_W %*% G_px) + P$yt * (P$A_W %*% G_py) +
    outer(P$psi_b, G_b) + outer(P$psi_sx, G_sx) + outer(P$psi_sr, G_sr)
  list(g = g, gc = gc)
}

## Corrected, cluster-robust score statistic. Returns the statistic and nu, the
## Satterthwaite degrees of freedom of the variance estimate (capped at G - 1).
.Q_c <- function(th, P, moms, beta = NULL) {
  C <- .gcorr(th, P, moms, beta)
  m <- colMeans(C$g)
  gc <- rowsum(sweep(C$gc, 2, colMeans(C$gc)), P$cl, reorder = FALSE)
  Om <- crossprod(gc) / P$n
  w <- if (.rc(Om) < 1e-12) NULL else tryCatch(solve(Om, m), error = function(e) NULL)
  if (is.null(w)) return(c(stat = NA_real_, nu = NA_real_))
  zz <- as.vector(gc %*% w)
  nu <- if (sum(zz^4) > 0) min(P$G - 1, 3 * sum(zz^2)^2 / sum(zz^4)) else P$G - 1
  c(stat = P$n * sum(m * w), nu = nu)
}

## Critical value and p-value. "hotelling": G/(G-1) nu m/(nu-m+1) F(m, nu-m+1);
## tends to chi2(m) with many balanced clusters. "chisq": chi2(m).
.crit <- function(m, nu, G, alpha, type) {
  if (type == "chisq") return(stats::qchisq(1 - alpha, m))
  if (is.na(nu) || nu - m + 1 < 1) return(Inf)
  G / (G - 1) * nu * m / (nu - m + 1) * stats::qf(1 - alpha, m, nu - m + 1)
}
.pval <- function(stat, m, nu, G, type) {
  if (is.na(stat) || m < 1) return(NA_real_)
  if (type == "chisq") return(stats::pchisq(stat, m, lower.tail = FALSE))
  if (is.na(nu) || nu - m + 1 < 1) return(NA_real_)
  stats::pf(stat * (G - 1) / G * (nu - m + 1) / (nu * m), m, nu - m + 1, lower.tail = FALSE)
}

## Sandwich standard error of the CUE slope (CUE with the uncorrected weight,
## corrected and clustered meat), as analytic_se of indep_iiv.R.
.analytic_se <- function(th, P, moms) {
  z <- .su(th, P); gi <- .gmat(z$s, z$u, moms)
  Wt <- tryCatch(solve(crossprod(sweep(gi, 2, colMeans(gi))) / P$n), error = function(e) NULL)
  if (is.null(Wt)) return(c(se_theta = NA_real_, se_beta = NA_real_))
  hh <- 1e-6
  gbar <- function(a) { zz <- .su(a, P); colMeans(.gmat(zz$s, zz$u, moms)) }
  Gd <- (gbar(th + hh) - gbar(th - hh)) / (2 * hh)
  den <- sum(Gd * (Wt %*% Gd))
  if (!is.finite(den) || den <= 0) return(c(se_theta = NA_real_, se_beta = NA_real_))
  gc <- .gcorr(th, P, moms)$gc
  psi_th <- -(1 / den) * as.vector(gc %*% (Wt %*% Gd))
  psi_b <- P$psi_b + (1 / cos(th)^2) * (P$sr / P$sx) * psi_th +
    tan(th) * (P$psi_sr / P$sx - P$sr * P$psi_sx / P$sx^2)
  meat <- function(v) { v <- rowsum(v - mean(v), P$cl, reorder = FALSE); sum(v^2) / P$n }
  c(se_theta = sqrt(meat(psi_th) / P$n), se_beta = sqrt(meat(psi_b) / P$n))
}
