#' DT-IIV estimation and inference without an external instrument
#'
#' Estimates the slope \eqn{\beta} of an endogenous regressor \eqn{x} in
#' \deqn{y = \beta x + W'\delta + u}
#' using distributional tendency (DT) conditions in place of an external
#' instrument. Write the regressor, after partialling out the controls, as
#' \eqn{x = \lambda u + s_0}, where \eqn{s_0} is its exogenous component. At a
#' tested slope \eqn{b} the residual \eqn{y - bx} and the part \eqn{s} of
#' \eqn{x} orthogonal to it are whitened, so that second moments carry no
#' information about \eqn{b}. The DT conditions are higher-moment restrictions
#' that hold only at the true slope, for example \eqn{E[u h(s)] = 0} under mean
#' independence. They need non-Gaussian data: skewness or excess kurtosis in
#' \eqn{s_0} (mean tendencies), in \eqn{u_0} (error-side tendencies) or a
#' skewed error (variance tendency).
#'
#' For each moment set the function reports
#' * a confidence set for \eqn{\beta}: the slopes at which a corrected,
#'   cluster-robust score test does not reject. It is robust to weak
#'   identification and can be an interval, a union of intervals or empty. An
#'   empty set means that the linear model with these restrictions is
#'   rejected at every slope;
#' * a continuously updated GMM (CUE) point estimate with a sandwich
#'   standard error, and the J test of the overidentifying restrictions.
#'
#' The search is over the rotation angle \eqn{\theta}, with
#' \eqn{\beta(\theta) = b_{OLS} + \tan\theta\, \sigma_r/\sigma_x}. The default
#' range \eqn{|\theta| \le \pi/4} is the fundamental domain of the paper: the
#' conditions are symmetric in \eqn{s} and \eqn{u}, so each solution has a
#' partner at \eqn{\theta \mp \pi/2}, and the domain picks the solution in
#' which the error is the component closer to the OLS residual.
#'
#' Critical values. With `crit = "hotelling"` (default) the statistic is
#' compared with \eqn{G/(G-1)\, \nu m/(\nu - m + 1)\, F_{m, \nu - m + 1}}, where
#' \eqn{G} is the number of clusters and \eqn{\nu \le G - 1} is a Satterthwaite
#' estimate of the degrees of freedom of the variance. With many balanced
#' clusters this tends to the \eqn{\chi^2_m} cutoff of `crit = "chisq"`.
#' Without clusters each observation is its own cluster.
#'
#' @param formula A two-sided formula `y ~ x` with one endogenous regressor.
#' @param data A data frame.
#' @param controls Optional one-sided formula of exogenous controls, for example
#'   `~ age + I(age^2) + factor(region)`. An intercept is always included.
#' @param cluster Optional one-sided formula naming the cluster variable (for
#'   example `~ state`), or a vector of cluster identifiers of length
#'   `nrow(data)`.
#' @param sets Moment sets to use: names from [dtiiv_moments()], or a named
#'   list whose elements are vectors of moment identifiers or lists of
#'   [dtiiv_moment()] objects.
#' @param level Confidence level of the sets.
#' @param crit Critical values: `"hotelling"` or `"chisq"`.
#' @param theta_range Range of rotation angles searched. The default is the
#'   fundamental domain \eqn{[-\pi/4, \pi/4]}.
#' @param ngrid Number of grid points in `theta_range` for the set and the
#'   CUE search. Crossings of the critical value are refined by root finding.
#' @param instruments Optional one-sided formula of external instruments, used
#'   only to report a 2SLS comparison.
#' @return An object of class `dtiiv`. Use [summary()], [confint()],
#'   [coef()], [plot()] and [dtiiv_test()] on it.
#' @references Dzhumashev, R. and Tursunalieva, A. Identification through
#'   distributional tendencies: intrinsic instruments for linear models.
#' @examples
#' set.seed(1)
#' d <- dtiiv_simulate(n = 2000, beta = 1)
#' fit <- dtiiv(y ~ x, data = d, controls = ~ w, sets = c("VT", "mean"))
#' fit
#' summary(fit)
#' @export
dtiiv <- function(formula, data, controls = NULL, cluster = NULL,
                  sets = c("VT", "two-sided", "four", "bounded", "mean", "error", "both"),
                  level = 0.95, crit = c("hotelling", "chisq"),
                  theta_range = c(-pi / 4, pi / 4), ngrid = 241, instruments = NULL) {
  cl_call <- match.call()
  crit <- match.arg(crit)
  stopifnot(is.data.frame(data), level > 0, level < 1, length(theta_range) == 2,
            theta_range[1] < theta_range[2], all(abs(theta_range) < pi / 2), ngrid >= 11)
  alpha <- 1 - level

  ## ---- data -------------------------------------------------------------------------
  if (!inherits(formula, "formula") || length(formula) != 3)
    stop("`formula` must be of the form y ~ x", call. = FALSE)
  rhs <- attr(stats::terms(formula), "term.labels")
  if (length(rhs) != 1) stop("`formula` must have exactly one regressor on the right: y ~ x. ",
                             "Put exogenous variables in `controls`.", call. = FALSE)
  cl_vec <- NULL; cl_name <- NULL
  if (!is.null(cluster)) {
    if (inherits(cluster, "formula")) {
      cl_name <- all.vars(cluster)
      if (length(cl_name) != 1) stop("`cluster` must name one variable", call. = FALSE)
    } else {
      if (length(cluster) != nrow(data)) stop("`cluster` must have length nrow(data)", call. = FALSE)
      cl_vec <- cluster
    }
  }
  vars <- unique(c(all.vars(formula), if (!is.null(controls)) all.vars(controls),
                   if (!is.null(instruments)) all.vars(instruments), cl_name))
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("variables not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  if (!is.null(cl_vec)) keep <- keep & !is.na(cl_vec)
  dd <- data[keep, , drop = FALSE]
  if (!is.null(cl_vec)) cl_vec <- cl_vec[keep]
  n_drop <- sum(!keep)
  mf <- stats::model.frame(formula, dd)
  y <- as.numeric(stats::model.response(mf))
  x <- stats::model.matrix(formula, mf)[, -1]
  if (is.matrix(x)) stop("the regressor must be numeric (one column)", call. = FALSE)
  x <- as.numeric(x)
  W <- if (is.null(controls)) matrix(nrow = nrow(dd), ncol = 0)
       else stats::model.matrix(controls, dd)[, -1, drop = FALSE]
  cluster_v <- if (!is.null(cl_name)) dd[[cl_name]] else cl_vec
  if (length(unique(y)) <= 2 && length(unique(x)) <= 2)
    stop("both y and x are binary. Every function of a binary variable is linear in it, ",
         "so the tendency conditions hold mechanically at OLS (mean tendencies) or at ",
         "zero (error-side tendencies) and carry no information about the slope. ",
         "DT-IIV needs a non-binary regressor or outcome.", call. = FALSE)
  if (length(unique(x)) <= 2)
    warning("x is binary: the mean tendencies hold mechanically at the OLS slope and cannot ",
            "correct for endogeneity.", call. = FALSE)

  P <- .prepare(y, x, W, cluster_v)
  if (P$clustered && P$G < 10)
    warning("only ", P$G, " clusters: cluster-robust inference is unreliable", call. = FALSE)
  S <- .resolve_sets(sets)

  ## ---- comparators ----------------------------------------------------------------------
  meat1 <- function(v) { v <- rowsum(v - mean(v), P$cl, reorder = FALSE); sum(v^2) }
  ols <- c(estimate = P$b, se = sqrt(meat1(P$xt * P$r)) / sum(P$xt^2))
  tsls <- NULL
  if (!is.null(instruments)) {
    Z <- stats::model.matrix(instruments, dd)[, -1, drop = FALSE]
    Zt <- qr.resid(P$qW, Z)
    xh <- as.vector(Zt %*% qr.coef(qr(Zt), P$xt))
    bi <- sum(xh * P$yt) / sum(xh * P$xt); e <- P$yt - P$xt * bi
    tsls <- c(estimate = bi, se = sqrt(meat1(xh * e)) / abs(sum(xh * P$xt)))
  }

  ## ---- per set ------------------------------------------------------------------------
  th <- seq(theta_range[1], theta_range[2], length.out = ngrid)
  res <- lapply(names(S), function(nm) .fit_set(P, S[[nm]], nm, th, alpha, crit))
  names(res) <- names(S)

  structure(list(
    call = cl_call, formula = formula, n = P$n, n_dropped = n_drop, G = P$G, clustered = P$clustered,
    k = P$k, level = level, crit = crit, theta_range = theta_range, ngrid = ngrid,
    ols = ols, tsls = tsls, sets = res, moments = S,
    vt_roots = .vt_roots(P, alpha, crit),
    diagnostics = .diagnostics(P),
    scale = c(b_ols = P$b, sx = P$sx, sr = P$sr),
    prep = P), class = "dtiiv")
}

## One moment set: S profile on the grid, confidence set, CUE, J test
.fit_set <- function(P, moms, nm, th, alpha, crit) {
  m <- length(moms)
  bg <- .beta_at(th, P)
  Qc <- t(vapply(bg, function(b) .Q_c(NA, P, moms, beta = b), c(stat = 0, nu = 0)))
  cv <- vapply(seq_along(bg), function(i) .crit(m, Qc[i, "nu"], P$G, alpha, crit), 0)
  pv <- vapply(seq_along(bg), function(i) .pval(Qc[i, "stat"], m, Qc[i, "nu"], P$G, crit), 0)
  f <- function(t) {
    q <- .Q_c(NA, P, moms, beta = .beta_at(t, P))
    q[["stat"]] - .crit(m, q[["nu"]], P$G, alpha, crit)
  }
  acc <- !is.na(Qc[, "stat"]) & Qc[, "stat"] <= cv
  pieces <- .components(acc, th, f)
  ci <- if (length(pieces)) do.call(rbind, lapply(pieces, function(p)
    data.frame(lower = .beta_at(p[1], P), upper = .beta_at(p[2], P),
               lower_at_edge = p[3] == 1, upper_at_edge = p[4] == 1)))
        else data.frame(lower = numeric(0), upper = numeric(0),
                        lower_at_edge = logical(0), upper_at_edge = logical(0))

  ## CUE: all local minima of the uncorrected objective, refined
  Qu <- vapply(th, .Q_u, 0, P = P, moms = moms)
  loc <- which(diff(sign(diff(c(Inf, Qu, Inf)))) == 2)
  if (!length(loc)) loc <- which.min(Qu)
  sols <- do.call(rbind, lapply(loc, function(i) {
    o <- stats::optimize(.Q_u, c(th[max(1, i - 1)], th[min(length(th), i + 1)]), P = P, moms = moms)
    c(theta = o$minimum, Q = o$objective)
  }))
  sols <- sols[order(sols[, "Q"]), , drop = FALSE]
  th_hat <- sols[1, "theta"]
  se <- .analytic_se(th_hat, P, moms)
  J <- J_nu <- J_p <- NA_real_
  if (m > 1) {
    lo <- max(min(th), th_hat - 0.05); hi <- min(max(th), th_hat + 0.05)
    oc <- stats::optimize(function(a) { q <- .Q_c(a, P, moms)[[1]]; if (is.finite(q)) q else .QBIG },
                          c(lo, hi))
    q <- .Q_c(oc$minimum, P, moms); J <- q[["stat"]]; J_nu <- q[["nu"]]
    J_p <- .pval(J, m - 1, J_nu, P$G, crit)
  }
  minima <- data.frame(beta = .beta_at(sols[, "theta"], P), theta = sols[, "theta"], Q = sols[, "Q"],
                       at_edge = abs(sols[, "theta"] - min(th)) < 1e-3 | abs(sols[, "theta"] - max(th)) < 1e-3)
  list(name = nm, m = m, moments = vapply(moms, `[[`, "", "id"),
       estimate = .beta_at(th_hat, P), se = unname(se["se_beta"]), theta = th_hat,
       at_edge = minima$at_edge[1], minima = minima,
       J = J, J_df = if (m > 1) m - 1 else NA_integer_, J_nu = J_nu, J_p = J_p,
       confset = ci,
       profile = data.frame(theta = th, beta = bg, stat = Qc[, "stat"], nu = Qc[, "nu"],
                            crit = cv, p = pv, accept = acc))
}

## Connected components of an acceptance vector; interior crossings refined by uniroot.
## Returns a list of c(theta_lo, theta_hi, lo_at_edge, hi_at_edge).
.components <- function(acc, th, f) {
  N <- length(th); out <- list(); i <- 1
  refine <- function(a, b) {
    fa <- f(a); fb <- f(b)
    if (is.finite(fa) && is.finite(fb) && sign(fa) != sign(fb))
      tryCatch(stats::uniroot(f, c(a, b), f.lower = fa, f.upper = fb, tol = 1e-10)$root,
               error = function(e) (a + b) / 2)
    else (a + b) / 2
  }
  while (i <= N) {
    if (acc[i]) {
      j <- i
      while (j < N && acc[j + 1]) j <- j + 1
      lo <- if (i > 1) refine(th[i - 1], th[i]) else th[1]
      hi <- if (j < N) refine(th[j], th[j + 1]) else th[N]
      out[[length(out) + 1]] <- c(lo, hi, as.numeric(i == 1), as.numeric(j == N))
      i <- j + 1
    } else i <- i + 1
  }
  out
}

## Roots of the sample variance tendency mean(u^2 s) on (-pi/2, pi/2), with the
## screen test of the remaining polynomial moments at each root.
.vt_roots <- function(P, alpha, crit) {
  reg <- .moment_registry()
  f <- function(a) { z <- .su(a, P); mean(z$u^2 * z$s) }
  th <- seq(-pi / 2 + 1e-4, pi / 2 - 1e-4, length.out = 4001)
  v <- vapply(th, f, 0)
  ix <- which(diff(sign(v)) != 0)
  if (!length(ix)) return(data.frame(beta = numeric(0), theta = numeric(0), in_domain = logical(0),
                                     rho_xu = numeric(0), screen_stat = numeric(0), screen_p = numeric(0)))
  rt <- vapply(ix, function(i) stats::uniroot(f, c(th[i], th[i + 1]), tol = 1e-12)$root, 0)
  scr <- unname(reg[c("vt", "u3s", "us3", "u2s2")])
  out <- lapply(rt, function(a) {
    q <- .Q_c(a, P, scr)
    b <- .beta_at(a, P); u <- P$yt - b * P$xt
    data.frame(beta = b, theta = a, in_domain = abs(a) <= pi / 4,
               rho_xu = stats::cor(P$xt, u),
               screen_stat = q[["stat"]], screen_p = .pval(q[["stat"]], 3, q[["nu"]], P$G, crit))
  })
  do.call(rbind, out)
}

## Shape of the data: skewness and excess kurtosis of the residualized regressor and
## the OLS residual, and the four third cross-moments with cluster-robust t.
.diagnostics <- function(P) {
  sk <- function(v) mean(v^3) / mean(v^2)^1.5
  ku <- function(v) mean(v^4) / mean(v^2)^2 - 3
  tt <- function(v) {
    m <- mean(v); vv <- rowsum(v - m, P$cl, reorder = FALSE)
    m / sqrt(sum(vv^2) / P$n^2)
  }
  A <- cbind(mu30 = P$xw^3, mu21 = P$xw^2 * P$rw, mu12 = P$xw * P$rw^2, mu03 = P$rw^3)
  list(shape = data.frame(variable = c("x (partialled)", "OLS residual"),
                          skewness = c(sk(P$xw), sk(P$rw)), excess_kurtosis = c(ku(P$xw), ku(P$rw))),
       third = data.frame(moment = colnames(A), estimate = colMeans(A), t = apply(A, 2, tt), row.names = NULL))
}

#' Test a hypothesized slope
#'
#' Computes the corrected score statistic of each moment set at given slopes.
#' The confidence sets of [dtiiv()] are the slopes at which this test does
#' not reject.
#'
#' @param object A `dtiiv` fit.
#' @param beta0 Numeric vector of slopes to test.
#' @param sets Names of moment sets in `object`; default all.
#' @return A data frame with the set, the slope, the statistic, its degrees of
#'   freedom, \eqn{\nu}, the p-value and the critical value at the fit's level.
#' @examples
#' set.seed(1)
#' d <- dtiiv_simulate(n = 2000, beta = 1)
#' fit <- dtiiv(y ~ x, data = d, controls = ~ w, sets = "mean")
#' dtiiv_test(fit, beta0 = c(0, 1))
#' @export
dtiiv_test <- function(object, beta0, sets = names(object$sets)) {
  stopifnot(inherits(object, "dtiiv"), is.numeric(beta0))
  P <- object$prep; alpha <- 1 - object$level
  out <- lapply(sets, function(nm) {
    moms <- object$moments[[nm]]; m <- length(moms)
    do.call(rbind, lapply(beta0, function(b) {
      q <- .Q_c(NA, P, moms, beta = b)
      data.frame(set = nm, beta0 = b, stat = q[["stat"]], df = m, nu = q[["nu"]],
                 p_value = .pval(q[["stat"]], m, q[["nu"]], P$G, object$crit),
                 crit = .crit(m, q[["nu"]], P$G, alpha, object$crit))
    }))
  })
  do.call(rbind, out)
}
