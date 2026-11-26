.fmt_set <- function(cs, digits = 4) {
  if (!nrow(cs)) return("empty")
  paste(vapply(seq_len(nrow(cs)), function(i)
    sprintf("%s%s, %s%s", if (cs$lower_at_edge[i]) "(" else "[",
            formatC(cs$lower[i], digits = digits, format = "g"),
            formatC(cs$upper[i], digits = digits, format = "g"),
            if (cs$upper_at_edge[i]) ")" else "]"), ""), collapse = " U ")
}

.set_table <- function(x, digits = 4) {
  data.frame(
    set = names(x$sets),
    m = vapply(x$sets, `[[`, 0L, "m"),
    estimate = vapply(x$sets, `[[`, 0, "estimate"),
    se = vapply(x$sets, `[[`, 0, "se"),
    J_p = vapply(x$sets, `[[`, 0, "J_p"),
    confidence_set = vapply(x$sets, function(s) .fmt_set(s$confset, digits), ""),
    row.names = NULL, stringsAsFactors = FALSE)
}

#' @export
print.dtiiv <- function(x, digits = 4, ...) {
  cat("DT-IIV: ", deparse(x$formula), "\n", sep = "")
  cat(sprintf("n = %d%s, %d controls (incl. intercept)%s\n", x$n,
              if (x$clustered) sprintf(" in %d clusters", x$G) else "", x$k,
              if (x$n_dropped) sprintf(", %d rows dropped for missing values", x$n_dropped) else ""))
  cat(sprintf("OLS %s (se %s)%s\n\n", formatC(x$ols[1], digits = digits, format = "g"),
              formatC(x$ols[2], digits = digits, format = "g"),
              if (is.null(x$tsls)) "" else sprintf(";  2SLS %s (se %s)",
                formatC(x$tsls[1], digits = digits, format = "g"), formatC(x$tsls[2], digits = digits, format = "g"))))
  tb <- .set_table(x, digits)
  tb$estimate <- formatC(tb$estimate, digits = digits, format = "g")
  tb$se <- ifelse(is.na(tb$se), "--", formatC(tb$se, digits = digits, format = "g"))
  tb$J_p <- ifelse(is.na(tb$J_p), "--", formatC(tb$J_p, digits = 3, format = "f"))
  names(tb) <- c("Set", "m", "CUE", "SE", "J p", sprintf("%g%% confidence set", 100 * x$level))
  print(tb, row.names = FALSE, right = FALSE)
  cat(sprintf("\nCritical values: %s.  ( ) marks an end at the edge of the search range.\n",
              if (x$crit == "hotelling") "Hotelling-type F with cluster degrees of freedom" else "chi-squared"))
  invisible(x)
}

#' Summarize a DT-IIV fit
#'
#' @param object A `dtiiv` fit.
#' @param x A `summary.dtiiv` object.
#' @param digits Significant digits.
#' @param ... Not used.
#' @return An object of class `summary.dtiiv`, printed with the shape
#'   diagnostics, the roots of the variance tendency and, for each moment set,
#'   the CUE estimate, the J test and the confidence set.
#' @export
summary.dtiiv <- function(object, ...) {
  structure(list(fit = object, table = .set_table(object)), class = "summary.dtiiv")
}

#' @rdname summary.dtiiv
#' @export
print.summary.dtiiv <- function(x, digits = 4, ...) {
  f <- x$fit
  print(f, digits = digits)
  cat("\nShape of the data (identification needs non-Gaussian components):\n")
  sh <- f$diagnostics$shape; sh$skewness <- round(sh$skewness, 3); sh$excess_kurtosis <- round(sh$excess_kurtosis, 3)
  print(sh, row.names = FALSE)
  th <- f$diagnostics$third; th$estimate <- round(th$estimate, 4); th$t <- round(th$t, 2)
  cat("Third moments of the whitened regressor (x) and OLS residual (r): mu_jk = E[x^j r^k]\n")
  print(th, row.names = FALSE)
  cat("\nRoots of the variance tendency E[u^2 s] = 0:\n")
  vr <- f$vt_roots
  if (!nrow(vr)) cat("  none\n") else {
    vr$beta <- signif(vr$beta, digits); vr$theta <- round(vr$theta, 4); vr$rho_xu <- round(vr$rho_xu, 3)
    vr$screen_stat <- round(vr$screen_stat, 2); vr$screen_p <- signif(vr$screen_p, 3)
    print(vr, row.names = FALSE)
    cat("  in_domain: |theta| <= pi/4.  screen: test of E[u^3 s], E[u s^3], E[u^2 s^2] - 1 at the root (3 df).\n")
  }
  cat("\nMoment sets:\n")
  for (s in f$sets) {
    cat(sprintf("  %-10s %s\n", s$name, paste(s$moments, collapse = ", ")))
    ## other local minima of the CUE objective that would themselves pass the J test
    oth <- s$minima[-1, , drop = FALSE]
    oth <- oth[oth$Q <= stats::qchisq(f$level, max(s$m - 1, 1)) & !oth$at_edge, , drop = FALSE]
    if (s$m == 1) oth <- oth[oth$Q < 1e-6, , drop = FALSE]
    if (nrow(oth))
      cat(sprintf("             other admissible CUE minima at beta = %s\n",
                  paste(signif(oth$beta, digits), collapse = ", ")))
    if (isTRUE(s$at_edge)) cat("             the CUE minimum is at the edge of the search range\n")
  }
  invisible(x)
}

#' Extract DT-IIV estimates
#'
#' @param object A `dtiiv` fit.
#' @param set Name of a moment set; default the first. For `coef`, `"all"`
#'   returns the CUE of every set.
#' @param parm Not used.
#' @param level Must equal the level of the fit; refit to change it.
#' @param ... Not used.
#' @return `coef`: the CUE slope (named vector). `confint`: a matrix with one
#'   row per interval of the confidence set (the set can have several pieces
#'   or none).
#' @export
coef.dtiiv <- function(object, set = names(object$sets)[1], ...) {
  if (identical(set, "all")) return(vapply(object$sets, `[[`, 0, "estimate"))
  .chk_set(object, set)
  stats::setNames(object$sets[[set]]$estimate, all.vars(object$formula)[2])
}

#' @rdname coef.dtiiv
#' @export
confint.dtiiv <- function(object, parm = NULL, level = object$level, set = names(object$sets)[1], ...) {
  if (!isTRUE(all.equal(level, object$level)))
    stop("the confidence sets are computed at level ", object$level, "; refit with dtiiv(level = ...)", call. = FALSE)
  .chk_set(object, set)
  cs <- object$sets[[set]]$confset
  out <- as.matrix(cs[, c("lower", "upper")])
  colnames(out) <- c("lower", "upper")
  rownames(out) <- if (nrow(out)) paste0("piece", seq_len(nrow(out))) else NULL
  out
}

.chk_set <- function(object, set) {
  if (!is.character(set) || length(set) != 1 || !set %in% names(object$sets))
    stop("`set` must be one of: ", paste(names(object$sets), collapse = ", "), call. = FALSE)
}

#' Plot the test statistic against the slope
#'
#' Plots, for each moment set, the corrected score statistic divided by its
#' critical value. The confidence set is where the curve lies below one.
#'
#' @param x A `dtiiv` fit.
#' @param sets Moment sets to plot; default all.
#' @param log Use a log scale for the ratio.
#' @param ... Passed to [graphics::plot()].
#' @return `x`, invisibly.
#' @export
plot.dtiiv <- function(x, sets = names(x$sets), log = TRUE, ...) {
  pr <- lapply(x$sets[sets], function(s) s$profile)
  rat <- lapply(pr, function(p) pmax(p$stat / p$crit, 1e-3))
  yl <- range(unlist(rat), 1, finite = TRUE)
  cols <- grDevices::hcl.colors(length(sets), "Dark 3")
  graphics::plot(pr[[1]]$beta, rat[[1]], type = "n", log = if (log) "y" else "", ylim = yl,
                 xlab = expression(beta), ylab = "statistic / critical value", ...)
  for (i in seq_along(sets)) graphics::lines(pr[[i]]$beta, rat[[i]], col = cols[i], lwd = 2)
  graphics::abline(h = 1, lty = 2)
  graphics::abline(v = x$ols[1], lty = 3, col = "grey40")
  if (!is.null(x$tsls)) graphics::abline(v = x$tsls[1], lty = 3, col = "grey60")
  graphics::legend("topright", legend = sets, col = cols, lwd = 2, bty = "n", cex = 0.8)
  invisible(x)
}

#' Simulate data for DT-IIV examples
#'
#' Generates \eqn{y = 2 + \beta x + 0.3 w + u_0} and
#' \eqn{x = 1 + 0.5 w + s_0 + \lambda u_0}, with \eqn{s_0}, \eqn{u_0} and
#' \eqn{w} independent and standardized, so that OLS is biased whenever
#' \eqn{\lambda \neq 0}.
#'
#' @param n Sample size.
#' @param beta True slope.
#' @param lambda Feedback of the error into the regressor.
#' @param s_dist,u_dist Distributions of \eqn{s_0} and \eqn{u_0}:
#'   `"chisq"` (skewed, \eqn{\chi^2_4}), `"t"` (heavy tails, \eqn{t_5}),
#'   `"lognormal"` (skewed and heavy-tailed) or `"normal"`.
#' @param clusters Number of equal-sized clusters; `NULL` for none.
#' @return A data frame with `y`, `x`, `w` and, with clusters, `cluster`.
#' @examples
#' d <- dtiiv_simulate(1000)
#' summary(lm(y ~ x + w, data = d))$coefficients["x", ]   # biased OLS
#' @export
dtiiv_simulate <- function(n = 2000, beta = 1, lambda = 0.4, s_dist = "chisq", u_dist = "lognormal",
                           clusters = NULL) {
  draw <- function(d) {
    v <- switch(d,
                chisq = stats::rchisq(n, 4),
                t = stats::rt(n, 5),
                lognormal = exp(stats::rnorm(n, 0, 0.5)),
                normal = stats::rnorm(n),
                stop("unknown distribution: ", d, call. = FALSE))
    (v - mean(v)) / stats::sd(v)
  }
  s0 <- draw(s_dist); u0 <- draw(u_dist); w <- stats::rnorm(n)
  x <- 1 + 0.5 * w + s0 + lambda * u0
  y <- 2 + beta * x + 0.3 * w + u0
  d <- data.frame(y = y, x = x, w = w)
  if (!is.null(clusters)) d$cluster <- rep(seq_len(clusters), length.out = n)
  d
}
