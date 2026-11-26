## Moment functions g(s, u) and their derivatives with respect to s and u.
##
## At a tested slope the residual and the regressor are whitened and rotated, which
## gives the error u and the exogenous component s (see dtiiv()). Each moment is a
## function g(s, u) whose mean is zero at the true slope under one of the
## restrictions of the paper. The derivatives gs = dg/ds and gu = dg/du enter the
## first-step correction of the score.

#' Define a moment function
#'
#' Creates a moment for use in [dtiiv()]. A moment is a function \eqn{g(s, u)}
#' whose mean is zero at the true slope, together with its derivatives with
#' respect to \eqn{s} and \eqn{u}. The derivatives are needed for the
#' first-step correction of the test statistic.
#'
#' @param id Short identifier, used in output.
#' @param g,gs,gu Functions of `(s, u)` returning a numeric vector: the moment,
#'   its derivative with respect to `s` and its derivative with respect to `u`.
#' @param label Text shown in printed output.
#' @param restriction The restriction under which the moment has mean zero
#'   (for documentation only).
#' @return An object of class `dtiiv_moment`.
#' @examples
#' # E[u * s^5] = 0 under mean independence
#' m5 <- dtiiv_moment("us5", g = function(s, u) u * s^5,
#'                    gs = function(s, u) 5 * u * s^4,
#'                    gu = function(s, u) s^5,
#'                    label = "E[u s^5]", restriction = "MI")
#' @export
dtiiv_moment <- function(id, g, gs, gu, label = id, restriction = NA_character_) {
  stopifnot(is.character(id), length(id) == 1, is.function(g), is.function(gs), is.function(gu))
  structure(list(id = id, label = label, restriction = restriction, g = g, gs = gs, gu = gu),
            class = "dtiiv_moment")
}

#' @export
print.dtiiv_moment <- function(x, ...) {
  cat(sprintf("<dtiiv moment> %s: %s  (%s)\n", x$id, x$label,
              if (is.na(x$restriction)) "user-defined" else x$restriction))
  invisible(x)
}

phi0 <- function(z) exp(-z^2 / 2)
## robust standardization of s: median and 1.4826 * MAD (rob_lib.R). The location and
## scale need no first-step correction: under MI, E[u k(s)] = 0 for every k.
.zs <- function(s) {
  m <- stats::median(s); tau <- 1.4826 * stats::median(abs(s - m))
  if (!is.finite(tau) || tau <= 0) tau <- stats::sd(s)
  list(z = (s - m) / tau, tau = tau)
}

.moment_registry <- function() {
  M <- list(
    dtiiv_moment("vt", function(s, u) u^2 * s, function(s, u) u^2, function(s, u) 2 * u * s,
                 "E[u^2 s]  (variance tendency)", "VT"),
    dtiiv_moment("us2", function(s, u) u * s^2, function(s, u) 2 * u * s, function(s, u) s^2,
                 "E[u s^2]", "MI"),
    dtiiv_moment("u3s", function(s, u) u^3 * s, function(s, u) u^3, function(s, u) 3 * u^2 * s,
                 "E[u^3 s]", "RMI"),
    dtiiv_moment("us3", function(s, u) u * s^3, function(s, u) 3 * u * s^2, function(s, u) s^3,
                 "E[u s^3]", "MI"),
    dtiiv_moment("u2s2", function(s, u) u^2 * s^2 - 1, function(s, u) 2 * u^2 * s, function(s, u) 2 * u * s^2,
                 "E[u^2 s^2] - 1", "CH"),
    dtiiv_moment("u2m1s3", function(s, u) (u^2 - 1) * s^3, function(s, u) 3 * (u^2 - 1) * s^2,
                 function(s, u) 2 * u * s^3, "E[(u^2 - 1) s^3]", "CH"),
    dtiiv_moment("utanhs", function(s, u) u * tanh(s), function(s, u) u * (1 - tanh(s)^2),
                 function(s, u) tanh(s), "E[u tanh(s)]", "MI"),
    dtiiv_moment("uphis", function(s, u) u * phi0(s), function(s, u) -u * s * phi0(s),
                 function(s, u) phi0(s), "E[u exp(-s^2/2)]", "MI"),
    dtiiv_moment("utanhz", function(s, u) u * tanh(.zs(s)$z),
                 function(s, u) { Z <- .zs(s); u * (1 - tanh(Z$z)^2) / Z$tau },
                 function(s, u) tanh(.zs(s)$z), "E[u tanh(z)], z = robust-standardized s", "MI"),
    dtiiv_moment("uphiz", function(s, u) u * phi0(.zs(s)$z),
                 function(s, u) { Z <- .zs(s); -u * Z$z * phi0(Z$z) / Z$tau },
                 function(s, u) phi0(.zs(s)$z), "E[u exp(-z^2/2)], z = robust-standardized s", "MI"),
    dtiiv_moment("stanhu", function(s, u) s * tanh(u), function(s, u) tanh(u),
                 function(s, u) s * (1 - tanh(u)^2), "E[s tanh(u)]", "RMI"),
    dtiiv_moment("stanhu2", function(s, u) s * tanh(u / 2), function(s, u) tanh(u / 2),
                 function(s, u) s * (1 - tanh(u / 2)^2) / 2, "E[s tanh(u/2)]", "RMI"),
    dtiiv_moment("sphiu", function(s, u) s * phi0(u), function(s, u) phi0(u),
                 function(s, u) -s * u * phi0(u), "E[s exp(-u^2/2)]", "RMI"),
    dtiiv_moment("suphiu", function(s, u) s * u * phi0(u), function(s, u) u * phi0(u),
                 function(s, u) s * (1 - u^2) * phi0(u), "E[s u exp(-u^2/2)]", "RMI"),
    dtiiv_moment("su3", function(s, u) s * u^3, function(s, u) u^3, function(s, u) 3 * s * u^2,
                 "E[s u^3]", "RMI")
  )
  stats::setNames(M, vapply(M, `[[`, "", "id"))
}

.set_registry <- list(
  "VT"        = "vt",
  "two-sided" = c("vt", "us2"),
  "one-sided" = c("us2", "us3"),
  "four"      = c("vt", "us2", "u3s", "us3"),
  "five"      = c("vt", "us2", "u3s", "us3", "u2s2"),
  "variance"  = c("vt", "us2", "u2s2", "u2m1s3"),
  "bounded"   = c("utanhz", "uphiz"),
  "mean"      = c("us2", "us3", "utanhs", "uphis"),
  "error"     = c("stanhu", "stanhu2", "sphiu", "suphiu"),
  "both"      = c("us2", "us3", "utanhs", "uphis", "stanhu", "stanhu2", "sphiu", "suphiu")
)

#' Moments and moment sets available in dtiiv
#'
#' Lists the built-in moments and the named moment sets that [dtiiv()]
#' accepts in its `sets` argument.
#'
#' The restrictions are those of the paper. VT: the variance of the error does
#' not depend on the exogenous component (variance tendency). MI: mean
#' independence, \eqn{E[u_0 \mid s_0] = 0}. RMI: reverse mean independence,
#' \eqn{E[s_0 \mid u_0] = 0}. CH: MI with a constant conditional variance.
#' Independence implies all of them.
#'
#' The named sets are:
#' * `"VT"`: the variance tendency alone (exactly identified);
#' * `"two-sided"`, `"one-sided"`, `"four"`, `"five"`, `"variance"`: the
#'   polynomial sets of the paper;
#' * `"bounded"`: the bounded mean tendencies of the paper,
#'   \eqn{E[u \tanh z]} and \eqn{E[u e^{-z^2/2}]} with \eqn{z} the
#'   standardization of \eqn{s} by its median and 1.4826 times its median
#'   absolute deviation; valid under MI and recommended for heavy-tailed
#'   regressors, with the default (Hotelling) critical values;
#' * `"mean"`: mean tendencies \eqn{E[u h(s)]} with
#'   \eqn{h(s) = (s^2, s^3, \tanh s, e^{-s^2/2})}, valid under MI;
#' * `"error"`: error-side tendencies \eqn{E[s a(u)]} with
#'   \eqn{a(u) = (\tanh u, \tanh(u/2), \phi(u), u\phi(u))}, valid under RMI;
#' * `"both"`: the eight mean and error-side moments together.
#'
#' @return A list with elements `moments` (a data frame) and `sets` (a named
#'   list of moment identifiers).
#' @examples
#' dtiiv_moments()$sets$mean
#' @export
dtiiv_moments <- function() {
  M <- .moment_registry()
  list(moments = data.frame(id = names(M), moment = vapply(M, `[[`, "", "label"),
                            restriction = vapply(M, `[[`, "", "restriction"), row.names = NULL),
       sets = .set_registry)
}

## Resolve the `sets` argument of dtiiv() into a named list of lists of moment objects.
.resolve_sets <- function(sets) {
  reg <- .moment_registry()
  one <- function(el, nm) {
    if (inherits(el, "dtiiv_moment")) el <- list(el)
    if (is.character(el)) {
      if (length(el) == 1 && el %in% names(.set_registry)) el <- .set_registry[[el]]
      bad <- setdiff(el, names(reg))
      if (length(bad)) stop("unknown moment(s) in set '", nm, "': ", paste(bad, collapse = ", "),
                            ". See dtiiv_moments().", call. = FALSE)
      return(unname(reg[el]))
    }
    if (is.list(el)) {
      out <- lapply(el, function(e) {
        if (inherits(e, "dtiiv_moment")) e
        else if (is.character(e) && length(e) == 1 && e %in% names(reg)) reg[[e]]
        else stop("set '", nm, "' must contain moment ids or dtiiv_moment() objects", call. = FALSE)
      })
      return(out)
    }
    stop("cannot interpret set '", nm, "'", call. = FALSE)
  }
  if (is.character(sets)) {
    bad <- setdiff(sets, names(.set_registry))
    if (length(bad)) stop("unknown set(s): ", paste(bad, collapse = ", "),
                          ". Use one of ", paste(names(.set_registry), collapse = ", "),
                          ", or a named list of moment ids.", call. = FALSE)
    sets <- stats::setNames(as.list(sets), sets)
  }
  if (is.null(names(sets)) || any(names(sets) == ""))
    stop("`sets` must be a character vector of set names or a named list", call. = FALSE)
  out <- Map(one, sets, names(sets))
  for (nm in names(out)) {
    ids <- vapply(out[[nm]], `[[`, "", "id")
    if (anyDuplicated(ids)) stop("set '", nm, "' repeats a moment", call. = FALSE)
  }
  out
}

## n x m matrices of g, gs, gu for a list of moments
.gmat <- function(s, u, moms, deriv = FALSE) {
  g <- vapply(moms, function(m) m$g(s, u), numeric(length(s)))
  if (!is.matrix(g)) g <- matrix(g, ncol = length(moms))
  if (!deriv) return(g)
  gs <- vapply(moms, function(m) m$gs(s, u), numeric(length(s)))
  gu <- vapply(moms, function(m) m$gu(s, u), numeric(length(s)))
  if (!is.matrix(gs)) { gs <- matrix(gs, ncol = length(moms)); gu <- matrix(gu, ncol = length(moms)) }
  list(g = g, gs = gs, gu = gu)
}
