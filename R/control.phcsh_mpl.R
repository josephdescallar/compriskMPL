#' Ancillary arguments for controling phcsh_mpl fits
#'
#' This is used to set various numeric parameters controlling a Cox model fit
#' using phcsh_mpl. Typically it would only be used in a call to phcsh_mpl.
#'
#' @param knots.perc.limit a vector which defines the percentile range of data
#' for the boundary knot and internal knots for mspline estimation
#'
#' @param gq.points defines the number of nodes to used for gaussian quadrature
#' estimation of numerical integrals.
#'
#' @param dgr number of degrees to be used in mspline estimation
#'
#' @param basis.intercept indicates whether intercept is to be included for
#' mspline estimation
#'
#' @param max.outer maximum number of iterations to perform for the marginal
#' likelihood method for penalty terms estimation.
#'
#' @param max.iter maximum number of iterations to perform for the quasi-Newton,
#' multiplicative iterative algorithm
#'
#' @param lambda initial values for the penalty term estimation
#'
#' @param n.basis vector which contains the number of basis functions to use by
#' risk
#'
#' @param iter.disp indicator which displays iteration count during model fit
#'
#' @param iknots.pos define the percentile placement of internal knots for
#' mspline estimation
#'
#' @param aps boolean which defines whether the automatic parameter selection
#' is used for penalty term selection. Default is TRUE.
#'
#' @param tmid boolean which indicates whether midpoint of interval censored
#' observations is used for mspline estimation
#'
#' @param inner.conv convergence criteria for quasi Newton - Multiplicative
#' iterative algorithm for model fit.
#'
#' @export
control.phcsh_mpl <- function(knots.perc.limit = c(0.05,0.95), gq.points = 15,
                              dgr = 2, basis.intercept = FALSE, max.outer = 10L,
                              max.iter = 3000L, lambda=NULL, n.basis = NULL,
                              iter.disp = FALSE, iknots.pos = NULL, aps = TRUE,
                              tmid = TRUE, inner.conv = 1e-8){

  out = list(knots.perc.limit = knots.perc.limit, gq.points = gq.points,
             dgr=dgr, basis.intercept = basis.intercept,
             max.outer = max.outer, max.iter = max.iter,
             lambda = lambda, n.basis = n.basis, iter.disp = iter.disp,
             iknots.pos = iknots.pos, aps = aps, tmid = tmid,
             inner.conv = inner.conv)
  out
}
