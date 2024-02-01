#' Fit Cause Specific Proportional Hazards Regression Model Via MPL
#'
#' Simultaneously estimate the regression coefficients and the baseline hazard
#' function of proportional hazard Cause specific hazards models for competing
#' risks with or without a cure fraction using maximum penalized likelihood
#' (MPL)
#'
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function (in the Surv, use type = "interval" with
#' interval censored data).
#'
#' @param risk a numeric vector which contains the risk/failure for the
#' i-th observation in the data frame.
#'
#' @param z an optional data frame which contains the covariate to be used when
#' a cure fraction is specified. If left NULL, then no cure fraction is assumed,
#' otherwise parameter estimation will occur for a cure fraction based on
#' the specified variables.
#'
#' @param data a data frame that includes the variables named in the formula
#' argument. If no dataset is indicated, variables will be taken from the global
#' environment.
#'
#' @param control object of class phcsh_mpl_control specifying control options
#' such smoothing parameter value, iteration numbers. Refer to
#' phcsh_mpl_control for defaults
#'
#' @param ... Other arguments. In coxph_mpl, these elements will be passed to
#' coxph_mpl_control. In print.coxph_mpl, these elements will be passed to the
#' print function.
#'
#' @return an object of class phcsh_mpl representing the fit. See
#' phcsh_mpl.object for defaults
#'
#' @export
phcsh_mpl <- function(formula, data, risk, z=NULL, control, ...){
  if(is.null(z)){
    fit <- fit.phcsh_mpl(formula, risk, data, control, ...)
  }
  else fit <- fit.phcshcf_mpl(formula, risk, z, data, control, ...)
  fit
}
