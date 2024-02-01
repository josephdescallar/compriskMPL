#' Predict function for a phcsh_mpl object
#'
#' Predict the baseline hazards, survival functions and cumulative incidence
#' functions for a specific covariate
#'
#' @param object an object inheriting from class phcsh_mpl, representing a
#' fitted cause specific proportional hazards model.
#'
#' @param covs vector of covariates for prediction of interest. If no
#' covariates are specificied then the reference values from object are used
#'
#' @param risk a number which specifies which risk to plot from phcsh_mpl
#' object
#'
#' @param sand use sandwich estimates of standard error
#'
#' @param n.points number of points to use for plot. Deafult is n = 1000.
#'
#' @export predict.phcsh_mpl
predict.phcsh_mpl <- function(object, covs=NULL, risk=1, n.points=1000,
                              sand=FALSE){
  r = risk
  psif <- function(x, bknots, iknots){
    splines2::mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
                      intercept = object$basis.intercept)
  }
  PSIf <- function(x, bknots, iknots){
    splines2::iSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
                      intercept = object$basis.intercept)
  }
  t.points = seq(object$b.knots[1], object$b.knots[2], length.out = n.points)
  pred.psi = psif(t.points, object$b.knots, object$i.knots[[r]])
  if(is.null(covs)){
    pred.h0r = pred.psi %*% object$"theta"[[r]]
  }
  else{
    pred.h0r = (pred.psi %*% object$"theta"[[r]])*as.vector(exp(t(covs) %*%
                                                                  object$beta[[r]]))
  }
  pred.h0r[pred.h0r < 1e-6] <- 1e-6
  if(sand == FALSE){
    VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                       object$theta.index[[r]]]
  }
  else if(sand == TRUE){
    VarCovMat.theta = object$sand[object$theta.index[[r]],
                                  object$theta.index[[r]]]
  }
  pred.h0r.var = diag(pred.psi %*% VarCovMat.theta %*% t(pred.psi))
  if(object$pos.def==1){
    pred.se = sqrt(pred.h0r.var) / (pred.h0r)
    pred.h0r.lower = as.vector(pred.h0r) * exp(-1.96*pred.se)
    pred.h0r.upper = as.vector(pred.h0r) * exp(1.96*pred.se)
  }
  #Baseline survival
  pred.PSI = PSIf(t.points, object$b.knots, object$i.knots[[r]])
  if(is.null(covs)){
    pred.S0r = exp(-(pred.PSI %*% object$"theta"[[r]]))
  }
  else {
    pred.S0r = exp(-(pred.PSI %*% object$"theta"[[r]])*as.vector(exp(t(covs)
                                                                     %*% object$beta[[r]])))
  }
  pred.S0r.logOR = log((1-pred.S0r) / pred.S0r)
  pred.S0r.prime2 = pred.S0r^2
  pred.chaz.var = diag(pred.PSI %*% VarCovMat.theta %*% t(pred.PSI))
  pred.S0r.var = pred.S0r.prime2 * pred.chaz.var
  pred.S0r.logOR.var = ((1/((1-pred.S0r)*pred.S0r))^2)*pred.S0r.var
  if(object$pos.def==1){
    pred.S0r.logOR.lower = pred.S0r.logOR + 1.96*sqrt(pred.S0r.logOR.var)
    pred.S0r.logOR.upper = pred.S0r.logOR - 1.96*sqrt(pred.S0r.logOR.var)
    pred.S0r.lower = exp(-pred.S0r.logOR.lower) / (1 +
                                                     exp(-pred.S0r.logOR.lower))
    pred.S0r.upper = exp(-pred.S0r.logOR.upper) / (1 +
                                                     exp(-pred.S0r.logOR.upper))
  }

  #Cumulative incidence function
  bma = (t.points - object$b.knots[1]) / 2
  bpa = (t.points + object$b.knots[1]) / 2
  t.gq.change = bma*matrix(object$nodes, nrow=length(bma), ncol =
                             object$gq.points, byrow=TRUE) + bpa
  pred.F0r.psi.gq = pred.F0r.PSI.gq = pred.F0r.h0qt.gq = list()
  pred.F0r.H0qt.gq = pred.F0r.S0qt.gq = pred.F0r.Integrand.gq = list()
  pred.F0r.dhdT.gq = pred.F0r.dSdT.gq = pred.F0r.dFdT.gqT = list()
  for(gq in 1:object$gq.points){
    pred.F0r.psi.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                  object$i.knots[[r]])
    pred.F0r.PSI.gq[[gq]] <- PSIf(t.gq.change[,gq], object$b.knots,
                                  object$i.knots[[r]])
    if(is.null(covs)){
      pred.F0r.h0qt.gq[[gq]] <- (psif(t.gq.change[,gq], object$b.knots,
                                      object$i.knots[[r]]) %*% object$"theta"[[r]])
    }
    else{
      pred.F0r.h0qt.gq[[gq]] <- (psif(t.gq.change[,gq], object$b.knots,
                                      object$i.knots[[r]]) %*% object$"theta"[[r]]) *
        (as.vector(exp(t(covs) %*% object$beta[[r]])))
    }
    pred.F0r.H0qt.gq.r = pred.F0r.S0qt.gq.r = list()
    for(q in 1:object$n.risk){
      if(is.null(covs)){
        pred.F0r.H0qt.gq.r[[q]] = PSIf(t.gq.change[,gq], object$b.knots,
                                       object$i.knots[[q]]) %*% object$"theta"[[q]]
      }
      else{
        pred.F0r.H0qt.gq.r[[q]] = PSIf(t.gq.change[,gq], object$b.knots,
                                       object$i.knots[[q]]) %*% object$"theta"[[q]] *
          (as.vector(exp(t(covs) %*% object$beta[[r]])))
      }
      pred.F0r.S0qt.gq.r[[q]] = exp(-pred.F0r.H0qt.gq.r[[q]])
    }

    pred.F0r.S0qt.gq[[gq]] = Reduce("*", pred.F0r.S0qt.gq.r)
    pred.F0r.Integrand.gq[[gq]] <- object$weights[gq]*pred.F0r.h0qt.gq[[gq]] *
      pred.F0r.S0qt.gq[[gq]]
    pred.F0r.dhdT.gq[[gq]] <- pred.F0r.psi.gq[[gq]]
    pred.F0r.dSdT.gq[[gq]] <- -pred.F0r.PSI.gq[[gq]] *
      as.vector(pred.F0r.S0qt.gq[[gq]])
    pred.F0r.dFdT.gqT[[gq]] <- (as.matrix(pred.F0r.dhdT.gq[[gq]] *
                                            as.vector(pred.F0r.S0qt.gq[[gq]]) + as.vector(pred.F0r.h0qt.gq[[gq]]) *
                                            as.matrix(pred.F0r.dSdT.gq[[gq]]))) * object$weights[gq]
  }
  pred.F0r <- bma*Reduce("+", pred.F0r.Integrand.gq)
  pred.F0r.dFdT <- bma*Reduce("+", pred.F0r.dFdT.gqT)
  pred.F0r.r.var <- diag(pred.F0r.dFdT %*% as.matrix(VarCovMat.theta) %*%
                           t(pred.F0r.dFdT))
  pred.F0r.logOR <- log(pred.F0r / (1-pred.F0r + 1e-12) + 1e-12)
  pred.F0r.logOR.var = ((1/((1-pred.F0r)*pred.F0r))^2)*pred.F0r.r.var
  pred.F0r.log <- log(pred.F0r + 1e-12)
  pred.F0r.log.var <- (1/pred.F0r^2)*pred.F0r.r.var
  if(object$pos.def==1){
    pred.F0r.logOR.lower = pred.F0r.logOR - 1.96*sqrt(pred.F0r.logOR.var)
    pred.F0r.logOR.upper = pred.F0r.logOR + 1.96*sqrt(pred.F0r.logOR.var)
    pred.F0r.lower <- exp(pred.F0r.logOR.lower) /
      (1 + exp(pred.F0r.logOR.lower))
    pred.F0r.upper <- exp(pred.F0r.logOR.upper) /
      (1 + exp(pred.F0r.logOR.upper))
  }
  predValues <- data.frame(cbind("t.points"=t.points, "pred.h0r"=pred.h0r,
                                 "pred.h0r.lower"=pred.h0r.lower, "pred.h0r.upper"=pred.h0r.upper,
                                 "pred.S0r"=pred.S0r, "pred.S0r.lower"=pred.S0r.lower,
                                 "pred.S0r.upper"=pred.S0r.upper, "pred.F0r"=pred.F0r,
                                 "pred.F0r.lower"=pred.F0r.lower, "pred.F0r.upper"=pred.F0r.upper))
  colnames(predValues) <- c("t", "pred.h0r", "pred.h0r.lower", "pred.h0r.upper",
                            "pred.S0r", "pred.S0r.lower", "pred.S0r.upper", "pred.F0r", "pred.F0r.lower",
                            "pred.F0r.upper")
  return(predValues)
}
