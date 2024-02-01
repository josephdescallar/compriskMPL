#' Summarise a phcsh_mpl object
#'
#' Additional information about the Cox proportional hazard model fit
#' represented by object is extracted and included in the returned object,
#' which is suitable for printing with the generic print function. The generic
#' coef function will extract the matrix of coefficients of interest with
#' standard errors, zz-statistics and p-values.
#'
#' @param object An object inheriting from class coxph_mpl, representing a
#' fitted Cox proportional hazard model.
#'
#' @param sand use sandwich estimates of standard error
#'
#' @export summary.phcsh_mpl
summary.phcsh_mpl <- function(object, sand = FALSE, grad=FALSE){
  if(sand == FALSE){
    col.names = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)",
                  "Gradient")
  }
  else if(sand == TRUE){
    col.names = c("Estimate", "Std. Error (sandwich)", "z-value", "Pr(>|z|)",
                  "Gradient")
  }
  var.names = dimnames(object$data$X)[[2]]
  risks = list()
  if(object$cf==1){
    matzG = cbind(object$"gamma",object$seG, object$"gamma" /
                    object$seG,2*(1-stats::pnorm(abs(object$"gamma" /
                                                       object$seG))), object$gamma.gradient)
    colnames(matzG) = col.names
  }
  for(r in 1:object$n.risk){
    if(sand == FALSE){
      matxB = cbind(object$"beta"[[r]],object$seB[[r]], object$"beta"[[r]] /
                      object$seB[[r]],2*(1-stats::pnorm(abs(object$"beta"[[r]] /
                                                              object$seB[[r]]))), object$beta.gradient[[r]])
      matxT = cbind(object$"theta"[[r]], object$seT[[r]], object$"theta"[[r]] /
                      object$seT[[r]],2*(1-stats::pnorm(abs(object$"theta"[[r]] /
                                                              object$seT[[r]]))), object$theta.gradient[[r]])
    }
    else if(sand == TRUE){
      matxB = cbind(object$"beta"[[r]],object$seB.sand[[r]], object$"beta"[[r]]
                    / object$seB.sand[[r]],2*(1-stats::pnorm(abs(object$"beta"[[r]] /
                                                                   object$seB.sand[[r]]))), object$beta.gradient[[r]])
      matxT = cbind(object$"theta"[[r]], object$seT.sand[[r]],
                    object$"theta"[[r]] / object$seT.sand[[r]],2 *
                      (1-stats::pnorm(abs(object$"theta"[[r]] / object$seT.sand[[r]]))),
                    object$theta.gradient[[r]])
    }
    dimnames(matxB) = list(paste(" ", var.names, sep = ""),col.names)
    colnames(matxT) = col.names
    if(grad==FALSE){
      risks[[r]] = list(Beta = matxB[, -ncol(matxB)],
                        Theta = matxT[, -ncol(matxT)])
    }
    else{
      risks[[r]] = list(Beta = matxB, Theta = matxT)
    }
  }
  if(grad==FALSE & object$cf==1){
    matzG=matzG[,-ncol(matzG)]
  }
  if(object$cf==1){
    out <- list(gamma=matzG, risks=risks)
  }
  else {out <- list(risks=risks)}
  out
}
