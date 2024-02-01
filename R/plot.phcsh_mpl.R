#' Plot a phcsh_mpl object
#'
#' Plot the bases used to estimate the baseline hazard parameters, as well as
#' the estimate and confidence interval of the baseline hazard, cumulative
#' baseline hazard and baseline survival functions (plots are selectable by
#' which).
#'
#' @param object an object inheriting from class coxph_mpl, representing a
#' fitted cause specific proportional hazards model.
#'
#' @param risk a number which specifies which risk  to plot from phcsh_mpl
#' object
#'
#' @param plots determines which type of plot to produce. bh refers to baseline
#' hazard plot. surv refers to survival plot, and cif refers to cumulative
#' incidence function plot. If none are specified, then all three plots are
#' displayed.
#'
#' @param sand use sandwich estimates of standard error
#'
#' @param n.points number of points to use for plot. Deafult is n = 1000.
#'
#' @export plot.phcsh_mpl
plot.phcsh_mpl <- function(object,risk=1, plots = c("bh", "surv", "cif"),
                           sand = FALSE, n.points = 1000, pred=NULL){
  r = risk

  if("bh" %in% plots) bh = 1
  else bh = 0
  if("surv" %in% plots) surv = 1
  else surv = 0
  if("cif" %in% plots) cif = 1
  else cif = 0
  psif <- function(x, bknots, iknots){
    splines2::mSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
                      intercept = object$basis.intercept)
  }
  PSIf <- function(x, bknots, iknots){
    splines2::iSpline(x, knots = iknots, Boundary = bknots, degree = object$dgr,
                      intercept = object$basis.intercept)
  }
  if(is.null(pred)){
    t.points = seq(object$b.knots[1], object$b.knots[2], length.out = n.points)
    plot.psi = psif(t.points, object$b.knots, object$i.knots[[r]])
    plot.h0r = plot.psi %*% object$"theta"[[r]]
    plot.h0r[plot.h0r < 1e-6] <- 1e-6
    if(sand == FALSE){
      VarCovMat.theta = object$VarCovMat[object$theta.index[[r]],
                                         object$theta.index[[r]]]
    }
    else if(sand == TRUE){
      VarCovMat.theta = object$sand[object$theta.index[[r]],
                                    object$theta.index[[r]]]
    }
    plot.h0r.var = diag(plot.psi %*% VarCovMat.theta %*% t(plot.psi))
    if(object$pos.def==1){
      plot.se = sqrt(plot.h0r.var) / (plot.h0r)
      plot.h0r.lower = as.vector(plot.h0r) * exp(-1.96*plot.se)
      plot.h0r.upper = as.vector(plot.h0r) * exp(1.96*plot.se)
    }
    #Baseline survival
    plot.PSI = PSIf(t.points, object$b.knots, object$i.knots[[r]])
    plot.S0r = exp(-(plot.PSI %*% object$"theta"[[r]]))
    plot.S0r.logOR = log((1-plot.S0r) / plot.S0r)
    plot.S0r.prime2 = plot.S0r^2
    plot.chaz.var = diag(plot.PSI %*% VarCovMat.theta %*% t(plot.PSI))
    plot.S0r.var = plot.S0r.prime2 * plot.chaz.var
    plot.S0r.logOR.var = ((1/((1-plot.S0r)*plot.S0r))^2)*plot.S0r.var
    if(object$pos.def==1){
      plot.S0r.logOR.lower = plot.S0r.logOR + 1.96*sqrt(plot.S0r.logOR.var)
      plot.S0r.logOR.upper = plot.S0r.logOR - 1.96*sqrt(plot.S0r.logOR.var)
      plot.S0r.lower = exp(-plot.S0r.logOR.lower) / (1 +
                                                       exp(-plot.S0r.logOR.lower))
      plot.S0r.upper = exp(-plot.S0r.logOR.upper) / (1 +
                                                       exp(-plot.S0r.logOR.upper))
    }
    #Cumulative incidence function
    bma = (t.points - object$b.knots[1]) / 2
    bpa = (t.points + object$b.knots[1]) / 2
    t.gq.change = bma*matrix(object$nodes, nrow=length(bma), ncol =
                               object$gq.points, byrow=TRUE) + bpa
    plot.F0r.psi.gq = plot.F0r.PSI.gq = plot.F0r.h0qt.gq = list()
    plot.F0r.H0qt.gq = plot.F0r.S0qt.gq = plot.F0r.Integrand.gq = list()
    plot.F0r.dhdT.gq = plot.F0r.dSdT.gq = plot.F0r.dFdT.gqT = list()
    for(gq in 1:object$gq.points){
      plot.F0r.psi.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                    object$i.knots[[r]])
      plot.F0r.PSI.gq[[gq]] <- PSIf(t.gq.change[,gq], object$b.knots,
                                    object$i.knots[[r]])
      plot.F0r.h0qt.gq[[gq]] <- psif(t.gq.change[,gq], object$b.knots,
                                     object$i.knots[[r]]) %*% object$"theta"[[r]]
      plot.F0r.H0qt.gq.r = plot.F0r.S0qt.gq.r = list()
      for(q in 1:object$n.risk){
        plot.F0r.H0qt.gq.r[[q]] = PSIf(t.gq.change[,gq], object$b.knots,
                                       object$i.knots[[q]]) %*% object$"theta"[[q]]
        plot.F0r.S0qt.gq.r[[q]] = exp(-plot.F0r.H0qt.gq.r[[q]])
      }
      plot.F0r.S0qt.gq[[gq]] = Reduce("*", plot.F0r.S0qt.gq.r)
      plot.F0r.Integrand.gq[[gq]] <- object$weights[gq]*plot.F0r.h0qt.gq[[gq]] *
        plot.F0r.S0qt.gq[[gq]]
      plot.F0r.dhdT.gq[[gq]] <- plot.F0r.psi.gq[[gq]]
      plot.F0r.dSdT.gq[[gq]] <- -plot.F0r.PSI.gq[[gq]] *
        as.vector(plot.F0r.S0qt.gq[[gq]])
      plot.F0r.dFdT.gqT[[gq]] <- (as.matrix(plot.F0r.dhdT.gq[[gq]] *
                                              as.vector(plot.F0r.S0qt.gq[[gq]]) +
                                              as.vector(plot.F0r.h0qt.gq[[gq]]) *
                                              as.matrix(plot.F0r.dSdT.gq[[gq]]))) *
        object$weights[gq]
    }
    plot.F0r <- bma*Reduce("+",plot.F0r.Integrand.gq)
    plot.F0r.dFdT <- bma*Reduce("+",plot.F0r.dFdT.gqT)
    plot.F0r.r.var <- diag(plot.F0r.dFdT %*% as.matrix(VarCovMat.theta) %*%
                             t(plot.F0r.dFdT))
    plot.F0r.logOR <- log(plot.F0r / (1-plot.F0r + 1e-12) + 1e-12)
    plot.F0r.logOR.var = ((1/((1-plot.F0r)*plot.F0r))^2)*plot.F0r.r.var
    plot.F0r.log <- log(plot.F0r + 1e-12)
    plot.F0r.log.var <- (1/plot.F0r^2)*plot.F0r.r.var
    if(object$pos.def==1){
      plot.F0r.logOR.lower = plot.F0r.logOR - 1.96*sqrt(plot.F0r.logOR.var)
      plot.F0r.logOR.upper = plot.F0r.logOR + 1.96*sqrt(plot.F0r.logOR.var)
      plot.F0r.lower <- exp(plot.F0r.logOR.lower) / (1 + exp(plot.F0r.logOR.lower))
      plot.F0r.upper <- exp(plot.F0r.logOR.upper) / (1 + exp(plot.F0r.logOR.upper))
    }
  }

  #Cumulative incidence function
  # bma = (t.points - object$b.knots[1]) / 2
  # bpa = (t.points + object$b.knots[1]) / 2
  # t.bma.gq = sweep(matrix(bma, nrow = length(t.points),
  # ncol = object$gq.points), MARGIN = 2, object$nodes, `*`) + bpa
  # t.bma.gq.l = split(t.bma.gq, rep(1:object$gq.points, each = length(t.points)))
  # t.bma.gq.w.psi = lapply(t.bma.gq.l, function(a) psif(a,
  #                  object$b.knots, object$i.knots[[risk]]))
  # t.bma.gq.w = sweep(matrix(bma, nrow = length(t.points),
  #             ncol = object$gq.points), MARGIN = 2, object$weights, `*`)
  # t.h0.q <- sapply(t.bma.gq.w.psi, function(a) a %*% object$theta[[risk]])
  # t.gq.PSI.r <- t.gq.H0.r <- t.gq.S0r.r <- list()
  # for(r in 1:object$n.risk){
  #   t.gq.PSI.r[[r]] = lapply(t.bma.gq.l, function(a) PSIf(a,
  #                      object$b.knots, object$i.knots[[r]]))
  #   t.gq.H0.r[[r]] = sapply(t.gq.PSI.r[[r]], function(a) a %*%
  #                              object$theta[[r]])
  #   t.gq.S0r.r[[r]] = exp(-t.gq.H0.r[[r]])
  # }
  # t.S.gq.q = Reduce("*", t.gq.S0r.r)
  # t.F.q = rowSums(t.h0.q * t.S.gq.q *
  #                         t.bma.gq.w)
  # t.S.q = rowSums(t.S.gq.q * t.bma.gq.w)
  # t.h0.q.sum = rowSums(t.h0.q * t.bma.gq.w)
  # t.H0.q <- PSIf(t.points, object$b.knots, object$i.knots[[risk]]) %*% object$theta[[risk]]

  else{
    t.points = pred$t
    plot.h0r = pred$pred.h0r
    plot.h0r.lower = pred$pred.h0r.lower
    plot.h0r.upper = pred$pred.h0r.upper
    plot.S0r = pred$pred.S0r
    plot.S0r.lower = pred$pred.S0r.lower
    plot.S0r.upper = pred$pred.S0r.upper
    plot.F0r = pred$pred.F0r
    plot.F0r.lower = pred$pred.F0r.lower
    plot.F0r.upper = pred$pred.F0r.upper
    plot.h0r[is.infinite(plot.h0r)] <- NA
    plot.h0r.lower[is.infinite(plot.h0r.lower)] <- NA
    plot.h0r.upper[is.infinite(plot.h0r.upper)] <- NA
    plot.S0r[is.infinite(plot.S0r)] <- NA
    plot.S0r.lower[is.infinite(plot.S0r.lower)] <- NA
    plot.S0r.upper[is.infinite(plot.S0r.upper)] <- NA
    plot.F0r[is.infinite(plot.F0r)] <- NA
    plot.F0r.lower[is.infinite(plot.F0r.lower)] <- NA
    plot.F0r.upper[is.infinite(plot.F0r.upper)] <- NA
  }
  plot.bh = function(h0,low=NULL,up=NULL){
    if(object$pos.def==1){
      plot(t.points, h0, xlab = "t", main = paste("Risk", risk),
           ylab = "baseline hazard", ylim = c(min(low, na.rm = TRUE),
                                              max(plot.h0r.upper, na.rm = TRUE)), type = "l")
      graphics::lines(t.points, low, lty = "dashed")
      graphics::lines(t.points, up, lty = "dashed")
    }
    else{
      plot(t.points, h0, xlab = "t", main = paste("Risk", risk), ylab =
             "baseline hazard", ylim = c(min(h0, na.rm = TRUE), max(h0, na.rm = TRUE)),
           type = "l")
    }
  }
  plot.surv = function(S0r,low=NULL,up=NULL){
    title = if(bh==0) paste("Risk", risk)
    else ""
    plot(t.points, S0r, xlab = "t", ylab = "Survival", ylim = c(0,1),
         main = title, type = "l")
    if(object$pos.def==1){
      graphics::lines(t.points, low, lty = "dashed")
      graphics::lines(t.points, up, lty = "dashed")
    }
  }
  plot.cif = function(F0r, low=NULL, up=NULL){
    title = if(surv==0 & bh == 0) paste("Risk", risk)
    else ""
    plot(t.points, F0r, xlab = "t", ylab = "CIF", ylim = c(0,1),
         main = title, type = "l")
    if(object$pos.def==1){
      graphics::lines(t.points, low, lty = "dashed")
      graphics::lines(t.points, up, lty = "dashed")
    }
  }
  graphics::par(mfrow = c(sum(bh, surv, cif),1))
  if(bh == 1) {
    if(object$pos.def==1)  plot.bh(plot.h0r, plot.h0r.lower, plot.h0r.upper)
    else plot.bh(plot.h0r)
  }
  if(surv == 1){
    if(object$pos.def==1) plot.surv(plot.S0r,plot.S0r.lower,plot.S0r.upper)
    else plot.surv(plot.S0r)
  }
  if(cif == 1){
    if(object$pos.def==1) plot.cif(plot.F0r, plot.F0r.lower, plot.F0r.upper)
    else plot.cif(plot.F0r)
  }
  plotValue <- list("t.points"=t.points, "plot.h0r"=plot.h0r,
                    "plot.h0r.lower"=plot.h0r.lower, "plot.h0r.upper"=plot.h0r.upper,
                    "plot.S0r"=plot.S0r, "plot.S0r.lower"=plot.S0r.lower,
                    "plot.S0r.upper"=plot.S0r.upper, "F0r"=plot.F0r,
                    "plot.F0r.lower"=plot.F0r.lower, "plot.F0r.upper"=plot.F0r.upper)
}
















