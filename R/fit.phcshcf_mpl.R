#' Function file for phcsh_mpl
#'
#' This is the function file for model fit for cause-specific hazards model
#' with a cure model using maximum penalized likelihood estimation.
#'
#' @export
fit.phcshcf_mpl <- function(formula, risk, z, data, control, ...){
  mc = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mc), 0)
  mc = mc[c(1, m)]
  if (m[1] == 0){
    stop("A formula argument is required")
  }
  else {
    "-"
  }
  mc[[1]] = as.name("model.frame")
  mc$formula = if (missing(data))
    stats::terms(formula)
  else stats::terms(formula, data = data)
  mf = eval(mc, parent.frame())
  mf_indx <- as.numeric(row.names(mf))
  y = stats::model.extract(mf, "response")
  if (nrow(mf) == 0)
    stop("No (non-missing) observations")
  mt = attr(mf, "terms")
  type = attr(y, "type")
  if (!inherits(y, "Surv")) {
    stop("Response must be a survival object")
  }
  if (attr(y, which = "type") == "right") {
    left = y[, 1]
    right = rep(NA, nrow(y))
    icase = which(y[, 2] == 1)
    right[icase] = y[icase, 1]
    y = survival::Surv(left, right, type = "interval2")
  }
  else if (type != "interval") {
    stop("\nPlease create the survival object using the option type='interval2'
         in the Surv function.\n")
  }
  n = nrow(mf)
  extraArgs <- list(...)
  if(length(extraArgs)){
    controlArgs <- names(formals(control.phcsh_mpl))
    m <- pmatch(names(extraArgs), controlArgs, nomatch = 0L)
    if (any(m == 0L))
      stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m ==
                                                                     0L]), domain = NA, call. = F)
  }
  if(missing(control))
    control <- control.phcsh_mpl(...)
  index = as.vector(row(mf)[,1])
  X = stats::model.matrix(mt, mf)
  X = X[, !apply(X, 2, function(x) all(x == x[1])), drop = FALSE]
  if (ncol(X) == 0) {
    X = matrix(0, n, 1)
    noX = TRUE
  }
  else {
    noX = FALSE
  }
  p = ncol(X)
  X.dataframe <- data.frame(X)
  xnames <- colnames(X.dataframe)
  mean_j = apply(X, 2, mean)
  XC = X - rep(mean_j, each = n)
  risk = risk[index]
  z = z[index,]
  Z = z[, !apply(z, 2, function(z) all(z == z[1])), drop = FALSE]
  Z = cbind("intercept" = rep(1,n),Z)
  # colnames(Z) <- z_names HARD CODED, FIX THIS!!
  #colnames(Z) <- c("z0","z1","z2") #HARD CODED, FIX THIS!!
  #factor(treatment.n) + factor(sex) + factor(Intracranial_Metastases_No) + factor(Extracranial_Dx) + factor(systemic) + factor(steroid)
  znames <- c(colnames(Z))
  o <- ncol(Z)
  mean_o = c(0,apply(Z[,-1], 2, mean))
  ZC = Z - rep(mean_o, each = n)
  data.export = data.frame(y,X,Z,risk)
  data.all = data.frame(y,X,risk,Z)
  data.all[,1][,2] = ifelse(data.all[,1][,3]==2,
                            data.all[,1][,1],data.all[,1][,2])
  data.all[,1][,1] = ifelse(data.all[,1][,3]==2, 0, data.all[,1][,1])
  l.risk <- sort(unique(risk))
  n.risk <- length(l.risk)
  tr = data.all[which(data.all[,1][,3]==0),1][,1]
  if(length(tr)==0) tr = NULL
  xr = as.matrix(data.all[which(data.all[,1][,3]==0), xnames])
  zr = as.matrix(data.all[which(data.all[,1][,3]==0), znames])
  n.r = length(tr)
  t.min = min(data.all[,1][,c(1,2)])
  t.max = max(data.all[,1][,c(1,2)])
  #b.knots <- c(t.min - 1e-12, t.max)
  b.knots <- c(t.min, t.max)
  knots.perc.limit = control$knots.perc.limit
  gq.points = control$gq.points
  dgr = control$dgr
  basis.intercept = control$basis.intercept
  gq <- statmod::gauss.quad(gq.points, "legendre")
  psif <- function(x, bknots, iknots){splines2::mSpline(x, knots = iknots,
                                                        Boundary = bknots, degree = dgr, intercept = basis.intercept)}
  PSIf <- function(x, bknots, iknots) {splines2::iSpline(x, knots = iknots,
                                                         Boundary = bknots, degree = dgr,intercept = basis.intercept)}
  te <- xe <- tl <- xl <- til <- tir <- xi <- n.e <- n.l <- n.i <- list()
  n.per.risk <- n.iknots <- i.knots <- perc.iknots <- t.all.r <- list()
  ti.rml <- ti.rpl <- ti.rml.gq <- ti.rml.gq.w <- ti.rml.gq.w.l <- list()
  ti.rml.gq.l <- ti.rml.gq.w.psi <- tmid <- list()
  n.basis <- R.mat <- oldbeta <- newbeta <- oldtheta <- newtheta <- list()
  tr.PSI <- te.psi <- oldlambda <- oldgamm <- newgamm <- ze <- zi <- list()
  for(q in 1:n.risk){
    te[[q]] = data.all[which(risk==q & data.all[,1][,3]==1),1][,1]
    xe[[q]] = as.matrix(data.all[which(risk==q & data.all[,1][,3]==1), xnames])
    ze[[q]] <- as.matrix(data.all[which(risk==q & data.all[,1][,3]==1), znames])
    n.e[[q]] = length(te[[q]])
    if(n.e[[q]] == 0) {
      te[[q]] = NA
      xe[[q]] = NA
    }
    til[[q]] = data.all[which(risk==q & (data.all[,1][,3]==3 |
                                           data.all[,1][,3]==2)),1][,1]
    tir[[q]] = data.all[which(risk==q & (data.all[,1][,3]==3 |
                                           data.all[,1][,3]==2)),1][,2]
    tmid[[q]] = rowMeans(cbind(til[[q]],tir[[q]]))
    xi[[q]] = as.matrix(data.all[which(risk==q & (data.all[,1][,3]==3 |
                                                    data.all[,1][,3]==2)), xnames])
    zi[[q]] = as.matrix(data.all[which(risk==q & (data.all[,1][,3]==3 |
                                                    data.all[,1][,3]==2)), znames])
    n.i[[q]] = length(til[[q]])
    if(n.i[[q]] == 0){
      til[[q]] = NA
      tir[[q]] = NA
      xi[[q]] = NA
    }
    ti.rml[[q]] = (tir[[q]] - til[[q]]) / 2
    ti.rpl[[q]] = (tir[[q]] + til[[q]]) / 2
    if(control$tmid==TRUE)
      t.all.r[[q]] = c(te[[q]], tmid[[q]])
    else
      t.all.r[[q]] = c(te[[q]], til[[q]], tir[[q]])
    n.per.risk[[q]] = n.e[[q]] + n.i[[q]]
    n.iknots[[q]] = max(floor(n.per.risk[[q]]^(1/3)) - 2,1)
    if(length(control$n.basis)!=0){
      n.basis[[q]] = control$n.basis[[q]]
      n.iknots[[q]] = n.basis[[q]] - 2
    }
    n.basis[[q]] = n.iknots[[q]] + dgr + basis.intercept
    if(n.iknots[[q]]==1)
      perc.iknots[[q]] = 0.5
    else if(n.iknots[[q]]==2)
      perc.iknots[[q]] = seq(knots.perc.limit[1], knots.perc.limit[2],
                             length.out = 4)[2:3]
    else
      perc.iknots[[q]] = seq(knots.perc.limit[1], knots.perc.limit[2],
                             length.out = n.iknots[[q]])
    i.knots[[q]]  = stats::quantile(unique(t.all.r[[q]]), prob =
                                      perc.iknots[[q]], names = FALSE, na.rm = TRUE, type = 3)
    if(!is.null(control$iknots.pos)){
      i.knots[[q]] = control$iknots.pos[[q]]
    }
    ti.rml.gq[[q]] = sweep(matrix(ti.rml[[q]], nrow = n.i[[q]],
                                  ncol = gq.points), MARGIN = 2, gq$nodes, `*`) + ti.rpl[[q]]
    ti.rml.gq.w[[q]] = sweep(matrix(ti.rml[[q]], nrow = n.i[[q]],
                                    ncol = gq.points), MARGIN = 2, gq$weights, `*`)
    ti.rml.gq.l[[q]] = split(ti.rml.gq[[q]], rep(1:gq.points,
                                                 each = n.i[[q]]))
    ti.rml.gq.w.l[[q]] = split(ti.rml.gq.w[[q]], rep(1:gq.points,
                                                     each = n.i[[q]]))
    ti.rml.gq.w.psi[[q]] = lapply(ti.rml.gq.l[[q]], function(a) psif(a,
                                                                     b.knots, i.knots[[q]]))
    Rtemp = matrix(0, nrow = n.basis[[q]], ncol = n.basis[[q]])
    xknots = c(rep(t.min, dgr), i.knots[[q]], rep(t.max, dgr))
    for(ii in 1:n.basis[[q]]){
      for(jj in 1:n.basis[[q]]){
        if(jj - ii<dgr){
          kntset = xknots[xknots >= xknots[jj] & xknots <= xknots[ii+dgr]]
          kntsum = 0
          for(kk in 1:(length(kntset) - 1)){
            kntsum = kntsum + splines2::mSpline(kntset[kk],
                                                knots = i.knots[[q]], degree = dgr, intercept =
                                                  basis.intercept, Boundary.knots = b.knots, derivs =
                                                  2)[ii]*splines2::mSpline(kntset[kk], knots = i.knots[[q]],
                                                                           degree = dgr, intercept = basis.intercept,
                                                                           Boundary.knots = b.knots, derivs = 2)[jj]*(kntset[kk + 1] -
                                                                                                                        kntset[kk])
          }
          Rtemp[ii, jj] = kntsum
        }
      }
    }
    R.mat[[q]] <- Rtemp
    R.mat[[q]][lower.tri(R.mat[[q]], diag = FALSE)] =
      t(R.mat[[q]])[lower.tri(R.mat[[q]], diag = FALSE)]
    oldbeta[[q]] = rep(0, p)
    newbeta[[q]] = rep(0, p)
    oldtheta[[q]] = rep(1, n.basis[[q]])
    newtheta[[q]] = rep(1, n.basis[[q]])
    tr.PSI[[q]] = if(length(tr)!=0) PSIf(tr, b.knots, i.knots[[q]])
    else NA
    te.psi[[q]] = if(n.e[[q]]!=0) psif(te[[q]], b.knots, i.knots[[q]])
    else NA
    if(is.null(control$lambda)) oldlambda[[q]] = 1e-5
    else oldlambda[[q]] = control$lambda[q]
  }
  valid=1
  oldgamm <- rep(0, o)
  newgamm <- rep(0, o)
  n.e.all = sum(unlist(n.e))
  n.i.all = sum(unlist(n.i))
  te.PSI.qr <- ti.gq.PSI.qr <- ti.gq.psi.qr <- te.psi.qr <- list()
  for(q in 1:n.risk){
    te.PSI.r <- ti.gq.PSI.r <- ti.gq.psi.r <- te.psi.r <- list()
    for(r in 1:n.risk){
      te.PSI.r[[r]] = if(n.e[[q]]!=0) PSIf(te[[q]], b.knots, i.knots[[r]])
      else NA
      ti.gq.PSI.r[[r]] = lapply(ti.rml.gq.l[[q]], function(a) PSIf(a,
                                                                   b.knots, i.knots[[r]]))
      ti.gq.psi.r[[r]] = lapply(ti.rml.gq.l[[q]], function(a) psif(a,
                                                                   b.knots, i.knots[[r]]))
      te.psi.r[[r]] = if(n.e[[q]]!=0) psif(te[[q]], b.knots, i.knots[[r]])
      else NA
    }
    te.PSI.qr[[q]] = te.PSI.r
    ti.gq.PSI.qr[[q]] = ti.gq.PSI.r
    ti.gq.psi.qr[[q]] = ti.gq.psi.r
    te.psi.qr[[q]] = te.psi.r
  }
  max.outer = control$max.outer
  max.iter = control$max.iter
  updbase <- function(theta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
                      ti.gq.PSI.qr){
    tr.H0.q = te.h0.q = ti.h0.q = ti.h0.q.l = te.H0.qr = ti.gq.H0.qr = list()
    ti.gq.S0r.qr = list()
    for(q in 1:n.risk){
      tr.H0.q[[q]] = tr.PSI[[q]] %*% theta[[q]]
      te.h0.q[[q]] = if(n.e[[q]]!=0) {te.psi[[q]] %*% theta[[q]]}
      else{NA}
      ti.h0.q[[q]] = sapply(ti.rml.gq.w.psi[[q]], function(a) a %*% theta[[q]])
      ti.h0.q.l[[q]] = split(ti.h0.q[[q]], rep(1:gq.points, each = n.i[[q]]))
      te.H0.r <- ti.gq.H0.r <- ti.gq.S0r.r <- list()
      for(r in 1:n.risk){
        te.H0.r[[r]] = if(n.e[[q]]!=0){te.PSI.qr[[q]][[r]] %*% theta[[r]]}
        else{NA}
        ti.gq.H0.r[[r]] = sapply(ti.gq.PSI.qr[[q]][[r]], function(a) a %*%
                                   theta[[r]])
        ti.gq.S0r.r[[r]] = if(n.i[[q]]!=0) exp(-ti.gq.H0.r[[r]])
      }
      te.H0.qr[[q]] = te.H0.r
      ti.gq.H0.qr[[q]] = ti.gq.H0.r
      ti.gq.S0r.qr[[q]] = ti.gq.S0r.r
    }
    out = list(tr.H0.q=tr.H0.q, te.h0.q=te.h0.q, ti.h0.q=ti.h0.q,
               ti.h0.q.l=ti.h0.q.l, te.H0.qr = te.H0.qr, ti.gq.H0.qr =
                 ti.gq.H0.qr, ti.gq.S0r.qr = ti.gq.S0r.qr)
    out
  }
  updllparms <- function(beta, xr, xe, tr.H0.q, xi, te.h0.q, te.H0.qr,
                         ti.h0.q, ti.gq.S0r.qr, ti.rml.gq.w, gamm, ze, zi, zr){
    xr.exb.q  = xe.exb.qr = tr.H.q = xi.exb.qr = te.h.q = te.H.qr = list()
    ti.h.q = ti.Sr.gq.qr = ti.S.gq.q = ti.F.q = ze.pi = zi.pi = list()
    zr.ezg <- pi.zr <- tr.H.r <- tr.Sr.r <- tr.Sr.pop.r <- te.S.r <- list()
    ze.ezg <- zi.ezg <- list()
    for(q in 1:n.risk){
      xr.exb.q[[q]] = exp(data.matrix(xr) %*% beta[[q]])
      tr.H.q[[q]] = tr.H0.q[[q]] * as.vector(xr.exb.q[[q]])
      xe.exb.r = xi.exb.r = te.H.r = ti.Sr.gq.r = list()
      for(r in 1:n.risk){
        xe.exb.r[[r]] = exp(data.matrix(xe[[q]]) %*% beta[[r]])
        xi.exb.r[[r]] = exp(data.matrix(xi[[q]]) %*% beta[[r]])
        te.H.r[[r]] = te.H0.qr[[q]][[r]] * as.vector(xe.exb.r[[r]])
        te.S.r[[r]] = exp(-te.H.r[[r]])
        ti.Sr.gq.r[[r]] = if(n.i[[q]]!=0)
          ti.gq.S0r.qr[[q]][[r]] ^ as.vector(xi.exb.r[[r]])
      }
      xe.exb.qr[[q]] = xe.exb.r
      xi.exb.qr[[q]] = xi.exb.r
      te.h.q[[q]] = te.h0.q[[q]] * as.vector(xe.exb.qr[[q]][[q]])
      te.H.qr[[q]] = te.H.r
      ti.h.q[[q]] = if(n.i[[q]]!=0)
        ti.h0.q[[q]] * as.vector(xi.exb.qr[[q]][[q]])
      ti.Sr.gq.qr[[q]] = ti.Sr.gq.r
      ti.S.gq.q[[q]] = Reduce("*", ti.Sr.gq.r)
      ti.F.q[[q]] = if(n.i[[q]]!=0) {rowSums(ti.h.q[[q]] * ti.S.gq.q[[q]] *
                                               ti.rml.gq.w[[q]])
      }
      else {NA}
      ze.ezg[[q]] <- exp(data.matrix(ze[[q]]) %*% gamm)
      zi.ezg[[q]] <- exp(data.matrix(zi[[q]]) %*% gamm)
      ze.pi[[q]] <- ze.ezg[[q]] / (1 + ze.ezg[[q]])
      zi.pi[[q]] <- zi.ezg[[q]] / (1 + zi.ezg[[q]])
      #zr.ezg[[q]] <- exp(data.matrix(zr) %*% gamm[[r]])
      #pi.zr[[q]] <- zr.ezg[[q]] / (1 + zr.ezg[[q]])
      tr.Sr.r[[q]] <- exp(-tr.H.q[[q]])
      #tr.Sr.pop.r[[q]] <- pi.zr[[q]]*tr.Sr.r[[q]] + (1 - pi.zr[[q]])
    }
    zr.ezg <- exp(data.matrix(zr) %*% gamm)
    zr.pi <- zr.ezg / (1+zr.ezg)
    tr.S <- Reduce("*", tr.Sr.r)
    tr.S.pop <- zr.pi*tr.S + (1 - zr.pi)
    out = list(xr.exb.q = xr.exb.q, xe.exb.qr = xe.exb.qr, tr.H.q = tr.H.q,
               xi.exb.qr = xi.exb.qr, te.h.q = te.h.q, te.H.qr = te.H.qr,
               ti.h.q = ti.h.q, ti.Sr.gq.qr = ti.Sr.gq.qr, ti.S.gq.q = ti.S.gq.q,
               ti.F.q = ti.F.q, zr.ezg=zr.ezg, tr.Sr.r=tr.Sr.r, tr.S=tr.S,
               zr.pi=zr.pi, tr.S.pop=tr.S.pop, ze.pi=ze.pi, zi.pi=zi.pi,
               ze.ezg=ze.ezg, zi.ezg=zi.ezg)
    out
  }
  gammascore.z <- data.matrix(Reduce("rbind",
                                     list(if(length(tr)!=0) zr, if(sum(unlist(n.e))!=0) Reduce("rbind", ze)),
                                     if(sum(unlist(n.i))!=0) Reduce("rbind", zi)))
  # gammascore.z <- data.matrix(if(length(tr)!=0) zr)
  gammascore.1 <- rep(1, nrow(gammascore.z))
  ze.all <- Reduce("rbind", ze)
  zi.all <- Reduce("rbind", zi)
  ze.all.t <- t(ze.all)
  zi.all.t <- t(zi.all)
  updscorebeta <- function(ti.h.q, ti.S.gq.q, ti.gq.H0.qr, xi.exb.qr,
                           ti.rml.gq.w){
    ti.A.qr = list()
    for(q in 1:n.risk){
      ti.A.r = list()
      for(r in 1:n.risk){
        ti.A.r[[r]] = if(n.i[[q]]!=0) {rowSums(ti.h.q[[q]] * ti.S.gq.q[[q]] *
                                                 ((q==r) - ti.gq.H0.qr[[q]][[r]] *
                                                    as.vector(xi.exb.qr[[q]][[r]])) * ti.rml.gq.w[[q]])}
        else {NA}
      }
      ti.A.qr[[q]] = ti.A.r
    }
    out = list(ti.A.qr = ti.A.qr)
    out
  }
  betascore.X <- betascore.1 <- betahess.X <- hess.AA.x <- list()
  theta.num.1 <- theta.den.1 <- list()
  for(q in 1:n.risk){
    betascore.X[[q]] = Reduce("rbind", list(if(length(tr)!=0) xr,
                                            if(n.e[[q]]!=0) xe[[q]], Reduce("rbind",
                                                                            mapply(function(a,b) if(a!=0) b, n.e, xe,
                                                                                   SIMPLIFY=FALSE)), Reduce("rbind", mapply(function(a,b)
                                                                                     if (a!=0) b, n.i, xi, SIMPLIFY = FALSE))))
    betascore.1[[q]] = rep(1, nrow(betascore.X[[q]]))
    betahess.X[[q]] = data.matrix(Reduce("rbind", list(if(length(tr)!=0) xr,
                                                       Reduce("rbind", mapply(function(a,b) if(a!=0) b, n.e, xe,
                                                                              SIMPLIFY=FALSE)), Reduce("rbind", mapply(function(a,b)
                                                                                if(a!=0) b, n.i, xi, SIMPLIFY = FALSE)))))
    hess.AA.x[[q]] = data.matrix(Reduce("rbind", list(if(length(tr)!=0) xr,
                                                      Reduce("rbind",mapply(function(a,b) if(a!=0) b, n.e, xe,
                                                                            SIMPLIFY = FALSE)), Reduce("rbind", mapply(function(a,b) if(a!=0) b, n.i,
                                                                                                                       xi, SIMPLIFY = FALSE)))))
    theta.num.1[[q]] = matrix(1, ncol = 1, nrow = n.e[[q]] + sum(unlist(n.i)))
    theta.den.1[[q]] = matrix(1, ncol = 1, nrow = n.r + sum(unlist(n.e)) +
                                sum(unlist(n.i)))
  }
  updscoretheta = function(theta, ti.S.gq.q, ti.gq.psi.qr, ti.rml.gq.w.psi,
                           ti.h0.q.l, xi.exb.qr){
    ti.S.gq.q.l = list()
    for(q in 1:n.risk){
      ti.S.gq.q.l[[q]] = if(n.i[[q]]!=0) {split(ti.S.gq.q[[q]],
                                                rep(1:ncol(ti.S.gq.q[[q]]), each =
                                                      nrow(ti.S.gq.q[[q]])))}
      else NA
    }
    ti.B1.qr = ti.B2.qr = list()
    for(q in 1:n.risk){
      ti.B1.r = ti.B2.r = list()
      for(r in 1:n.risk){
        if(q==r){
          ti.B1.r[[r]] = if(n.i[[q]]!=0) {Reduce("+", mapply(function(a,b,c) a *
                                                               as.vector(b) * as.vector(c), ti.gq.psi.qr[[q]][[r]],
                                                             ti.S.gq.q.l[[q]], ti.rml.gq.w.l[[q]],
                                                             SIMPLIFY = FALSE))}
          else NA
        }
        else{
          ti.B1.r[[r]] = if(n.i[[q]]!=0) {matrix(0, nrow = n.i[[q]],
                                                 ncol = n.basis[[r]])}
          else NA
        }
        ti.B2.r[[r]] = if(n.i[[q]]!=0) {Reduce("+", mapply(function(a,b,c,d) a *
                                                             b * c * d * as.vector(xi.exb.qr[[q]][[q]]),
                                                           ti.h0.q.l[[q]], ti.gq.PSI.qr[[q]][[r]], ti.S.gq.q.l[[q]],
                                                           ti.rml.gq.w.l[[q]], SIMPLIFY = FALSE))}
        else NA
      }
      ti.B1.qr[[q]] = ti.B1.r
      ti.B2.qr[[q]] = ti.B2.r
    }
    out = list(ti.S.gq.q.l = ti.S.gq.q.l, ti.B1.qr = ti.B1.qr,
               ti.B2.qr = ti.B2.qr)
    out
  }
  hess.AA.x.diag = as.matrix(Matrix::bdiag(hess.AA.x))
  updhessian <- function(ze.all, ze.pi.all, zr, tr.S, zr.pi.all, zi.all.t,
                         zi.all, tr.H.q, xr, tr.PSI, xr.exb, ti.S.gq.q, ti.h.q, gq.points, n.i, ti.gq.H0.qr,
                         xi.exb.qr, ti.rml.gq.w, ti.gq.psi.qr, ti.gq.PSI.qr,
                         te.H.qr, ti.F.q, ti.A.qr, te.PSI.qr, ti.B1.qr,
                         ti.B2.qr, te.psi.qr, theta, xr.exb.q, xe.exb.qr, lambda,
                         R.mat){
    CC <- -1*ze.all.t %*% diag(c(ze.pi.all*(1-ze.pi.all))) %*% ze.all +
      t(zr) %*% diag(c(((tr.S-1)*((zr.pi.all*(1-zr.pi.all)^3)/((1-zr.pi.all +
                                                                  zr.pi.all*tr.S)^2))))) %*% zr - zi.all.t %*% diag(c((zi.pi.all*(1-zi.pi.all)))) %*% zi.all
    CA <- CB <- list()
    for(r in 1:n.risk){
      # CA[[r]] <- -1*t(zr) %*% diag(c((tr.S*tr.H.q[[r]]*zr.pi.all*(1-zr.pi.all)) /
      #            (1-zr.pi.all+zr.pi.all*tr.S)^2)) %*% xr
      CA[[r]] <- -1*t(zr) %*% diag(c(((1-zr.pi.all+zr.pi.all*tr.S)*zr.pi.all*(1-zr.pi.all)*tr.S*tr.H.q[[r]]
                                      - ((tr.S-1)*zr.pi.all*(1-zr.pi.all)*zr.pi.all*tr.S*tr.H.q[[r]])) / (1-zr.pi.all+zr.pi.all*tr.S)^2)) %*% xr
      CB[[r]] <- -1*t(zr) %*% diag(c((tr.S*zr.pi.all*(1-zr.pi.all)*xr.exb.q[[r]]) /
                                       (1-zr.pi.all+zr.pi.all*tr.S)^2)) %*% tr.PSI[[r]]
    }


    ti.h.q.l = ti.gq.H.qr =  list()
    for(q in 1:n.risk){
      ti.h.q.l[[q]] = if(n.i[[q]]!=0) {split(ti.h.q[[q]], rep((1:gq.points),
                                                              each = n.i[[q]]))}
      else NA
      ti.gq.H.r = list()
      for(r in 1:n.risk){ if(n.i[[q]]!=0){ti.gq.H.r[[r]] =
        ti.gq.H0.qr[[q]][[r]] * as.vector(xi.exb.qr[[q]][[r]])}
        else NA
      }
      ti.gq.H.qr[[q]] = ti.gq.H.r
    }
    ti.AA.qrjtk = ti.AB.qrjtu = ti.BB1.qrutz = ti.BB2.qrutz = list()
    for(q in 1:n.risk){
      ti.AA.rjtk = ti.AB.rjtu = ti.BB1.rutz = ti.BB2.rutz = list()
      for(r in 1:n.risk){
        ti.AA.jtk = ti.AB.jtu = ti.BB1.utz = ti.BB2.utz = list()
        for(t in 1:n.risk){
          if(q==r & r==t & q == t){
            ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
                                                         ti.S.gq.q[[q]] * (1 - 3*ti.gq.H.qr[[q]][[r]] +
                                                                             ti.gq.H.qr[[q]][[r]]^2)) * ti.rml.gq.w[[q]])}
            else {NA}
          }
          else if(q==r & r!=t & q!=t){
            ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
                                                         ti.S.gq.q[[q]] * ti.gq.H.qr[[q]][[t]] *(ti.gq.H.qr[[q]][[r]] - 1))
                                                      * ti.rml.gq.w[[q]])}
            else {NA}
          }
          else if((q!=r & q==t & r!=t) | (q!=r & q!=t & r==t)){
            ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
                                                         ti.S.gq.q[[q]] * ti.gq.H.qr[[q]][[r]] *(ti.gq.H.qr[[q]][[t]] - 1))
                                                      * ti.rml.gq.w[[q]])}
            else {NA}
          }
          else if(q!=r & q!=t & r!=t){
            ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
                                                         ti.S.gq.q[[q]] * ti.gq.H.qr[[q]][[r]] * ti.gq.H.qr[[q]][[t]])
                                                      * ti.rml.gq.w[[q]])}
            else {NA}
          }
          ti.AB.ju = list()
          for(w in 1:gq.points){
            if(q==r & q==t & r==t){
              ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
                (ti.S.gq.q[[q]][,w] * ti.gq.psi.qr[[q]][[t]][[w]] +
                   ti.h.q[[q]][,w]* ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
                   ti.gq.PSI.qr[[q]][[t]][[w]] - 2 * ti.h.q[[q]][,w] *
                   ti.S.gq.q[[q]][,w] * ti.gq.PSI.qr[[q]][[t]][[w]] -
                   ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
                   ti.gq.psi.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
              }
              else {NA}
            }
            else if(q==r & q!=t & r!=t){
              ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
                (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                   ti.gq.H.qr[[q]][[r]][,w] * ti.gq.PSI.qr[[q]][[t]][[w]] -
                   ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                   ti.gq.PSI.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
              }
              else {NA}
            }
            else if(q!=r & q==t & r!=t){
              ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
                (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                   ti.gq.H.qr[[q]][[r]][,w] * ti.gq.PSI.qr[[q]][[t]][[w]] -
                   ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
                   ti.gq.psi.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
              }
              else {NA}
            }
            else if(q!=r & q!=t & r==t){
              ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
                (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                   ti.gq.H.qr[[q]][[r]][,w] * ti.gq.PSI.qr[[q]][[t]][[w]] -
                   ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                   ti.gq.PSI.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
              }
              else {NA}
            }
            else if(q!=r & q!=t & r!=t){
              ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
                (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                   ti.gq.PSI.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
              }
              else {NA}
            }
          }
          ti.AB.jtu[[t]] = Reduce("+", ti.AB.ju)
          ti.BB1.uz = ti.BB2.uz = list()
          for(u in 1:n.basis[[r]]){
            ti.BB1.z = ti.BB2.z = list()
            for(z in 1:n.basis[[t]]){
              ti.BB1 = ti.BB2 = list()
              for(w in 1:gq.points){
                if(n.i[[q]]!=0){
                  if(q==r & q==t & r==t){
                    ti.BB1[[w]] = -ti.S.gq.q[[q]][,w] *
                      ti.gq.psi.qr[[q]][[r]][[w]][,u] *
                      ti.gq.PSI.qr[[q]][[t]][[w]][,z] *
                      ti.rml.gq.w[[q]][,w]
                    ti.BB2[[w]] = (ti.S.gq.q[[q]][,w] *
                                     ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                                     ti.gq.psi.qr[[q]][[t]][[w]][,z] - ti.h.q[[q]][,w]
                                   * ti.S.gq.q[[q]][,w] * ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                                     ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                  else if(q==r & q!=t & r!=t){
                    ti.BB1[[w]] = -ti.S.gq.q[[q]][,w] *
                      ti.gq.psi.qr[[q]][[r]][[w]][,u] *
                      ti.gq.PSI.qr[[q]][[t]][[w]][,z] *
                      ti.rml.gq.w[[q]][,w]
                    ti.BB2[[w]] = (- ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                                     ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                                     ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                  else if(q!=r & q==t & r!=t){
                    ti.BB1[[w]] = 0
                    ti.BB2[[w]] = (ti.S.gq.q[[q]][,w] *
                                     ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                                     ti.gq.psi.qr[[q]][[t]][[w]][,z] - ti.h.q[[q]][,w]
                                   * ti.S.gq.q[[q]][,w] * ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                                     ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                  else if((q!=r & q!=t & r==t) | (q!=r & q!=t & r!=t)){
                    ti.BB1[[w]] = 0
                    ti.BB2[[w]] = (-ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                                     ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                                     ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                }
                else{
                  ti.BB1[[w]] = NA
                  ti.BB2[[w]] = NA
                }
              }
              ti.BB1.z[[z]] = Reduce("+", ti.BB1)
              ti.BB2.z[[z]] = Reduce("+", ti.BB2)
            }
            ti.BB1.uz[[u]] = ti.BB1.z
            ti.BB2.uz[[u]] = ti.BB2.z
          }
          ti.BB1.utz[[t]] = ti.BB1.uz
          ti.BB2.utz[[t]] = ti.BB2.uz
        }
        ti.AA.rjtk[[r]] = ti.AA.jtk
        ti.AB.rjtu[[r]] = ti.AB.jtu
        ti.BB1.rutz[[r]] = ti.BB1.utz
        ti.BB2.rutz[[r]] = ti.BB2.utz
      }
      ti.AA.qrjtk[[q]] = ti.AA.rjtk
      ti.AB.qrjtu[[q]] = ti.AB.rjtu
      ti.BB1.qrutz[[q]] = ti.BB1.rutz
      ti.BB2.qrutz[[q]] = ti.BB2.rutz
    }
    te.h0.q = list()
    for(q in 1:n.risk){
      te.h0.q[[q]] = te.psi.qr[[q]][[q]] %*% theta[[q]]
    }
    hess.AA.rt.mat = hess.AB.rt.mat = hess.AB.x = list()
    for(r in 1:n.risk){
      hess.AA.t.mat = hess.AB.t.mat = list()
      for(t in 1:n.risk){
        hess.AA.t.mat[[t]] = Matrix::Diagonal(x = c(if(r==t){if(length(tr)!=0)
          -(tr.H.q[[t]]*tr.S*zr.pi.all +  zr.pi.all*(1-zr.pi.all)*tr.S*tr.H.q[[t]]^2/ (1-zr.pi.all+zr.pi.all*tr.S)^2)} else{rep(0, n.r)},
          if(r==t){unlist(mapply(function(a,b) if(a!=0)
            -b[[t]], n.e, te.H.qr, SIMPLIFY = FALSE))}
          else{rep(0, sum(unlist(n.e)))},
          unlist(mapply(function(a,b,c,d) if(d!=0)
          {(a*b[[r]][[t]] - c[[r]]*c[[t]]) / (a^2)}, ti.F.q,
          ti.AA.qrjtk, ti.A.qr, n.i, SIMPLIFY = FALSE))))
      }
      hess.AA.rt.mat[[r]] = Reduce("rbind", hess.AA.t.mat)
    }
    hess.AA.mat = Reduce("cbind", hess.AA.rt.mat)
    hess.AA = t(hess.AA.x.diag) %*% hess.AA.mat %*% hess.AA.x.diag
    AA.rjtk = list()
    AA.mat = matrix(0, ncol = p*n.risk, nrow = p*n.risk)
    AA.count1 = 1
    AA.count2 = 1
    for(r in 1:n.risk){
      AA.jtk = list()
      for(j in 1:p){
        AA.tk = list()
        for(t in 1:n.risk){
          AA.k = list()
          for(k in 1:p){
            if(r==t){
              xj.r <- xr[,j]
              xk.r <- xr[,k]
              AA.tr = as.vector(-(tr.H.q[[r]]*zr.pi.all*tr.S*((1 - zr.pi.all+zr.pi.all*tr.S)*(1 - tr.H.q[[t]]) + zr.pi.all*tr.S*tr.H.q[[t]])) / (1 - zr.pi.all+zr.pi.all*tr.S)^2)
              xj.e = xj.i = xk.e = xk.i = AA.te = AA.ti = vector()
              for(q in 1:n.risk){
                if(n.e[[q]] != 0) {
                  xj.e = c(xj.e, xe[[q]][, j])
                  xk.e = c(xk.e, xe[[q]][, k])
                  AA.te = c(AA.te, -te.H.qr[[q]][[t]])
                }
                if(n.i[[q]] != 0){
                  xj.i = c(xj.i, xi[[q]][, j])
                  xk.i = c(xk.i, xi[[q]][, k])
                  AA.ti = c(AA.ti, as.vector((ti.F.q[[q]] *
                                                ti.AA.qrjtk[[q]][[r]][[t]] - ti.A.qr[[q]][[r]] *
                                                ti.A.qr[[q]][[t]]) / ti.F.q[[q]]^2))
                }
              }
              xj = c(xj.r, xj.e, xj.i)
              xk = c(xk.r, xk.e, xk.i)
              AA.elem = c(AA.tr, AA.te, AA.ti)
              AA.temp = t(xj) %*% diag(AA.elem) %*% xk
            }
            else{
              xj.r <- xr[,j]
              xk.r <- xr[,k]
              AA.tr = as.vector(-(tr.H.q[[r]]*zr.pi.all*tr.S*(-(1 - zr.pi.all+zr.pi.all*tr.S)*tr.H.q[[t]] + zr.pi.all*tr.S*tr.H.q[[t]])) / (1 - zr.pi.all+zr.pi.all*tr.S)^2)
              xj.i = xk.i = AA.ti = vector()
              for(q in 1:n.risk){
                if(n.i[[q]] != 0){
                  xj.i = c(xj.i, xi[[q]][, j])
                  xk.i = c(xk.i, xi[[q]][, k])
                  AA.ti = c(AA.ti, as.vector((ti.F.q[[q]] *
                                                ti.AA.qrjtk[[q]][[r]][[t]] - ti.A.qr[[q]][[r]] *
                                                ti.A.qr[[q]][[t]]) / ti.F.q[[q]]^2))
                }
              }
              xj = c(xj.r, xj.i)
              xk = c(xk.r, xk.i)
              AA.elem = c(AA.tr, AA.ti)
              AA.temp = t(xj) %*% diag(AA.elem) %*% xk
            }
            AA.mat[AA.count1, AA.count2] = AA.temp
            AA.count2 = AA.count2 + 1
          }
        }
        AA.count1 = AA.count1 + 1
        AA.count2 = 1
      }
    }
    hess.AB.rjtu = list()
    for(r in 1:n.risk){
      hess.AB.jtu = list()
      for(j in 1:p){
        hess.AB.tu = list()
        for(t in 1:n.risk){
          hess.AB.u = list()
          for(u in 1:n.basis[[t]]){
            hess.AB.u[[u]] = sum(c(if(r==t){if(length(tr)!=0)
              -(zr.pi.all*tr.S*tr.PSI[[t]][,u]*xr[,j]*xr.exb.q[[t]]*(r==t) + zr.pi.all*(1-zr.pi.all)*tr.S*tr.H.q[[r]]*tr.PSI[[t]][,u]) / (1-zr.pi.all+zr.pi.all*tr.S)^2} else{rep(0,n.r)},
              if(r==t){unlist(mapply(function(a,b,c,d) if(d!=0)
                -a[[t]][,u]*b[,j]*c[[t]], te.PSI.qr, xe, xe.exb.qr, n.e,
                SIMPLIFY = FALSE))},
              unlist(mapply(function(a,b,c,d,e,f,g,h) if(h!=0) {f[,j]*g[[t]]*(
                (a*b[[r]][[t]][,u]) - (c[[r]]*(d[[t]][,u] - e[[t]][,u]))) / a^2},
                ti.F.q, ti.AB.qrjtu, ti.A.qr, ti.B1.qr, ti.B2.qr, xi, xi.exb.qr,
                n.i, SIMPLIFY = FALSE))))
          }
          hess.AB.tu[[t]] = Reduce("rbind",hess.AB.u)
        }
        hess.AB.jtu[[j]] = Reduce("rbind",hess.AB.tu)
      }
      hess.AB.rjtu[[r]] = Reduce("cbind",hess.AB.jtu)
    }
    AB.mat = matrix(0, nrow = p*n.risk, ncol = sum(unlist(n.basis)))
    AB.count1 = 1
    AB.count2 = 1
    for(r in 1:n.risk){
      for(j in 1:p){
        for(t in 1:n.risk){
          for(z in 1:n.basis[[t]]){
            if(r==t){
              AB.exb.r <- xr.exb.q[[t]]
              #if(n.r != 0) AB.tr <- -(zr.pi.all*tr.S*tr.PSI[[t]][,z]*((1 - zr.pi.all+zr.pi.all*tr.S)*(1-tr.H.q[[r]]) - tr.H.q[[r]]*tr.S*zr.pi.all)) / (1 - zr.pi.all+zr.pi.all*tr.S)^2
              #if(n.r != 0) AB.tr <- -(zr.pi.all*tr.S*(tr.PSI[[t]][,z]*(1 - zr.pi.all+zr.pi.all*tr.S)*(1-tr.H.q[[r]])) + tr.H.q[[r]]*tr.S*zr.pi.all*tr.PSI[[t]][,z]) / (1 - zr.pi.all+zr.pi.all*tr.S)^2
              if(n.r != 0) AB.tr <- -((1-zr.pi.all+zr.pi.all*tr.S)*zr.pi.all*tr.S*tr.PSI[[t]][,z]*(1-tr.H.q[[r]]) + (zr.pi.all^2)*(tr.S^2)*tr.H.q[[r]]*tr.PSI[[t]][,z]) / (1 - zr.pi.all+zr.pi.all*tr.S)^2
              else AB.tr = vector()
              xj.r <- vector()
              xj.r <- c(xj.r, xr[,j])
              xj.e <- xj.i <- AB.exb.e <- AB.exb.i <- AB.te <- AB.ti <- vector()
              for(q in 1:n.risk){
                if(n.e[[q]] != 0){
                  xj.e = c(xj.e, xe[[q]][,j])
                  AB.exb.e = c(AB.exb.e, xe.exb.qr[[q]][[t]])
                  AB.te = c(AB.te, -te.PSI.qr[[q]][[t]][,z])
                }
                if(n.i[[q]] != 0){
                  xj.i = c(xj.i, xi[[q]][,j])
                  AB.exb.i = c(AB.exb.i, xi.exb.qr[[q]][[t]])
                  AB.ti = c(AB.ti, ((ti.F.q[[q]]*ti.AB.qrjtu[[q]][[r]][[t]][,z])
                                    - (ti.A.qr[[q]][[r]]*(ti.B1.qr[[q]][[t]][,z] -
                                                            ti.B2.qr[[q]][[t]][,z]))) / ti.F.q[[q]]^2)
                }
              }
              xj = c(xj.r, xj.e, xj.i)
              AB.elem = c(AB.tr, AB.te, AB.ti)
              AB.exb = c(AB.exb.r, AB.exb.e, AB.exb.i)
              AB.temp = t(xj) %*% diag(AB.elem) %*% AB.exb
            }
            else{
              #if(n.r != 0) AB.tr <- -(zr.pi.all*tr.S*tr.PSI[[t]][,z]*((-zr.pi.all+zr.pi.all*tr.S)*(-tr.H.q[[r]]) - tr.H.q[[r]]*tr.S*zr.pi.all)) / (1 - zr.pi.all+zr.pi.all*tr.S)^2
              #if(n.r != 0) AB.tr <- -(zr.pi.all*tr.S*(tr.PSI[[t]][,z]*(1 - zr.pi.all+zr.pi.all*tr.S)*(1-tr.H.q[[r]])) + tr.H.q[[r]]*tr.S*zr.pi.all*tr.PSI[[t]][,z]) / (1 - zr.pi.all+zr.pi.all*tr.S)^2
              if(n.r != 0) AB.tr <- -((1-zr.pi.all+zr.pi.all*tr.S)*zr.pi.all*tr.S*tr.PSI[[t]][,z]*(-tr.H.q[[r]]) + (zr.pi.all^2)*(tr.S^2)*tr.H.q[[r]]*tr.PSI[[t]][,z]) / (1 - zr.pi.all+zr.pi.all*tr.S)^2
              else AB.tr = vector()
              xj.r <- vector()
              xj.r <- c(xj.r, xr[,j])
              xj.i = AB.exb.i = AB.ti  = vector()
              for(q in 1:n.risk){
                if(n.i[[q]] != 0){
                  xj.i = c(xj.i,xi[[q]][,j])
                  AB.exb.i = c(AB.exb.i, xi.exb.qr[[q]][[t]])
                  AB.ti = c(AB.ti,((ti.F.q[[q]]*ti.AB.qrjtu[[q]][[r]][[t]][,z]) -
                                     (ti.A.qr[[q]][[r]]*(ti.B1.qr[[q]][[t]][,z] -
                                                           ti.B2.qr[[q]][[t]][,z]))) / ti.F.q[[q]]^2)
                }
              }
              xj = c(xj.r, xj.i)
              AB.elem = c(AB.tr, AB.ti)
              exb = c(xr.exb.q[[t]], AB.exb.i)
              AB.temp = t(xj) %*% diag(AB.elem) %*% exb
            }
            AB.mat[AB.count1, AB.count2] = AB.temp
            AB.count2 = AB.count2 + 1
          }
        }
        AB.count1 = AB.count1 + 1
        AB.count2 = 1
      }
    }
    hess.AB = Reduce("cbind", hess.AB.rjtu)
    hess.BB.rtuz.mat = list()
    for(r in 1:n.risk){
      hess.BB.tuz.mat = list()
      for(u in 1:n.basis[[r]]){
        hess.BBtz.mat = list()
        for(t in 1:n.risk){
          hess.BBz.mat = list()
          for(z in 1:n.basis[[t]]){
            hess.BBz.mat[[z]] = sum(c(
              if(r==t & n.e[[r]]!=0)
              {-te.psi.qr[[r]][[r]][,u]*te.psi.qr[[t]][[t]][,z]/te.h0.q[[r]]^2}
              else{rep(0, n.e[[r]])}, unlist(mapply(function(a,b,c,d,e,f,g)
                if(g!=0) {f[[r]] * f[[t]] * (((a * (b[[r]][[t]][[u]][[z]] -
                                                      c[[r]][[t]][[u]][[z]])) - (d[[r]][,u] - e[[r]][,u]) * (d[[t]][,z] -
                                                                                                               e[[t]][,z]))) / a^2}, ti.F.q, ti.BB1.qrutz, ti.BB2.qrutz, ti.B1.qr,
                ti.B2.qr, xi.exb.qr, n.i, SIMPLIFY = FALSE))))
          }
          hess.BBtz.mat[[t]] = Reduce("rbind", hess.BBz.mat)
        }
        hess.BB.tuz.mat[[u]] = Reduce("rbind", hess.BBtz.mat)
      }
      hess.BB.rtuz.mat[[r]] = Reduce("cbind", hess.BB.tuz.mat)
    }
    hess.BB.mat = Reduce("cbind", hess.BB.rtuz.mat)
    BB.mat = matrix(0, nrow = sum(unlist(n.basis)), ncol = sum(unlist(n.basis)))
    BB.count1 = 1
    BB.count2 = 1
    for(r in 1:n.risk){
      for(u in 1:n.basis[[r]]){
        for(t in 1:n.risk){
          for(z in 1:n.basis[[t]]){
            BB.exb.r.e = BB.exb.t.e = BB.te = vector()
            if(r==t){
              if(n.r != 0){
                BB.exb.r.r <- xr.exb.q[[r]]
                BB.exp.t.r <- xr.exb.q[[t]]
                BB.tr <- -zr.pi.all*(1-zr.pi.all)*tr.PSI[[r]][,u]*tr.S*tr.PSI[[t]][,z] / (1-zr.pi.all+zr.pi.all*tr.S)^2
              }
              if(n.e[[r]] != 0){
                BB.exb.r.e = rep(1, n.e[[r]])
                BB.exb.t.e = rep(1, n.e[[r]])
                BB.te = -(te.psi.qr[[r]][[r]][,u]*te.psi.qr[[r]][[r]][,z]) /
                  te.h0.q[[r]]^2
              }
              BB.exb.r.i = BB.exb.t.i = BB.ti = vector()
              for(q in 1:n.risk){
                if(n.i[[q]] != 0){
                  BB.exb.r.i = c(BB.exb.r.i, xi.exb.qr[[q]][[r]])
                  BB.exb.t.i = c(BB.exb.t.i, xi.exb.qr[[q]][[t]])
                  BB.ti = c(BB.ti, (((ti.F.q[[q]] *
                                        (ti.BB1.qrutz[[q]][[r]][[t]][[u]][[z]] -
                                           ti.BB2.qrutz[[q]][[r]][[t]][[u]][[z]])) -
                                       (ti.B1.qr[[q]][[r]][,u] - ti.B2.qr[[q]][[r]][,u]) *
                                       (ti.B1.qr[[q]][[t]][,z] - ti.B2.qr[[q]][[t]][,z]))) /
                              ti.F.q[[q]]^2)
                }
              }
            }
            else{
              BB.exb.r.i = BB.exb.t.i = BB.ti = vector()
              if(n.r != 0){
                BB.exb.r.r <- xr.exb.q[[r]]
                BB.exp.t.r <- xr.exb.q[[t]]
                BB.tr <- -zr.pi.all*(1-zr.pi.all)*tr.PSI[[r]][,u]*tr.S*tr.PSI[[t]][,z] / (1-zr.pi.all+zr.pi.all*tr.S)^2
              }
              for(q in 1:n.risk){
                if(n.i[[q]] != 0){
                  BB.exb.r.i = c(BB.exb.r.i, xi.exb.qr[[q]][[r]])
                  BB.exb.t.i = c(BB.exb.t.i, xi.exb.qr[[q]][[t]])
                  BB.ti = c(BB.ti, (((ti.F.q[[q]] *
                                        (ti.BB1.qrutz[[q]][[r]][[t]][[u]][[z]] -
                                           ti.BB2.qrutz[[q]][[r]][[t]][[u]][[z]])) -
                                       (ti.B1.qr[[q]][[r]][,u] - ti.B2.qr[[q]][[r]][,u]) *
                                       (ti.B1.qr[[q]][[t]][,z] - ti.B2.qr[[q]][[t]][,z]))) /
                              ti.F.q[[q]]^2)
                }
              }
            }
            BB.exb.r = c(BB.exb.r.r, BB.exb.r.e, BB.exb.r.i)
            BB.exb.t = c(BB.exp.t.r, BB.exb.t.e, BB.exb.t.i)
            BB.elem = c(BB.tr, BB.te, BB.ti)
            BB.temp = t(BB.exb.r) %*% diag(BB.elem) %*% BB.exb.t
            BB.mat[BB.count1, BB.count2] = BB.temp
            BB.count2 = BB.count2 + 1
          }
        }
        BB.count1 = BB.count1 + 1
        BB.count2 = 1
      }
    }

    # AA <- list()
    # for(r in 1:n.risk){
    #   A <- list()
    #   for(t in 1:n.risk){
    #     A[[t]] <- -t(xr) %*% diag(c((tr.H.q[[t]]*zr.pi.all*tr.S*(r==t) + zr.pi.all*(1 -
    #               zr.pi.all)*tr.S*tr.H.q[[t]]*tr.H.q[[r]]) / (1 -
    #               zr.pi.all+zr.pi.all*tr.S)^2)) %*% xr +
    #               -t(xe[[r]]) %*% diag(c(te.H.qr[[r]][[t]]*(r==t))) %*% xe[[r]] +
    #               -t(xi[[r]]) %*% diag(c(ti.F.q[[r]]*ti.AA.qrjtk[[r]][[r]][[t]] -
    #               ti.A.qr[[r]][[r]]))
    #   }
    #   AA[[r]] <- A
    # }
    lambdaR = as.matrix(Matrix::bdiag(mapply(function(a,b)  as.vector(a)*b,
                                             oldlambda, R.mat, SIMPLIFY = FALSE)))
    hess.BB = hess.BB.mat - 2*lambdaR
    BB.mat.pen = BB.mat - 2*lambdaR
    CC.mat <- CC
    CA.mat <- Reduce("cbind", CA)
    CB.mat <- Reduce("cbind", CB)
    hess.old = cbind(rbind(hess.AA, hess.AB), rbind(t(hess.AB), hess.BB))
    hess = cbind(rbind(CC.mat, t(CA.mat), t(CB.mat)), rbind(CA.mat, AA.mat, t(AB.mat)), rbind(CB.mat, AB.mat, BB.mat.pen))
    hess.no.pen = cbind(rbind(CC.mat, t(CA.mat), t(CB.mat)), rbind(CA.mat, AA.mat, t(AB.mat)), rbind(CB.mat, AB.mat, BB.mat))
    rlist <- list(CC, CA, CB, "hess"=hess)
    return(rlist)
  }
  betascore.mat <- betascore <- betahess.mat <- betahess <- betainc <- list()
  num.mat <- dJ <- num.exb <- thetascore.num <- den.mat <- den.exb <- list()
  thetascore.den <- thetascore.half <- thetascore <- list()
  for(outer in 1:max.outer){
    for(iter in 1:max.iter){
      if(control$iter.disp==TRUE)
        cat("Outer iteration", outer, "of", max.outer, ": Inner iteration",
            iter, "of",max.iter, "\r")
      prev.oldgamm <- oldgamm
      prev.oldbeta <- oldbeta
      prev.oldtheta <- oldtheta
      base <- updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
                      ti.gq.PSI.qr)
      llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                           base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                           oldgamm, ze, zi, zr)
      R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, oldtheta, R.mat,
                     oldlambda, SIMPLIFY = FALSE)
      llik0 = sum(log(unlist(llparms$te.h.q)+1e-12),
                  -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                  -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                  log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
      # gammascore.mat <- diag(c(if(n.r!=0) ((llparms$tr.S - 1) *
      #                   llparms$zr.pi*(1 - llparms$zr.pi)) / (1-llparms$zr.pi +
      #                   llparms$zr.pi*llparms$tr.S), if(sum(unlist(n.e))!=0) 1 /
      #                   (1 + unlist(llparms$ze.ezg)), if(sum(unlist(n.i))!=0) 1 /
      #                   (1 + unlist(llparms$zi.ezg))))
      # gammascore <- t(gammascore.z) %*% gammascore.mat %*%
      #               data.matrix(gammascore.1)
      ze.pi.all <- Reduce("rbind", llparms$ze.pi)
      zi.pi.all <- Reduce("rbind", llparms$zi.pi)
      zr.pi.all <- llparms$zr.pi
      tr.S <- llparms$tr.S
      gammascore <- ze.all.t %*% (1-ze.pi.all) +
        t(zr) %*% ((zr.pi.all*(1-zr.pi.all)*(tr.S-1))/(1-zr.pi.all+zr.pi.all*tr.S)) +
        zi.all.t %*% (1-zi.pi.all)
      # gammahess.mat <- diag(c(if(n.r!=0) -((llparms$tr.S-1)*(llparms$zr.pi*(1 -
      #             llparms$zr.pi)^3)) /
      #               (1 - llparms$zr.pi + llparms$zr.pi*llparms$tr.S)^2,
      #                     if(sum(unlist(n.e))!=0) -unlist(llparms$ze.ezg)/
      #                       (1 + unlist(llparms$ze.ezg))^2,
      #                     if(sum(unlist(n.i))!=0) -unlist(llparms$zi.ezg)/
      #                       (1 + unlist(llparms$zi.ezg))^2))
      # gammahess.mat <- diag(c(if(n.r!=0) -((llparms$tr.S - 1)*((llparms$tr.S.pop^2)*(( (llparms$zr.pi*(1 - llparms$zr.pi)^2) - ((llparms$zr.pi^2)*(1 - llparms$zr.pi)) + (llparms$zr.pi^2*(1 - llparms$zr.pi)^2) )))) / llparms$tr.S.pop^3,
      #                  if(sum(unlist(n.e))!=0) -unlist(llparms$ze.ezg)/
      #                   (1 + unlist(llparms$ze.ezg))^2,
      #                   if(sum(unlist(n.i))!=0) -unlist(llparms$zi.ezg)/
      #                   (1 + unlist(llparms$zi.ezg))^2))
      gammahess <- -1*ze.all.t %*% diag(c(ze.pi.all*(1-ze.pi.all))) %*% ze.all +
        #  t(zr) %*% diag(c(((tr.S-1)*((zr.pi.all*(1-zr.pi.all)^3)/((1-zr.pi.all + zr.pi.all*tr.S)^2))))) %*% zr -
        t(zr) %*% diag(c(((tr.S-1)*((zr.pi.all*(1-zr.pi.all)^3 - zr.pi.all^3*(1-zr.pi.all)*tr.S)/((1-zr.pi.all + zr.pi.all*tr.S)^2))))) %*% zr -
        zi.all.t %*% diag(c((zi.pi.all*(1-zi.pi.all)))) %*% zi.all
      # gammahess <- solve(-(t(gammascore.z) %*% gammahess.mat %*%
      #                           gammascore.z))
      gammainc <- solve(-gammahess) %*% gammascore
      newgamm = oldgamm + as.vector(gammainc)
      llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                           base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                           newgamm, ze, zi, zr)
      llik1 = sum(log(unlist(llparms$te.h.q)+1e-12),
                  -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                  -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                  log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
      ome = 0.6
      while(llik1 <= llik0){
        newgamm = oldgamm + (ome * as.vector(gammainc))
        llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                             base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                             newgamm, ze, zi, zr)
        llik1 = sum(log(unlist(llparms$te.h.q)+1e-12),
                    -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                    -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                    log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
        if(ome >= 1e-2) ome = ome * 0.6
        else if (ome < 1e-2 & ome >= 1e-5) ome = ome * 5e-2
        else if (ome < 1e-5 & ome > 1e-20) ome = ome * 1e-5
        else break
      }
      llik0 <- llik1
      oldgamm <- newgamm
      R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, oldtheta, R.mat,
                     oldlambda, SIMPLIFY = FALSE)
      llik1 = sum(log(unlist(llparms$te.h.q)+1e-12),
                  -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                  -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                  log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
      for(r in 1:n.risk){
        base <- updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
                        ti.gq.PSI.qr)
        llparms <- updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                              base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                              oldgamm, ze, zi, zr)
        betaparms <- updscorebeta(llparms$ti.h.q, llparms$ti.S.gq.q,
                                  base$ti.gq.H0.qr, llparms$xi.exb.qr, ti.rml.gq.w)
        betascore.mat[[r]] = diag(c(if(length(tr)!=0) -llparms$zr.pi *
                                      llparms$tr.H.q[[r]] * llparms$tr.S /
                                      llparms$tr.S.pop, rep(1, n.e[[r]]),
                                    unlist(mapply(function(a,b)
                                      if(a!=0) -b[[r]], n.e, llparms$te.H.qr, SIMPLIFY =
                                        FALSE)), unlist(mapply(function(a,b,c)
                                          if(a!=0) b[[r]] / c, n.i, betaparms$ti.A.qr,
                                          llparms$ti.F.q, SIMPLIFY = FALSE))))
        betascore[[r]] = t(betascore.X[[r]]) %*% betascore.mat[[r]] %*%
          betascore.1[[r]]
        betahess.mat[[r]] = diag(c(if(length(tr)!=0) (llparms$zr.pi *
                                                        llparms$tr.H.q[[r]] * llparms$tr.S *
                                                        (llparms$tr.S.pop * (llparms$tr.H.q[[r]]
                                                                             - 1) - llparms$zr.pi[[r]] * llparms$tr.S *
                                                           llparms$tr.H.q[[r]])) / llparms$tr.S.pop^2,
                                   unlist(mapply(function(a,b)
                                     if(a!=0) b[[r]], n.e, llparms$te.H.qr, SIMPLIFY =
                                       FALSE)),  unlist(mapply(function(a,b,c) if(a!=0)
                                         (b[[r]] / c) ^ 2, n.i, betaparms$ti.A.qr,
                                         llparms$ti.F.q, SIMPLIFY = FALSE))))
        betahess[[r]] = solve(t(betahess.X[[r]]) %*% betahess.mat[[r]] %*%
                                betahess.X[[r]])
        betainc[[r]] = betahess[[r]] %*% betascore[[r]]
        newbeta[[r]] = oldbeta[[r]] + as.vector(betainc[[r]])
        llparms <- updllparms(newbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                              base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                              oldgamm, ze, zi, zr)
        llik2 = sum(log(unlist(llparms$te.h.q)+1e-12),
                    -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                    -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                    log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
        ome = 0.6
        while(llik2 <= llik1){
          newbeta[[r]] = oldbeta[[r]] + (ome * as.vector(betainc[[r]]))
          llparms = updllparms(newbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                               base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                               oldgamm, ze, zi, zr)
          llik2 =  sum(log(unlist(llparms$te.h.q)+1e-12),
                       -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                       -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                       log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
          if(ome >= 1e-2) ome = ome * 0.6
          else if (ome < 1e-2 & ome >= 1e-5) ome = ome * 5e-2
          else if (ome < 1e-5 & ome > 1e-20) ome = ome * 1e-5
          else break
        }
        llik1 = llik2
        oldbeta = newbeta
      }
      R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, oldtheta, R.mat,
                     oldlambda, SIMPLIFY = FALSE)
      llik2 =  sum(log(unlist(llparms$te.h.q)+1e-12),
                   -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                   -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                   log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)), na.rm = TRUE)
      for(r in 1:n.risk){
        base <- updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
                        ti.gq.PSI.qr)
        llparms <- updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                              base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                              oldgamm, ze, zi, zr)
        thetaparms <- updscoretheta(oldtheta, llparms$ti.S.gq.q, ti.gq.psi.qr,
                                    ti.rml.gq.w.l, base$ti.h0.q.l, llparms$xi.exb.qr)
        num.mat[[r]] = Reduce("rbind", list(if(n.e[[r]]!=0) te.psi[[r]] /
                                              as.vector(base$te.h0.q[[r]]), Reduce("rbind",
                                                                                   mapply(function(a,b,c) if(a!=0) b[[r]] / as.vector(c),
                                                                                          n.i, thetaparms$ti.B1.qr, llparms$ti.F.q,
                                                                                          SIMPLIFY = FALSE))))
        dJ[[r]] = 2*(oldtheta[[r]] %*% R.mat[[r]])
        num.exb[[r]] = Reduce("rbind", list(matrix(1, nrow = n.e[[r]],
                                                   ncol = 1), Reduce("rbind", mapply(function(a,b) if(a!=0)
                                                   {b[[r]]}, n.i, llparms$xi.exb.qr,
                                                   SIMPLIFY = FALSE))))
        thetascore.num[[r]] = (as.vector(t(theta.num.1[[r]]) %*%
                                           (num.mat[[r]]*as.vector(num.exb[[r]])))) -
          as.vector(oldlambda[[r]]) * pmin(dJ[[r]], 0) +
          1e-12
        den.mat[[r]] = Reduce("rbind", list(if(length(tr)!=0) (tr.PSI[[r]] *
                                                                 as.vector(llparms$tr.S) *
                                                                 as.vector(llparms$zr.pi)) /
                                              as.vector(llparms$tr.S.pop), Reduce("rbind",
                                                                                  mapply(function(a,b) if(a!=0) b[[r]],
                                                                                         n.e, te.PSI.qr, SIMPLIFY = FALSE)),
                                            Reduce("rbind", mapply(function(a,b,c) if(a!=0)
                                            {b[[r]] / c}, n.i, thetaparms$ti.B2.qr, llparms$ti.F.q,
                                            SIMPLIFY = FALSE))))
        den.exb[[r]] = Reduce("rbind",  list(if(length(tr)!=0)
          llparms$xr.exb.q[[r]], Reduce("rbind",
                                        mapply(function(a,b) if(a!=0) b[[r]], n.e,
                                               llparms$xe.exb.qr, SIMPLIFY = FALSE)),  Reduce ("rbind",
                                                                                               mapply(function(a,b) if(a!=0) b[[r]], n.i,
                                                                                                      llparms$xi.exb.qr, SIMPLIFY = FALSE))))
        thetascore.den[[r]] = as.vector(t(theta.den.1[[r]]) %*%
                                          (den.mat[[r]] * as.vector(den.exb[[r]]))) +
          as.vector(oldlambda[[r]]) * pmax(dJ[[r]], 0) +
          1e-12
        thetascore.half[[r]] = oldtheta[[r]]*(thetascore.num[[r]] /
                                                thetascore.den[[r]])
        thetascore[[r]] = thetascore.num[[r]] - thetascore.den[[r]]
        newtheta[[r]] = as.vector(oldtheta[[r]] + (thetascore.half[[r]] -
                                                     oldtheta[[r]]))
        base = updbase(newtheta, tr.PSI, te.psi, ti.rml.gq.w.psi,
                       te.PSI.qr, ti.gq.PSI.qr)
        llparms <- updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                              base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                              oldgamm, ze, zi, zr)
        R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, newtheta, R.mat,
                       oldlambda, SIMPLIFY = FALSE)
        llik3 =  sum(log(unlist(llparms$te.h.q)+1e-12),
                     -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                     -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                     log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)),
                     na.rm = TRUE)
        ome = 0.6
        while(llik3 <= llik2){
          newtheta[[r]] = as.vector(oldtheta[[r]] + ome *
                                      (thetascore.half[[r]] - oldtheta[[r]]))
          base = updbase(newtheta, tr.PSI, te.psi, ti.rml.gq.w.psi,
                         te.PSI.qr, ti.gq.PSI.qr)
          llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                               base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                               oldgamm, ze, zi, zr)
          R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, newtheta, R.mat,
                         oldlambda, SIMPLIFY = FALSE)
          llik3 =  sum(log(unlist(llparms$te.h.q)+1e-12),
                       -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                       -sum(unlist(R.pen)), log(unlist(llparms$ze.pi)),
                       log(unlist(llparms$zi.pi)), log(unlist(llparms$tr.S.pop)),
                       na.rm = TRUE)
          if (ome >= 1e-2) ome = ome * 0.6
          else if (ome < 1e-2 & ome >= 1e-5) ome = ome * 5e-2
          else if (ome < 1e-5 & ome > 1e-20) ome = ome * 1e-5
          else break
        }
        llik2 = llik3
        oldtheta = newtheta
      }
      if(((max(mapply(function(a,b) abs(a - b), oldbeta, prev.oldbeta)) <
           control$inner.conv) & (max(mapply(function(a,b) max(abs(a - b)),
                                             oldtheta, prev.oldtheta)) < control$inner.conv) &
          (max(mapply(function(a,b) max(abs(a - b)), oldgamm, prev.oldgamm)) < control$inner.conv))){
        break
      }
    }#inner
    base <- updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
                    ti.gq.PSI.qr)
    llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                         base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                         oldgamm, ze, zi, zr)
    betaparms <- updscorebeta(llparms$ti.h.q, llparms$ti.S.gq.q,
                              base$ti.gq.H0.qr, llparms$xi.exb.qr, ti.rml.gq.w)
    thetaparms <- updscoretheta(oldtheta, llparms$ti.S.gq.q, ti.gq.psi.qr,
                                ti.rml.gq.w.l, base$ti.h0.q.l, llparms$xi.exb.qr)
    ze.pi.all <- Reduce("rbind", llparms$ze.pi)
    zi.pi.all <- Reduce("rbind", llparms$zi.pi)
    zr.pi.all <- llparms$zr.pi
    tr.S <- llparms$tr.S
    tr.H.q <- llparms$tr.H.q

    hessupd = updhessian(ze.all, ze.pi.all, zr, tr.S, zr.pi.all, zi.all.t,
                         zi.all, tr.H.q, xr, tr.PSI, llparms$xr.exb, llparms$ti.S.gq.q, llparms$ti.h.q, gq.points,
                         n.i, base$ti.gq.H0.qr, llparms$xi.exb.qr, ti.rml.gq.w,
                         ti.gq.psi.qr, ti.gq.PSI.qr, llparms$te.H.qr,
                         llparms$ti.F.q, betaparms$ti.A.qr, te.PSI.qr,
                         thetaparms$ti.B1.qr, thetaparms$ti.B2.qr, te.psi.qr, oldtheta,
                         llparms$xr.exb.q, llparms$xe.exb.qr, oldlambda, R.mat)
    hessian = hessupd$hess
    theta.gradient = active.theta = active.gradient = active = list()
    active.r = nonactive.r.index = nonactive.n = list()
    for(r in 1:n.risk){
      theta.gradient[[r]] = thetascore.num[[r]] - thetascore.den[[r]]
      active.theta[[r]] = oldtheta[[r]] < 1e-2
      active.gradient[[r]] = theta.gradient[[r]] < -0.01
      active[[r]] = active.theta[[r]]==1 & active.gradient[[r]]==1
      active.r[[r]] = which(active[[r]] %in% 1)
      nonactive.n[[r]] = n.basis[[r]] - sum(active[[r]])
      nonactive.r.index[[r]] = seq(from = if(r==1) 1
                                   else utils::tail(nonactive.r.index[[r-1]],n=1)+1,
                                   by = 1, length.out = nonactive.n[[r]])
    }
    active.index = which(unlist(active)) + n.risk*p + o
    gamma.index = seq(1, o, 1)
    beta.index = seq(o+1,(n.risk*p + o),1)
    if(length(active.index)!=0){
      G.QR.inv = solve(-hessian[-active.index, -active.index])
      pos.def = if(min(eigen(-hessian[-active.index,
                                      -active.index])$values) < 0) 0
      #pos.def = if(min(sqrt(diag(G.QR.inv))) <= 0) 0
      else 1
      neg.var = if(min(diag(G.QR.inv)) <= 0) 0
      else 1
    }
    else{
      G.QR.inv = solve(-hessian)
      pos.def = if(min(eigen(-hessian)$values) < 0) 0
      #pos.def = if(min(sqrt(diag(G.QR.inv))) <= 0) 0
      else 1
      neg.var = if(min(diag(G.QR.inv)) <= 0) 0
      else 1
    }
    if(neg.var==0) valid = 0
    if(neg.var==0 & control$aps == TRUE){
      if(control$iter.disp == TRUE){
        cat("\nWarning: Estimation terminated early at outer iteration", outer,
            "\n")
        cat("Negative variances, penalty value cannot be
         computed\n")
      }
      break
    }
    G.QR.inv.theta = G.QR.inv[c(-gamma.index, -beta.index), c(-gamma.index, -beta.index)]
    G.QR.inv.r = R.sigma = oldsigma = vr = R.sigma.r = newsigma = list()
    newdf = newlambda = list()
    for(r in 1:n.risk){
      G.QR.inv.r[[r]] = G.QR.inv.theta[nonactive.r.index[[r]],
                                       nonactive.r.index[[r]]]
      oldsigma[[r]] = 1 / (2*oldlambda[[r]])
      R.sigma[[r]] = R.mat[[r]] / as.vector(oldsigma[[r]])
      R.sigma.r[[r]] = if(sum(active[[r]])==0) R.sigma[[r]]
      else R.sigma[[r]][-active.r[[r]], -active.r[[r]]]
      vr[[r]] = sum(diag(G.QR.inv.r[[r]] %*% R.sigma.r[[r]]))
      newdf[[r]] = n.basis[[r]] - vr[[r]]
      if(newdf[[r]] < 0) newdf[[r]] = 1
      newsigma[[r]] = (t(oldtheta[[r]]) %*% R.mat[[r]] %*% oldtheta[[r]]) /
        newdf[[r]]
      newlambda[[r]] = 1 / (2*newsigma[[r]])
    }
    if(control$aps == TRUE){
      if(outer != 1){
        if(abs(max(unlist(olddf) - unlist(newdf))) < 1) break
        else{
          oldlambda = newlambda
          olddf = newdf
        }
      }
      else {
        oldlambda = newlambda
        olddf = newdf
      }
    }
    else break
  }#outer
  if(control$iter.disp == TRUE & pos.def==0){
    cat("\n")
    cat("Warning: Hessian matrix is not positive definite\n")
  }
  n.parms = n.risk*p + sum(unlist(n.basis))
  if(length(active.index)!=0){
    VarCovMat.temp = solve(-hessian[-active.index, -active.index])
    nr = nrow(VarCovMat.temp)
    VarCovMat.temp.1 = rbind(VarCovMat.temp, 0)[replace(rep(nr + 1L, nr +
                                                              length(active.index)), -active.index, seq_len(nr)), ]
    VarCovMat = cbind(VarCovMat.temp.1, 0)[, replace(rep(nr + 1L, nr +
                                                           length(active.index)), -active.index, seq_len(nr))]
  }
  else{
    VarCovMat = solve(-hessian)
  }
  se = sqrt(diag(VarCovMat))
  seB.index = seT.index = theta.index = seB = seT = seB.sand = seT.sand = list()
  for(r in 1:n.risk){
    if(r==1){
      seG.index = seq(1,o,1)
      seB.index[[r]] = seq(o+1,p+o,1)
      seT.index[[r]] = seq(o+n.risk*p+1, by=1, length.out = n.basis[[r]])
    }
    else{
      seB.index[[r]] = seq(utils::tail(seB.index[[r-1]]+1,1), by = 1,
                           length.out = p)
      seT.index[[r]] = seq(utils::tail(seT.index[[r-1]]+1,1), by = 1,
                           length.out = n.basis[[r]])
    }
    seG = se[seG.index]
    seB[[r]] = se[seB.index[[r]]]
    seT[[r]] = se[seT.index[[r]]]
  }
  fit <- list("gamma"=oldgamm, "beta" = oldbeta,"theta" = oldtheta,oldlambda,outer,iter)
  fit$z.names = znames
  fit$data = list(X = X)
  fit$n.risk = n.risk
  fit$seG = seG
  fit$seB = seB
  fit$seT = seT
  fit$gamma.gradient = gammascore
  fit$beta.gradient = betascore
  fit$theta.gradient = mapply(function(a,b) as.vector(a) - as.vector(b),
                              thetascore.num, thetascore.den, SIMPLIFY = FALSE)
  fit$VarCovMat = VarCovMat
  fit$b.knots = b.knots
  fit$i.knots = i.knots
  fit$basis.intercept = basis.intercept
  fit$dgr = dgr
  fit$theta.index = seT.index
  fit$gq.points = gq.points
  fit$nodes = gq$nodes
  fit$weights = gq$weights
  fit$pos.def = pos.def
  fit$n.basis = n.basis
  fit$valid = valid
  fit$lambda = oldlambda
  fit$n.r = n.r
  fit$n.e = n.e
  fit$n.i = n.i
  fit$hessian = hessian
  fit$cf = 1
  class(fit) = "phcsh_mpl"
  fit
}










