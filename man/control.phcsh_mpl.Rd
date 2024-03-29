% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/control.phcsh_mpl.R
\name{control.phcsh_mpl}
\alias{control.phcsh_mpl}
\title{Ancillary arguments for controling phcsh_mpl fits}
\usage{
control.phcsh_mpl(
  knots.perc.limit = c(0.05, 0.95),
  gq.points = 15,
  dgr = 2,
  basis.intercept = FALSE,
  max.outer = 10L,
  max.iter = 3000L,
  lambda = NULL,
  n.basis = NULL,
  iter.disp = FALSE,
  iknots.pos = NULL,
  aps = TRUE,
  tmid = TRUE,
  inner.conv = 1e-08
)
}
\arguments{
\item{knots.perc.limit}{a vector which defines the percentile range of data
for the boundary knot and internal knots for mspline estimation}

\item{gq.points}{defines the number of nodes to used for gaussian quadrature
estimation of numerical integrals.}

\item{dgr}{number of degrees to be used in mspline estimation}

\item{basis.intercept}{indicates whether intercept is to be included for
mspline estimation}

\item{max.outer}{maximum number of iterations to perform for the marginal
likelihood method for penalty terms estimation.}

\item{max.iter}{maximum number of iterations to perform for the quasi-Newton,
multiplicative iterative algorithm}

\item{lambda}{initial values for the penalty term estimation}

\item{n.basis}{vector which contains the number of basis functions to use by
risk}

\item{iter.disp}{indicator which displays iteration count during model fit}

\item{iknots.pos}{define the percentile placement of internal knots for
mspline estimation}

\item{aps}{boolean which defines whether the automatic parameter selection
is used for penalty term selection. Default is TRUE.}

\item{tmid}{boolean which indicates whether midpoint of interval censored
observations is used for mspline estimation}

\item{inner.conv}{convergence criteria for quasi Newton - Multiplicative
iterative algorithm for model fit.}
}
\description{
This is used to set various numeric parameters controlling a Cox model fit
using phcsh_mpl. Typically it would only be used in a call to phcsh_mpl.
}
