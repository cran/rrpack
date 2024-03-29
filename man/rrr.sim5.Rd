% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.R
\name{rrr.sim5}
\alias{rrr.sim5}
\title{Simulation model 5}
\usage{
rrr.sim5(
  n = 40,
  p = 100,
  q = 50,
  nrank = 5,
  rx = 10,
  s2n = 1,
  rho_X = 0,
  rho_E = 0,
  nout = 10,
  vout = NULL,
  voutsd = 2,
  nlev = 10,
  vlev = 10,
  vlevsd = NULL,
  SigmaX = "CorrCS",
  SigmaE = "CorrCS"
)
}
\arguments{
\item{n}{sample size}

\item{p}{number of predictors}

\item{q}{numbers of responses}

\item{nrank}{model rank}

\item{rx}{rank of the design matrix}

\item{s2n}{signal to noise ratio}

\item{rho_X}{correlation parameter for predictors}

\item{rho_E}{correlation parameter for errors}

\item{nout}{number of outliers; should be smaller than n}

\item{vout}{control mean-shifted value of outliers}

\item{voutsd}{control mean-shifted magnitude of outliers}

\item{nlev}{number of high-leverage outliers}

\item{vlev}{control value of leverage}

\item{vlevsd}{control magnitude of leverage}

\item{SigmaX}{correlation structure of predictors}

\item{SigmaE}{correlation structure of errors}
}
\value{
similated model and data
}
\description{
Generate data from a mean-shifted reduced-rank regression model
}
\references{
She, Y. and Chen, K. (2017) Robust reduced-rank regression. \emph{Biometrika}, 104 (3), 633--647.
}
