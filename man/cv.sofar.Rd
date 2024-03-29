% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sofar.R
\name{cv.sofar}
\alias{cv.sofar}
\title{Sparse orthognal factor regression tuned by cross validation}
\usage{
cv.sofar(
  Y,
  X,
  nrank = 1,
  su = NULL,
  sv = NULL,
  nfold = 5,
  norder = NULL,
  modstr = list(),
  control = list(),
  screening = FALSE
)
}
\arguments{
\item{Y}{response matrix}

\item{X}{covariate matrix}

\item{nrank}{an integer specifying the desired rank/number of factors}

\item{su}{a scaling vector for U such that \eqn{U^{T}U = diag(s_{u})}}

\item{sv}{a scaling vector for V such that \eqn{V^{T}V = diag(s_{v})}}

\item{nfold}{number of fold; used for cv.sofar}

\item{norder}{observation orders to constrct data folds; used for cv.sofar}

\item{modstr}{a list of internal model parameters controlling the model
fitting}

\item{control}{a list of internal computation parameters controlling
optimization}

\item{screening}{If TRUE, marginal screening via lasso is performed before
sofar fitting.}
}
\description{
Sparse orthognal factor regression tuned by cross validation
}
\details{
The model parameters can be specified through argument \code{modstr}.
The available elements include
\itemize{
    \item{mu}: parameter in the augmented Lagrangian function.
    \item{mugamma}: increament of mu along iterations to speed up
        computation.
    \item{WA}: weight matrix for A.
    \item{WB}: weight matrix for B.
    \item{Wd}: weight matrix for d.
    \item{wgamma}: power parameter in constructing adaptive weights.
}
The model fitting can be controled through argument \code{control}.
The avilable elements include
\itemize{
   \item{nlam}: number of lambda triplets to be used.
   \item{lam.min.factor}: set the smallest lambda triplets as a fraction of
       the estimation lambda.max triplets.
   \item{lam.max.factor}: set the largest lambda triplets as a multiple of
       the estimation lambda.max triplets.
   \item{lam.AB.factor}: set the relative penalty level between A/B and D.
   \item{penA,penB,penD}: if TRUE, penalty is applied.
   \item{lamA}: sequence of tuning parameters for A.
   \item{lamB}: sequence of tuning parameters for B.
   \item{lamD}: sequence of tuning parameters for d.
   \item{methodA}: penalty for penalizing A.
   \item{methodB}: penalty for penalizing B.
   \item{epsilon}: convergence tolerance.
   \item{maxit}: maximum number of iterations.
   \item{innerEpsilon}: convergence tolerance for inner subroutines.
   \item{innerMaxit}: maximum number of iterations for inner subroutines.
   \item{sv.tol}: tolerance for singular values.
}
}
