% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/r4.R
\name{rrs.fit}
\alias{rrs.fit}
\title{Fitting reduced-rank ridge regression with given rank and shrinkage penalty}
\usage{
rrs.fit(Y, X, nrank = min(ncol(Y), ncol(X)), lambda = 1, coefSVD = FALSE)
}
\arguments{
\item{Y}{a matrix of response (n by q)}

\item{X}{a matrix of covariate (n by p)}

\item{nrank}{an integer specifying the desired rank}

\item{lambda}{tunging parameter for the ridge penalty}

\item{coefSVD}{logical indicating the need for SVD for the
coeffient matrix int the output}
}
\value{
S3 \code{rrr} object, a list consisting of
  \item{coef}{coefficient of rrs}
  \item{coef.ls}{coefficient of least square}
  \item{fitted}{fitted value of rrs}
  \item{fitted.ls}{fitted value of least square}
  \item{A}{right singular matrix}
  \item{Ad}{sigular value vector}
  \item{nrank}{rank of the fitted rrr}
}
\description{
Fitting reduced-rank ridge regression with given rank and shrinkage penalty
}
\examples{
library(rrpack)
Y <- matrix(rnorm(400), 100, 4)
X <- matrix(rnorm(800), 100, 8)
rfit <- rrs.fit(Y, X)
}
\references{
Mukherjee, A. and Zhu, J. (2011) Reduced rank ridge regression and its
kernal extensions.

Mukherjee, A., Chen, K., Wang, N. and Zhu, J. (2015) On the degrees of
freedom of reduced-rank estimators in multivariate
regression. \emph{Biometrika}, 102, 457--477.
}
