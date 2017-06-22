### objects from mrrr.R ========================================================
##' @export
print.mrrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.cv.mrrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$fit$rank), "\n")

  invisible(x)
}

##' @export
print.summary.mrrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Coefficients:\n")
  print(x$coef)

  cat("\nDispersion:\n")
  print(x$dispersion)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.cv.mrrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Coefficients:\n")
  print(x$fit$coef)

  cat("\nDispersion:\n")
  print(x$fit$dispersion)

  cat("\nEstimated Rank:", sprintf("%d", x$fit$rank), "\n")

  invisible(x)
}


### objects from r4.R ==========================================================
##' @export
print.r4 <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.rrs.fit <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.r4 <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Coefficients:\n")
  print(x$coef)

  cat("\nEstimated sparse mean shifts:\n")
  print(x$s)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}


### for objects from rrr.R =====================================================
##' @export
print.rrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.cv.rrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.rrr.fit <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Coefficients:\n")
  print(x$coef)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.rrr.leverage <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nLeverages:\n")
  print(x$leverage)

  cat("\nDegree of freedom:", sprintf("%d", x$df), "\n")

  invisible(x)
}

##' @export
print.rrr.cookD <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nCook's distance:\n")
  print(x$cookD)

  cat("\nDegree of freedom:", sprintf("%d", x$df), "\n")
  print(x$df)

  invisible(x)
}

##' @export
print.summary.rrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Coefficients:\n")
  print(x$coef)

  cat("\nEstimated left singular matrix:\n")
  print(x$U)

  cat("\nEstimated singular value matrix:\n")
  print(x$D)

  cat("\nEstimated right singular matrix:\n")
  print(x$V)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.cv.rrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Coefficients:\n")
  print(x$coef)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}


### objects from sofar.R =======================================================
##' @export
print.sofar <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.cv.sofar <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.sofar <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated left singular matrix:\n")
  print(x$U)

  cat("\nEstimated singular value matrix:\n")
  print(x$D)

  cat("\nEstimated right singular matrix:\n")
  print(x$V)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")
  invisible(x)
}

##' @export
print.summary.cv.sofar <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated left singular matrix:\n")
  print(x$U)

  cat("\nEstimated singular value matrix:\n")
  print(x$D)

  cat("\nEstimated right singular matrix:\n")
  print(x$V)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")
  invisible(x)
}


### objects form srrr.R ========================================================
##' @export
print.srrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.cv.srrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.srrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated left singular matrix:\n")
  print(x$U)

  cat("\nEstimated singular value matrix:\n")
  print(x$D)

  cat("\nEstimated right singular matrix:\n")
  print(x$V)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.cv.srrr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated left singular matrix:\n")
  print(x$U)

  cat("\nEstimated singular value matrix:\n")
  print(x$D)

  cat("\nEstimated right singular matrix:\n")
  print(x$V)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}


### objects from ssvd.R ========================================================
##' @export
print.surr <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nDegree of freedom:", sprintf("%d", x$df), "\n")

  invisible(x)
}

##' @export
print.rssvd <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}

##' @export
print.summary.rssvd <- function(x, ...)
{
  cat("Call:\n")
  print.default(x$call)

  cat("\nEstimated left singular matrix:\n")
  print(x$U)

  cat("\nEstimated singular value matrix:\n")
  print(x$D)

  cat("\nEstimated right singular matrix:\n")
  print(x$V)

  cat("\nEstimated Rank:", sprintf("%d", x$rank), "\n")

  invisible(x)
}
