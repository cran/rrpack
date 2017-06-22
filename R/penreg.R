## Auxiliary for controlling penreg.fit
##
## Auxiliary function used by penreg_Rcpp to fine tune the fitting
## of penalized regression with lasso penalty.
##
## @param epsilon positive convergence tolerance epsilon; the iterations
##   converge when |new - old | / (old + 0.1) < epsilon.
## @param maxit integer giving the maximal number of iterations.
## @param trace logical indicating if tracing the objective is needed.
penreg.control <-
  function(epsilon = 1e-5,
           maxit = 1000,
           trace = FALSE) {
    if (!is.numeric(epsilon) || epsilon <= 0)
      stop("value of 'epsilon' must be > 0")
    maxit <- as.integer(maxit)
    if (!is.integer(maxit) || maxit <= 0)
      stop("maximum number of iterations must be > 0")
    list(epsilon = epsilon,
         maxit = maxit,
         trace = trace)
  }


penreg.native <- function(Y,
                          X,
                          lambda,
                          beta0,
                          X.standardized = FALSE,
                          control = list()) {
  stopifnot(NROW(Y) == NROW(X))
  stopifnot(NCOL(X) == length(beta0))
  
  control <- do.call("surr.control", control)
  
  ss <- 1
  if (!X.standardized) {
    ## no centering?
    ss <- sqrt(colSums(X ^ 2))
    X <- X / matrix(ss,
                    nrow = nrow(X),
                    ncol = ncol(X),
                    byrow = TRUE)
    vlambda <- lambda / ss
  }
  
  beta0 <- beta0 * ss
  r <- Y - X %*% beta0
  p <- NCOL(X)
  beta <- bchg <- beta0
  for (iter in 1L:control$maxit) {
    ## loop over all covariates
    for (j in 1L:p) {
      z <- sum(X[, j] * r) + beta[j]
      bj <- softThres(z, vlambda[j])
      if (control$trace)
        cat("iter = ", iter, ", j = ", j, ", bj = ", bj, "\n")
      beta[j] <- bj
      bchg[j] <- beta[j] - beta0[j]
      if (bchg[j] != 0) {
        r <- r  - bchg[j] * X[, j]
      }
    }
    ## check for convergence
    ## what if beta == 0?
    if (sqrt(sum(bchg ^ 2) / (sum(beta0 ^ 2) + 0.1)) < control$epsilon)
      break
    beta0 <- beta
  }
  ## transform beta back to original scale
  if (!X.standardized)
    beta <- beta / ss
  beta
}
