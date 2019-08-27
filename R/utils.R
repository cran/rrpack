## lnorm
lnorm <- function(a, p = 1)
{
  (sum(abs(a) ^ p)) ^ (1. / p)
}

## Vectorization a matrix by row
mstack <- function(S)
{
  as.vector(t(S))
}


## Autoregressive covariance structure
CorrAR <- function(p, rho)
{
  Sigma <- matrix(nrow = p, ncol = p, NA)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      Sigma[i,j] <- rho ^ (abs(i - j))
    }
  }
  Sigma
}


CorrCS <- function(p, rho)
{
  Sigma <- matrix(nrow = p, ncol = p, rho)
  diag(Sigma) <- 1
  Sigma
}

## positive part
positivePart <- function(x) {
    x[x < 0] <- 0
    x
}
