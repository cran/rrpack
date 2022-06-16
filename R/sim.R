##' Simulation model 1
##'
##' Similar to the the RSSVD simulation model in Chen, Chan, Stenseth (2012),
##' JRSSB.
##'
##' @param n,p,q model dimensions
##' @param nrank model rank
##' @param s2n signal to noise ratio
##' @param sigma error variance. If specfied, then s2n has no effect
##' @param rho_X correlation parameter in the generation of predictors
##' @param rho_E correlation parameter in the generation of random errors
##'
##' @return similated model and data
##'
##'@references
##'
##' Chen, K., Chan, K.-S. and Stenseth, N. C. (2012) Reduced rank stochastic
##' regression with a sparse singular value decomposition.  \emph{Journal of the
##' Royal Statistical Society: Series B}, 74, 203--221.
##'
##' @importFrom stats runif rnorm arima.sim
##' @export
rrr.sim1 <-
  function(n = 50,
           p = 25,
           q = 25,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0)
  {
    Sigma <- CorrAR
    
    U <- matrix(ncol = nrank, nrow = p)
    V <- matrix(ncol = nrank, nrow = q)
    
    U[, 1] <- c(sample(c(-1, 1), 5, replace = T), rep(0, p - 5))
    U[, 2] <-
      c(rep(0, 3), U[4, 1], -U[5, 1], sample(c(-1, 1), 3, replace = T), rep(0, p -
                                                                              8))
    U[, 3] <- c(rep(0, 8), sample(c(-1, 1), 2, replace = T), rep(0, p - 10))
    U[, 1] <- U[, 1] / lnorm(U[, 1], 2)
    U[, 2] <- U[, 2] / lnorm(U[, 2], 2)
    U[, 3] <- U[, 3] / lnorm(U[, 3], 2)
    
    V[, 1] <- c(sample(c(1, -1), 5, replace = T) * runif(5, 0.5, 1), rep(0, q -
                                                                           5))
    V[, 2] <-
      c(rep(0, 5),
        sample(c(1, -1), 5, replace = T) * runif(5, 0.5, 1),
        rep(0, q - 10))
    V[, 3] <-
      c(rep(0, 10),
        sample(c(1, -1), 5, replace = T) * runif(5, 0.5, 1),
        rep(0, q - 15))
    V[, 1] <- V[, 1] / lnorm(V[, 1], 2)
    V[, 2] <- V[, 2] / lnorm(V[, 2], 2)
    V[, 3] <- V[, 3] / lnorm(V[, 3], 2)
    
    D <- diag(c(20, 15, 10))
    
    Xsigma <- Sigma(p, rho_X)
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    ## X <- diag(p)
    
    UU <- matrix(nrow = n, ncol = q, rnorm(n * q, 0, 1))
    if (rho_E != 0) {
      for (t in 1:n)
        UU[t, ] <- arima.sim(list(order = c(1, 0, 0), ar = rho_E), n = q)
    }
    
    C <- U %*% D %*% t(V)
    C3 <- U[, 3] %*% t(V[, 3]) * D[3, 3]
    
    Y3 <- X %*% C3
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C + UU
    
    list(
      Y = Y,
      X = X,
      C = C,
      U = U,
      V = V,
      D = D,
      Xsigma = Xsigma
    )
    
  }



##' Simulation model 2
##'
##' Similar to the the SRRR simulation model in Chen and Huang (2012), JASA
##'
##' @param n sample size
##' @param p number of predictors
##' @param p0 number of relevant predictors
##' @param q number of responses
##' @param q0 number of relevant responses
##' @param nrank model rank
##' @param s2n signal to noise ratio
##' @param sigma error variance. If specfied, then s2n has no effect
##' @param rho_X correlation parameter in the generation of predictors
##' @param rho_E correlation parameter in the generation of random errors
##'
##' @return similated model and data
##'
##' @references
##'
##' Chen, L. and Huang, J.Z. (2012) Sparse reduced-rank regression for
##' simultaneous dimension reduction and variable selection. \emph{Journal of
##' the American Statistical Association}, 107:500, 1533--1545.
##'
##' @importFrom stats rnorm
##' @export
rrr.sim2 <-
  function(n = 100,
           p = 50,
           p0 = 10,
           q = 50,
           q0 = 10,
           nrank = 3,
           s2n = 1,
           sigma = NULL,
           rho_X = 0.5,
           rho_E = 0) {
    Sigma = CorrCS
    
    A1 <- matrix(ncol = nrank, nrow = q0, rnorm(q0 * nrank))
    A0 <- matrix(ncol = nrank, nrow = q - q0, 0)
    A <- rbind(A1, A0)
    B1 <- matrix(ncol = nrank, nrow = p0, rnorm(p0 * nrank))
    B0 <- matrix(ncol = nrank, nrow = p - p0, 0)
    B <- rbind(B1, B0)
    C <- B %*% t(A)
    
    Xsigma <- Sigma(p, rho_X)
    X <- MASS::mvrnorm(n, rep(0, p), Xsigma)
    
    UU <- MASS::mvrnorm(n, rep(0, q), Sigma(q, rho_E))
    
    #    ###Their definition is tr(CXC)/tr(E), which seems to be wrong
    #    sigma <- sqrt(sum(diag(t(C)%*%Sigma(p,rho_X)%*%C))/sum(diag(t(UU)%*%UU))/s2n)
    #    UU <- UU*sigma
    
    svdC <- svd(C)
    C3 <- svdC$u[, nrank] %*% t(svdC$v[, nrank]) * svdC$d[nrank]
    Y3 <- X %*% C3
    ##sigma <- sqrt(var(as.numeric(Y3))/var(as.numeric(UU))/s2n)
    ##the same
    if (is.null(sigma)) {
      sigma <- sqrt(sum(as.numeric(Y3) ^ 2) / sum(as.numeric(UU) ^ 2) / s2n)
    }
    UU <- UU * sigma
    
    Y <- matrix(nrow = n, ncol = q, NA)
    Y <- X %*% C + UU
    
    list(
      Y = Y,
      X = X,
      C = C,
      A = A,
      B = B,
      U = B,
      V = A,
      sigma = sigma,
      Xsigma = Xsigma
    )
    
  }




##' Simulation model 3
##'
##' Generate data from a mixed-response reduced-rank regression model
##'
##' @param n sample size
##' @param p number of predictors
##' @param q.mix numbers of Gaussian, Bernolli and Poisson responses
##' @param nrank model rank
##' @param intercept a vector of intercept
##' @param mis.prop missing proportion
##' @return similated model and data
##'
##' @references
##' Chen, K., Luo, C., and Liang, J. (2017) Leveraging mixed and incomplete outcomes
##' through a mixed-response reduced-rank regression. \emph{Technical report}.
##' @importFrom stats rnorm runif gaussian binomial poisson rbinom rpois
##' @export
rrr.sim3 <-
  function(n = 100,
           p = 30,
           q.mix = c(5, 20, 5),
           nrank = 2,
           intercept = rep(0.5, 30),
           mis.prop = 0.2) {
    q <- sum(q.mix)
    ##Gaussian, binomial, poisson
    q1 <- q.mix[1]
    q2 <- q.mix[2]
    q3 <- q.mix[3]
    X <- matrix(rnorm(n * p), n, p)
    
    u0 <- matrix(nrow = p, ncol = nrank, rnorm(p * nrank, 0, 1))
    u0 <- qr.Q(qr(u0))
    v0 <-
      matrix(
        nrow = q,
        ncol = nrank,
        runif(q * nrank, 0.5, 1) * sample(x = c(1, -1), q * nrank, replace = TRUE)
      )
    C <- u0 %*% t(v0)
    ##intercept
    ##intercept <- rep(0.5,q)
    C0 <- rbind(intercept, C)
    X0 <- cbind(1, X)
    
    MU <- X0 %*% C0
    family <- list(gaussian(), binomial(), poisson())
    familygroup <- c(rep(1, q1), rep(2, q2), rep(3, q3))
    cfamily <- unique(familygroup)
    nfamily <- length(cfamily)
    
    Y <- matrix(nrow = n, ncol = q, 0)
    sigma <- 1  ##sd(MU[, familygroup == 1])
    #disp.true[repnum, imis] <- sigma^2
    if (sum(familygroup == 1) > 0) {
      Y[, familygroup == 1] <-
        MU[, familygroup == 1] + matrix(nrow = n, ncol = q1, rnorm(n * q1, 0, sigma))
    }
    if (sum(familygroup == 2) > 0) {
      prob <- as.matrix(family[[2]]$linkinv(MU[, familygroup == 2]))
      Y[, familygroup == 2] <-
        apply(prob, 2, function(a)
          rbinom(n = n, size = 1, a))
    }
    if (sum(familygroup == 3) > 0) {
      prob <- as.matrix(family[[3]]$linkinv(MU[, familygroup == 3]))
      Y[, familygroup == 3] <-
        apply(prob, 2, function(a)
          rpois(n = n, lambda = a))
    }
    ###setup missing response.
    N <- n * q
    if (is.numeric(mis.prop) & mis.prop < 1 & mis.prop > 0) {
      N_mis <- round(N * mis.prop)
      index_mis <- sample(1:N, size = N_mis, replace = FALSE)
      Y_mis <- Y
      Y_mis[index_mis] <- NA
    } else{
      index_mis <- NULL
      Y_mis <- Y
    }
    
    list(
      Y = Y,
      Y.mis = Y_mis,
      index.miss = index_mis,
      X = X,
      C = C,
      family = family[q.mix != 0],
      familygroup = familygroup
    )
  }





##' Simulation model 4
##'
##' Generate data from a mean-shifted reduced-rank regression model
##'
##' @param n sample size
##' @param p number of predictors
##' @param q numbers of responses
##' @param nrank model rank
##' @param s2n signal to noise ratio
##' @param rho_X correlation parameter for predictors
##' @param rho_E correlation parameter for errors
##' @param nout number of outliers; should be smaller than n
##' @param vout control mean-shifted value of outliers
##' @param voutsd control mean-shifted magnitude of outliers
##' @param nlev number of high-leverage outliers
##' @param vlev control value of leverage
##' @param vlevsd control magnitude of leverage
##' @param SigmaX correlation structure of predictors
##' @param SigmaE correlation structure of errors
##' @return similated model and data
##'
##' @references
##' She, Y. and Chen, K. (2017) Robust reduced-rank regression. \emph{Biometrika}, 104 (3), 633--647.
##' @importFrom stats rnorm sd
##' @export
rrr.sim4 <- function(n = 100,
                     p = 12,
                     q = 8,
                     nrank = 3,
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
                     SigmaE = "CorrCS") {
  CorrAR <- function(p, rho)
  {
    Sigma <- matrix(nrow = p, ncol = p, NA)
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        Sigma[i, j] <- rho ^ (abs(i - j))
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
  
  SigmaX <- switch(SigmaX,
                   "CorrCS" = CorrCS,
                   "CorrAR" = CorrAR)
  
  SigmaE <- switch(SigmaE,
                   "CorrCS" = CorrCS,
                   "CorrAR" = CorrAR)
  
  ##simulate data
  B0 <- matrix(ncol = nrank, nrow = p, rnorm(nrank * p))
  B1 <- matrix(ncol = nrank, nrow = q, rnorm(nrank * q))
  B <- B0 %*% t(B1)
  
  X <- MASS::mvrnorm(n, rep(0, p), SigmaX(p, rho_X))
  E <- MASS::mvrnorm(n, rep(0, q), SigmaE(q, rho_E))
  
  ##set signal to noise ratio
  mindXB <- round(svd(X %*% B)$d, 2)[nrank]
  sigma <- mindXB / sqrt(sum(E ^ 2)) / s2n
  
  Ymean <- X %*% B
  Y <- Ymean + sigma * E
  
  if (nout != 0) {
    if (is.null(vout)) {
      Ysd <- apply(Ymean, 2, sd) * sample(c(-1, 1), q, replace = TRUE)
      C <- voutsd * matrix(nrow = nout,
                           ncol = q,
                           byrow = T,
                           Ysd)
      Y[1:nout, ] <-  Y[1:nout, ] + C
    } else{
      Vout <- vout * sample(c(-1, 1), q, replace = TRUE)
      C <- matrix(nrow = nout,
                  ncol = q,
                  byrow = T,
                  Vout)
      Y[1:nout, ] <- Y[1:nout, ] + C
    }
  }
  if (nlev != 0) {
    if (is.null(vlev)) {
      Xsd <- apply(X, 2, sd) * sample(c(-1, 1), p, replace = TRUE)
      Xlev <- vlevsd * matrix(nrow = nlev,
                              ncol = p,
                              byrow = T,
                              Xsd)
      X[1:nlev, ] <-  Xlev
    } else{
      Vlev <- vlev * sample(c(-1, 1), p, replace = TRUE)
      Xlev <- matrix(nrow = nlev,
                     ncol = p,
                     byrow = T,
                     Vlev)
      X[1:nlev, ] <- Xlev
    }
  }
  
  
  list(X = X,
       Y = Y,
       B = B,
       C = C)
  
}



##' Simulation model 5
##'
##' Generate data from a mean-shifted reduced-rank regression model
##'
##' @param n sample size
##' @param p number of predictors
##' @param q numbers of responses
##' @param nrank model rank
##' @param rx rank of the design matrix
##' @param s2n signal to noise ratio
##' @param rho_X correlation parameter for predictors
##' @param rho_E correlation parameter for errors
##' @param nout number of outliers; should be smaller than n
##' @param vout control mean-shifted value of outliers
##' @param voutsd control mean-shifted magnitude of outliers
##' @param nlev number of high-leverage outliers
##' @param vlev control value of leverage
##' @param vlevsd control magnitude of leverage
##' @param SigmaX correlation structure of predictors
##' @param SigmaE correlation structure of errors
##' @return similated model and data
##'
##' @references
##' She, Y. and Chen, K. (2017) Robust reduced-rank regression. \emph{Biometrika}, 104 (3), 633--647.
##' @importFrom stats sd
##' @export
rrr.sim5 <-
  function(n = 40,
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
           SigmaE = "CorrCS") {
    CorrAR <- function(p, rho)
    {
      Sigma <- matrix(nrow = p, ncol = p, NA)
      for (i in seq_len(p)) {
        for (j in seq_len(p)) {
          Sigma[i, j] <- rho ^ (abs(i - j))
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
    
    SigmaX <- switch(SigmaX,
                     "CorrCS" = CorrCS,
                     "CorrAR" = CorrAR)
    
    SigmaE <- switch(SigmaE,
                     "CorrCS" = CorrCS,
                     "CorrAR" = CorrAR)
    
    ##simulate data
    B0 <- matrix(ncol = nrank, nrow = p, rnorm(nrank * p))
    B1 <- matrix(ncol = nrank, nrow = q, rnorm(nrank * q))
    B <- B0 %*% t(B1)
    
    Xsigma <- SigmaX(p, rho_X)
    a.eig <- eigen(Xsigma)
    Xsigma.sqrt <-
      a.eig$vectors %*% diag(sqrt(a.eig$values)) %*% solve(a.eig$vectors)
    
    X0 <- matrix(ncol = rx, nrow = n, rnorm(rx * n))
    X1 <- matrix(ncol = rx, nrow = p, rnorm(rx * p))
    X <- X0 %*% t(X1) %*% Xsigma.sqrt
    
    E <- MASS::mvrnorm(n, rep(0, q), SigmaE(q, rho_E))
    
    ##set signal to noise ratio
    mindXB <- round(svd(X %*% B)$d, 2)[nrank]
    sigma <- mindXB / sqrt(sum(E ^ 2)) / s2n
    
    Ymean <- X %*% B
    Y <- Ymean + sigma * E
    
    if (nout != 0) {
      if (is.null(vout)) {
        Ysd <- apply(Ymean, 2, sd) * sample(c(-1, 1), q, replace = TRUE)
        C <- voutsd * matrix(nrow = nout,
                             ncol = q,
                             byrow = T,
                             Ysd)
        Y[1:nout, ] <-  Y[1:nout, ] + C
      } else{
        Vout <- vout * sample(c(-1, 1), q, replace = TRUE)
        C <- matrix(nrow = nout,
                    ncol = q,
                    byrow = T,
                    Vout)
        Y[1:nout, ] <- Y[1:nout, ] + C
      }
    }
    if (nlev != 0) {
      if (is.null(vlev)) {
        Xsd <- apply(X, 2, sd) * sample(c(-1, 1), p, replace = TRUE)
        Xlev <- vlevsd * matrix(nrow = nlev,
                                ncol = p,
                                byrow = T,
                                Xsd)
        X[1:nlev, ] <-  Xlev
      } else{
        Vlev <- vlev * sample(c(-1, 1), p, replace = TRUE)
        Xlev <- matrix(nrow = nlev,
                       ncol = p,
                       byrow = T,
                       Vlev)
        X[1:nlev, ] <- Xlev
        X[1:nlev, ] <- Xlev
      }
    }
    
    
    list(X = X,
         Y = Y,
         B = B,
         C = C)
    
  }
