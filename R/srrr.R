##' Row-sparse reduced-eank regresssion
##'
##' Row-sparse reduced-rank regresssion for a prespecified rank; produce a
##' solution path for selecting predictors
##'
##' Model parameters controlling the model fitting can be specified through
##' argument \code{modstr}. The available elements include
##' \itemize{
##'     \item{lamA}: tuning parameter sequence.
##'     \item{nlam}: number of tuning parameters; no effect if \code{lamA} is
##'         specified.
##'     \item{minLambda}: minimum lambda value, no effect if \code{lamA} is
##'         specified.
##'     \item{maxLambda}: maxmum lambda value, no effect if lamA is specified.
##'     \item{WA}: adaptive weights. If NULL, the weights are constructed from
##'         RRR.
##'     \item{wgamma}: power parameter for constructing adaptive weights.
##' }
##' Similarly, the computational parameters controlling optimization can be
##' specified through argument \code{control}. The available elements include
##' \itemize{
##'     \item{epsilon}: epsilonergence tolerance.
##'     \item{maxit}: maximum number of iterations.
##'     \item{inner.eps}: used in inner loop.
##'     \item{inner.maxit}: used in inner loop.
##' }
##'
##' @usage
##' srrr(Y, X, nrank = 2, method = c("glasso", "adglasso"),
##'      ic.type = c("BIC", "BICP", "AIC", "GCV", "GIC"),
##'      A0 = NULL, V0 = NULL, modstr = list(),
##'      control = list(), screening = FALSE)
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nrank prespecified rank
##' @param method group lasso or adaptive group lasso
##' @param ic.type information criterion
##' @param A0 initial value
##' @param V0 initial value
##' @param modstr a list of model parameters controlling the model fitting
##' @param control a list of parameters for controlling the fitting process
##' @param screening If TRUE, marginal screening via glm is performed before
##'     srrr fitting.
##'
##' @return A list of fitting results
##'
##' @references
##'
##' Chen, L. and Huang, J. Z. (2012) Sparse reduced-rank regression for
##' simultaneous dimension reduction and variable selection. \emph{Journal of
##' the American Statistical Association}. 107:500, 1533--1545.
##'
##' @examples
##' library(rrpack)
##' p <- 100; n <- 100; nrank <- 3
##' mydata <- rrr.sim2(n, p, p0 = 10,q = 50, q0 = 10, nrank = 3,
##'                    s2n = 1, sigma = NULL, rho_X = 0.5, rho_E = 0)
##' fit1 <- with(mydata, srrr(Y, X, nrank = 3))
##' summary(fit1)
##' coef(fit1)
##' plot(fit1)
##' @export
srrr <- function(Y,
                 X,
                 nrank = 2,
                 method = c("glasso", "adglasso"),
                 ic.type = c("BIC", "BICP", "AIC", "GCV", "GIC"),
                 A0 = NULL,
                 V0 = NULL,
                 modstr = list(),
                 control = list(),
                 screening = FALSE)
{
  ## record function call
  Call <- match.call()

  method <- match.arg(method)
  ic <- match.arg(ic.type)

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

  if (screening) {
    ## Initial screening from lasso
    ini.screen <- glm.screen(Y, X)
    p.index <- ini.screen$p.index
    q.index <- ini.screen$q.index

    ## Reduce the problem
    Yorg <- Y
    Xorg <- X
    porg <- p
    qorg <- q

    Y <- Y[, q.index]
    X <- X[, p.index]
    q <- length(q.index)
    p <- length(p.index)
  }

  control <- do.call("srrr.control", control)
  modstr <- do.call("srrr.modstr", modstr)



  nlam <- modstr$nlam
  lamA <- modstr$lamA
  minLambda <- modstr$minLambda
  maxLambda <- modstr$maxLambda
  WA <- modstr$WA
  wgamma <- modstr$wgamma

  epsilon <- control$epsilon
  maxit <- control$maxit
  inner.eps <- control$inner.eps
  inner.maxit <- control$inner.maxit

  ICpath <- array(NA, dim = c(nlam, 5))
  Apath <- array(NA, dim = c(nlam, p, nrank))
  Vpath <- array(NA, dim = c(nlam, q, nrank))

  if (is.null(V0) | is.null(A0)) {
    ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    V0 <- ini$coefSVD$v
    A0 <- ini$coefSVD$u * ini$coefSVD$d
  } else {
    ini <- NULL
    C0 <- A0 %*% t(V0)
  }

  ## Adaptive weights
  if (is.null(WA)) {
    if (is.null(ini))
      ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    WA <- switch(
      method,
      ##lasso=matrix(nrow=p,ncol=nrank,1),
      ##adlasso=(ini$U%*%ini$D)^{-wgamma},
      glasso = rep(1, p),
      adglasso = apply(ini$coef, 1, function(a)
        sqrt(sum(a ^ 2))) ^ {
          -wgamma
        }
    )
  }

  ## Tuning sequence
  if (is.null(lamA)) {
    if (is.null(ini))
      ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    if (is.null(maxLambda))
      maxLambda <- max(apply(ini$coef, 1, function(a)
        sqrt(sum(a ^ 2))) / WA) * 2
    if (is.null(minLambda))
      minLambda <- maxLambda * 1e-4
    if (is.null(nlam))
      nlam <- 20
    lamA <- c(0, exp(seq(
      log(minLambda),
      log(maxLambda),
      length = nlam - 1
    )))
  }

  YX <- crossprod(Y, X)
  ## initial
  Apath[1, ,] <- A0
  Vpath[1, ,] <- V0

  h <- 1
  for (l in seq_len(nlam)) {
    A0 <- as.matrix(Apath[ifelse(h > 1, h - 1, 1), , ])
    V0 <- as.matrix(Vpath[ifelse(h > 1, h - 1, 1), , ])
    h <- h + 1

    ##cat("l =", l,"\n")
    fit <- srrr_Rcpp(Y,
                     X,
                     method,
                     A0,
                     V0,
                     nrank,
                     lamA[l],
                     epsilon,
                     maxit,
                     inner.eps,
                     inner.maxit,
                     WA)

    A0 <- fit$A
    V0 <- fit$V
    Apath[l, , ] <- A0
    Vpath[l, , ] <- V0
    ICpath[l, 1] <- fit$BIC
    ICpath[l, 2] <- fit$BICP
    ICpath[l, 3] <- fit$AIC
    ICpath[l, 4] <- fit$GCV
    ICpath[l, 5] <- fit$GIC
    ##print(paste("lam1=",lam1[l1],"   lam2=",lam2[l2]))
  }

  minid <- apply(ICpath, 2, which.min)

  ##Which IC to use?
  ICid <- switch(
    ic,
    BIC  = minid[1],
    BICP = minid[2],
    AIC  = minid[3],
    GCV  = minid[4],
    GIC  = minid[5]
  )
  ## Note that minid are reported based on the original tuning sequences.

  A <- as.matrix(Apath[ICid, ,])
  d <- apply(A, 2, function(a)
    sqrt(sum(a ^ 2)))
  if (screening) {
    U.reduced <- A %*% diag(1 / d, nrow = nrank, ncol = nrank)
    V.reduced <- Vpath[ICid, ,]
    coef.reduced <- A %*% t(V.reduced)
    U <- matrix(nrow = porg, ncol = nrank, 0)
    U[p.index,] <- U.reduced
    V <- matrix(nrow = qorg, ncol = nrank, 0)
    V[q.index,] <- V.reduced
    D <- diag(d, nrow = nrank, ncol = nrank)
    coefMat <- matrix(nrow = porg, ncol = qorg, 0)
    coefMat[p.index, q.index] <- coef.reduced
  } else {
    D <- diag(d, nrow = nrank, ncol = nrank)
    U <- A %*% diag(1 / d, nrow = nrank, ncol = nrank)
    V <- Vpath[ICid, ,]
    coefMat <- A %*% t(V)
    p.index <- NULL
    q.index <- NULL
  }

  out <- list(
    call = Call,
    Y = Y,
    X = X,
    p.index = p.index,
    q.index = q.index,
    A.path = Apath,
    V.path = Vpath,
    ic.path = ICpath,
    lambda = lamA,
    minid = minid,
    coef = coefMat,
    U = U,
    V = V,
    D = D,
    rank = nrank
  )
  class(out) <- "srrr"
  out
}



##' Row-sparse reduced-rank regression tuned by cross validation
##'
##' Row-sparse reduced-rank regression tuned by cross validation
##'
##' Model parameters controlling the model fitting can be specified through
##' argument \code{modstr}. The available elements include
##' \itemize{
##'     \item{lamA}: tuning parameter sequence.
##'     \item{nlam}: number of tuning parameters; no effect if \code{lamA} is
##'         specified.
##'     \item{minLambda}: minimum lambda value, no effect if \code{lamA} is
##'         specified.
##'     \item{maxLambda}: maxmum lambda value, no effect if lamA is specified.
##'     \item{WA}: adaptive weights. If NULL, the weights are constructed from
##'         RRR.
##'     \item{wgamma}: power parameter for constructing adaptive weights.
##' }
##' Similarly, the computational parameters controlling optimization can be
##' specified through argument \code{control}. The available elements include
##' \itemize{
##'     \item{epsilon}: epsilonergence tolerance.
##'     \item{maxit}: maximum number of iterations.
##'     \item{inner.eps}: used in inner loop.
##'     \item{inner.maxit}: used in inner loop.
##' }
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nrank prespecified rank
##' @param method group lasso or adaptive group lasso
##' @param nfold fold number
##' @param norder for constructing the folds
##' @param A0 initial value
##' @param V0 initial value
##' @param modstr a list of model parameters controlling the model
##'     fitting
##' @param control a list of parameters for controlling the fitting process
##' @return A list of fitting results
##' @references
##'
##' Chen, L. and Huang, J.Z. (2012) Sparse reduced-rank regression for
##' simultaneous dimension reduction and variable selection. \emph{Journal of
##' the American Statistical Association}. 107:500, 1533--1545.
##'
##' @export
cv.srrr <-
  function(Y,
           X,
           nrank = 1,
           method = c("glasso", "adglasso"),
           nfold = 5,
           norder = NULL,
           A0 = NULL,
           V0 = NULL,
           modstr = list(),
           control = list())
  {
    ## record function call
    Call <- match.call()

    method <- match.arg(method)

    p <- ncol(X)
    q <- ncol(Y)
    n <- nrow(Y)

    if (n != nrow(X))
      stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

    control <- do.call("srrr.control", control)
    modstr <- do.call("srrr.modstr", modstr)

    nlam <- modstr$nlam
    lamA <- modstr$lamA
    minLambda <- modstr$minLambda
    maxLambda <- modstr$maxLambda
    WA <- modstr$WA
    wgamma <- modstr$wgamma

    epsilon <- control$epsilon
    maxit <- control$maxit
    inner.eps <- control$inner.eps
    inner.maxit <- control$inner.maxit

    ## ICpath <- array(dim=c(nlam,4),NA)
    ## Apath <- array(dim=c(nlam,p,nrank),NA)
    ## Vpath <- array(dim=c(nlam,q,nrank),NA)

    if (is.null(V0) | is.null(A0)) {
      ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
      V0 <- ini$coefSVD$v
      A0 <- ini$coefSVD$u * ini$coefSVD$d
    } else{
      ini <- NULL
      C0 <- A0 %*% t(V0)
    }

    ## Adaptive weights
    if (is.null(WA)) {
      if (is.null(ini))
        ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
      WA <- switch(
        method,
        ## lasso=matrix(nrow=p,ncol=nrank,1),
        ## adlasso=(ini$U%*%ini$D)^{-wgamma},
        glasso = rep(1, p),
        adglasso = apply(ini$coef, 1, function(a)
          sqrt(sum(a ^ 2))) ^ {
            -wgamma
          }
      )
    }

    ## Tuning sequence
    if (is.null(lamA)) {
      if (is.null(ini))
        ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
      if (is.null(maxLambda))
        maxLambda <- max(apply(ini$coef, 1, function(a)
          sqrt(sum(a ^ 2))) / WA) * 2
      if (is.null(minLambda))
        minLambda <- maxLambda * 1e-4
      if (is.null(nlam))
        nlam <- 20
      lamA <- c(0, exp(seq(
        log(minLambda),
        log(maxLambda),
        length = nlam - 1
      )))
    }

    ndel <- round(n / nfold)
    if (is.null(norder))
      norder <- sample(n)

    cr_path <- array(NA, dim = c(nlam, nfold))

    for (f in seq_len(nfold)) {
      ##determine cr sample
      if (f != nfold) {
        iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
      } else {
        iddel <- norder[(1 + ndel * (f - 1)):n]
      }
      ndel <- length(iddel)
      nf <- n - ndel
      ## idkeep <- seq_len(n)[-iddel]

      Xf <- X[-iddel, ]
      Xfdel <- X[iddel, ]
      Yf <- Y[-iddel, ]
      Yfdel <- Y[iddel, ]

      ## YXf <- crossprod(Yf,Xf)
      fitf <- srrr(
        Yf,
        Xf,
        nrank = nrank,
        method,
        ic.type = "BIC",
        A0 = A0,
        V0 = V0,
        list(
          lamA = lamA,
          nlam = nlam,
          minLambda = minLambda,
          maxLambda = maxLambda,
          WA = WA,
          wgamma = wgamma
        ),
        list(
          epsilon = epsilon,
          maxit = maxit,
          inner.eps = inner.eps,
          inner.maxit = inner.maxit
        )
      )

      for (ll in seq_len(nlam)) {
        cr_path[ll, f] <- sum((Yfdel - Xfdel %*%
                                 fitf$A.path[ll, ,] %*%
                                 t(fitf$V.path[ll, , ])) ^ 2)
      }
    }

    index <- order(colSums(cr_path))
    ## FIXME: can rowMeans be used?
    crerr <- rowSums(cr_path[, index]) / length(index) * nfold
    minid <- which.min(crerr)

    ## Refit with all data
    fit <- srrr(
      Y,
      X,
      nrank = nrank,
      method = method,
      ic.type = "BIC",
      A0 = A0,
      V0 = V0,
      list(
        lamA = lamA,
        nlam = nlam,
        minLambda = minLambda,
        maxLambda = maxLambda,
        WA = WA,
        wgamma = wgamma
      ),
      list(
        epsilon = epsilon,
        maxit = maxit,
        inner.eps = inner.eps,
        inner.maxit = inner.maxit
      )
    )

    A <- as.matrix(fit$A.path[minid, ,])
    d <- apply(A, 2, function(a)
      sqrt(sum(a ^ 2)))
    U <- A %*% diag(1 / d, nrow = nrank, ncol = nrank)
    V <- fit$V.path[minid, , ]
    coef <- A %*% t(V)

    ## output
    out <- list(
      call = Call,
      Y = Y,
      X = X,
      A.path = fit$A.path,
      V.path = fit$V.path,
      ic.path = fit$ic.path,
      cv.path = crerr,
      cverror = crerr[minid],
      lambda = lamA,
      minid = minid,
      A = A,
      coef = coef,
      U = U,
      V = V,
      D = diag(d, nrow = nrank, ncol = nrank),
      rank = nrank
    )
    class(out) <- "cv.srrr"
    out
  }


### internal functions =========================================================

### Screening and generate initial values.
##' @importFrom stats coef
##' @importFrom glmnet cv.glmnet
glm.screen <- function(Y,
                       X,
                       offset = NULL,
                       family = list(gaussian()),
                       familygroup = rep(1, ncol(Y)),
                       intercept = FALSE,
                       standardize = FALSE)
{
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  ## 1. Use glmnet to get initial values.
  ## 2. Decompose to get initial values
  coef.lasso <- matrix(nrow = p, ncol = q, 0)
  for (i in seq_len(q)) {
    fiti <- cv.glmnet(
      X,
      Y[, i],
      offset = offset[, q],
      family = family[familygroup[i]][[1]]$family,
      intercept = intercept,
      standardize = standardize
    )
    ## summary(fiti)
    coef.lasso[, i] <- as.vector(coef(fiti))[-1]
  }
  p.index <- which(apply(coef.lasso, 1, lnorm) != 0)
  q.index <- which(apply(coef.lasso, 2, lnorm) != 0)

  list(p.index = p.index,
       q.index = q.index)
}


## Internal function for specifying model parameters
##
## a list of internal model parameters controlling the model fitting
##
## @param lamA tuning parameter sequence
## @param nlam number of tuning parameters; no effect if lamA is specified.
## @param minLambda minimum lambda value, no effect if lamA is specified.
## @param maxLambda maxmum lambda value, no effect if lamA is specified.
## @param WA adaptive weights. If NULL, the weights are constructed from RRR
## @param wgamma power parameter for constructing adaptive weights.
##
## @return a list of model parameters.
srrr.modstr <- function(lamA = NULL,
                        nlam = 100,
                        minLambda = NULL,
                        maxLambda = NULL,
                        WA = NULL,
                        wgamma = 2)
{
  list(
    lamA = lamA,
    nlam = nlam,
    minLambda = minLambda,
    maxLambda = maxLambda,
    WA = WA,
    wgamma = wgamma
  )
}



## Internal function for specifying computation parameters
##
## a list of internal computational parameters controlling optimization
##
## @param epsilon epsilonergence tolerance
## @param maxit maximum number of iterations
## @param inner.eps used in inner loop
## @param inner.maxit used in inner loop
##
## @return a list of computational parameters.
srrr.control <- function(epsilon = 1e-4,
                         maxit = 200L,
                         inner.eps = 1e-4,
                         inner.maxit = 50L)
{
  list(
    epsilon = epsilon,
    maxit = maxit,
    inner.eps = inner.eps,
    inner.maxit = inner.maxit
  )
}
