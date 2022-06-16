##' Reduced-rank regression with a sparse singular value decomposition
##'
##' Reduced-rank regression with a sparse singular value decomposition using the
##' iterative exclusive extraction algorithm.
##'
##' The model fitting can be controled through argument \code{control}.
##' The available elements include
##' \itemize{
##'     \item{maxit}: maximum number of iterations.
##'     \item{epsilon}: convergence tolerance.
##'     \item{innerMaxit}: maximum number of iterations for inner steps.
##'     \item{innerEpsilon}: convergence tolerance for inner steps.
##'     \item{nlambda}: number of tuning parameters.
##'     \item{adaptive}: if Ture, use adaptive penalization.
##'     \item{gamma0}: power parameter for constructing adaptive weights.
##'     \item{minLambda}: multiplicate factor to determine the minimum lambda.
##'     \item{niter.eea}: the number of iterations in the iterative exclusive
##' extraction algorithm.
##'     \item{df.tol}: tolerance.
##' }
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nrank integer specification of the desired rank
##' @param ic.type character specifying which information criterion to use to
##'     select the best: `BIC', `BICP', and `AIC'
##' @param orthX logical indicating if X is orthogonal, in which case a faster
##'     algorithm is used
##' @param control a list of parameters controlling the fitting process
##' @param screening If TRUE, marginal screening via glm is performed before
##'     srrr fitting.
##'
##' @return S3 \code{rssvd.path} object, a list consisting of
##'   \item{Upath}{solution path of U}
##'   \item{Vpath}{solution path of V}
##'   \item{Dpath}{solution path of D}
##'   \item{U}{estimated left singular matrix that is orthogonal}
##'   \item{V}{estimated right singular matrix that is orthogonal}
##'   \item{D}{estimated singular values such that C=UDVt}
##'   \item{rank}{estimated rank}
##' @examples
##' library(rrpack)
##' ## Simulate data from a sparse factor regression model
##' p <- 50; q <- 50; n <- 100; nrank <- 3
##' mydata <- rrr.sim1(n, p, q, nrank, s2n = 1, sigma = NULL,
##'                    rho_X = 0.5, rho_E = 0.3)
##' fit1 <- with(mydata, rssvd(Y, X, nrank = nrank + 1))
##' summary(fit1)
##' plot(fit1)
##'@references
##'
##' Chen, K., Chan, K.-S. and Stenseth, N. C. (2012) Reduced rank stochastic
##' regression with a sparse singular value decomposition.  \emph{Journal of the
##' Royal Statistical Society: Series B}, 74, 203--221.
##'
##' @export
rssvd <-
  function(Y,
           X,
           nrank,
           ic.type = c("BIC", "BICP", "AIC"),
           orthX = FALSE,
           control = list(),
           screening = FALSE) {
    Call <- match.call()
    control <- do.call("rssvd.control", control)
    tol <- control$df.tol
    ic <- match.arg(ic.type)


    ## IC <- match(ic, c("BIC", "BICP", "AIC"))
    p <- ncol(X)
    q <- ncol(Y)
    n <- nrow(Y)

    if (n != nrow(X))
      stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

    if (screening == TRUE) {
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

    ini <- NULL
    ranknow <- nrank
    iter.eea <- 1
    Ufit <- matrix(nrow = p, ncol = ranknow, 0)
    Vfit <- matrix(nrow = q, ncol = ranknow, 0)
    Dfit <- rep(0, ranknow)

    while (iter.eea <= control$niter.eea) {
      fiti <-  rssvd.eea(
        Y,
        X,
        nrank = ranknow,
        ic.type = ic,
        orthX = orthX,
        ini = ini,
        control = control,
        screening = FALSE
      )

      ## Current solution
      dorder <- order(fiti$D, decreasing = TRUE)
      Ufit[, 1:ranknow] <- fiti$U[, dorder]
      Vfit[, 1:ranknow] <- fiti$V[, dorder]
      Dfit[1:ranknow] <- fiti$D[dorder]

      ## New rank
      ranknew <- sum(fiti$D > tol)

      ## if zero solution happens
      if (is.na(ranknew) | ranknew == 0) {
        ## ranknew <- 1
        #Ufit <- matrix(nrow = p, ncol = ranknow, 0)
        #Vfit <- matrix(nrow = q, ncol = ranknow, 0)
        #Dfit <- rep(0, ranknow)
        ## out of loop
        iter.eea <- control$niter.eea
      } else {
        Ufit <- Ufit[, 1:ranknew]
        Vfit <- Vfit[, 1:ranknew]
        Dfit <- Dfit[1:ranknew]
      }

      ## if ranknew == 1, then only need to iterate once if ranknow!=ranknew.
      if (ranknew == 1 & ranknow != ranknew) {
        ini <- rssvd.init(Y, X, nrank = ranknew, control = control)
        ini <- rssvd.init(
          Y,
          X,
          nrank = ranknew,
          U0 = Ufit,
          V0 = Vfit,
          D0 = Dfit,
          Uw = ini$U0,
          Vw = ini$V0,
          control = control
        )
        ranknow <- ranknew
        ##out of loop in the next step. one more iteration
        iter.eea <- control$niter.eea - 1
      }

      if (ranknew > 1 & iter.eea < control$niter.eea) {
        ##update weights
        if (ranknow != ranknew)
          ini <- rssvd.init(Y, X, nrank = ranknew,
                            control = control)
        ini <- rssvd.init(
          Y,
          X,
          nrank = ranknew,
          U0 = Ufit,
          V0 = Vfit,
          D0 = Dfit,
          Uw = ini$U0,
          Vw = ini$V0,
          control = control
        )
        ranknow <- ranknew
      }

      iter.eea <- iter.eea + 1

    }

    if (screening) {
      U <- matrix(nrow = porg, ncol = ranknew, 0)
      U[p.index,] <- Ufit
      V <- matrix(nrow = qorg, ncol = ranknew, 0)
      V[q.index,] <- Vfit
    } else {
      U <- Ufit
      V <- Vfit
      p.index <- NULL
      q.index <- NULL
    }

    ret <- list(
      call = Call,
      Y = Y,
      X = X,
      p.index = p.index,
      q.index = q.index,
      U.path = fiti$U.path,
      V.path = fiti$V.path,
      D.path = fiti$D.path,
      ic.path = fiti$ic.path,
      U = U,
      V = V,
      D = Dfit,
      rank = ranknew
    )
    class(ret) <- "rssvd"
    ret
  }



## Not used.
surr.fit <- function(Y,
                     X,
                     lambda,
                     U = NULL,
                     V = NULL,
                     D = NULL,
                     WU = NULL,
                     WV = NULL,
                     Xtran = NULL,
                     control = list())
{
  Call <- match.call()
  ## set up control parameters
  control <- do.call("rssvd.control", control)
  innerControl <- penreg.control(epsilon = control$innerEpsilon,
                                 maxit   = control$innerMaxit)
  ## innerControl2 <- penreg.control(epsilon = control$innerEpsilon,
  ##                                 maxit = 1)

  ## dimensions
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  ## set up initial values
  if (is.null(U) || is.null(V) || is.null(WU) || is.null(WV)) {
    ini <- rssvd.init(Y, X, nrank = 1, control = control)

    WU <- ini$Wuel
    WV <- ini$Wvel
    U  <- ini$U0
    V  <- ini$V0
    D  <- 1
  }

  U <- U / WU
  V <- V / WV
  D <- D / lnorm(U, 2) / lnorm(V, 2)
  U  <- U / lnorm(U, 2)
  V <- V / lnorm(V, 2)

  ## Transform X by weight matrix
  if (is.null(Xtran))
    Xtran <- X * as.vector(WU)
  Yvec <- as.vector(Y)  ## faster than c(Y)

  ## while(j<=max.iter && conv[j]>=converge && sum(U==0)!=p){
  flag <- FALSE
  for (j in seq_len(control$maxit)) {
    U_p <- as.matrix(U)
    V_p <- as.matrix(V)
    D_p <- D

    ll0 <- lnorm(tcrossprod(U, V), 2)
    ## U-step, adaptive lasso
    X_U <- kron_RcppArma(as.matrix(WV * V), Xtran)
    ## U <- as.matrix(.Call("pcd", X_U, Yvec, lambda,
    ##                      as.vector(U_p), inner.conv,
    ##                      as.integer(inner.iter)))
    U <-
      penreg_Rcpp(Yvec, X_U, lambda, as.vector(U_p), innerControl)
    D <- lnorm(U)
    ## FIXME: what is all(U == 0)?
    if (D == 0 || is.nan(D)) {
      flag <- TRUE
      conv <- control$epsilon * 2
      break
    }
    U <- U / D

    ## V-step, adaptive lasso
    X_V <- kron_RcppArma(diag(as.vector(WV)), Xtran %*% U)
    ## V <- as.matrix(.Call("pcd", X_V, Yvec, lambda,
    ##                      as.vector(V_p), inner.conv,
    ##                      as.integer(inner.iter)))
    V <-
      penreg_Rcpp(Yvec, X_V, lambda, as.vector(V_p), innerControl)
    D <- lnorm(V)
    if (D == 0 || is.nan(D)) {
      flag <- TRUE
      conv <- control$epsilon * 2
      break
    }
    V <- V / D

    ## Convergence check
    conv <-
      lnorm(tcrossprod(U, V) - tcrossprod(U_p, V_p), 2) / ll0
    if (conv < control$epsilon || all(U == 0))
      break
  }

  if (!flag) {
    residual <- Yvec - X_V %*% V * D # D is a scalar
    sse <- c(crossprod(residual, residual))
    ## incorrect definition
    ## obj <- sse/2/q/n + lambda*D

    converged <- ifelse(conv <= control$epsilon, TRUE, FALSE)

    U <- round(U, 8)
    V <- round(V, 8)
    dfu0 <- sum(U != 0)
    dfv0 <- sum(V != 0)

    V <- WV * V * D
    U <- WU * U
    DV <- lnorm(V, 2)
    DU <- lnorm(U, 2)
    V <- V / DV
    U <- U / DU
    D <- DV * DU

    BIC <-  log(sse) + log(q  *  n) / q / n * (dfv0 + dfu0 - 1)
    BICP <- log(sse) + 2 * log(p * q) / q / n * (dfv0 + dfu0 - 1)
    AIC <-  log(sse) + 2 / q / n * (dfv0 + dfu0 - 1)

  } else {
    residual <- Yvec
    sse <- c(crossprod(residual, residual))
    ## obj <- sse/2/q/n + lambda*D

    converged <- ifelse(conv <= control$epsilon, TRUE, FALSE)

    U <- as.vector(rep(0, p))
    V <- as.vector(rep(0, q))
    D <- 0
    dfu0 <- 0
    dfv0 <- 0

    BIC <-  log(sse)
    BICP <- log(sse)
    AIC <-  log(sse)
  }

  ret <- list(
    call = Call,
    sse = sse,
    df = c(dfu0, dfv0),
    ic = c(BIC, BICP, AIC),
    U = U,
    V = V,
    D = D,
    converged = converged
  )
  class(ret) <- "surr"
  ret
}



rssvd.eea <- function(Y,
                      X,
                      nrank,
                      minLambda = 1e-04,
                      ic.type = c("BICP", "BIC", "AIC"),
                      orthX = FALSE,
                      ini = NULL,
                      control = list(),
                      screening = FALSE)
{
  Call <- match.call()
  control <- do.call("rssvd.control", control)
  ic <- match.arg(ic.type)
  IC <- match(ic, c("BICP", "BIC", "AIC"))

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

  if (screening == TRUE) {
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

  if (is.null(ini)) {
    if (orthX) {
      ini <- rssvd.init.orth(Y, X, nrank = nrank, control = control)
    } else{
      ini <- rssvd.init(Y, X, nrank = nrank, control = control)
    }
  }
  ## Information criteria
  ## BIC, BICP, AIC
  ## Default is BIC
  nlambda <- control$nlambda
  IC_path <- array(dim = c(nlambda, nrank, 3))

  ## Solution path
  U_path <- array(dim = c(p, nlambda, nrank))
  V_path <- array(dim = c(q, nlambda, nrank))
  D_path <- matrix(nrow = nlambda, ncol = nrank)

  ## Final solution selected by IC
  U <- matrix(nrow = p, ncol = nrank)
  V <- matrix(nrow = q, ncol = nrank)
  D <- rep(NA, nrank)
  Lam <- rep(NA, nrank)

  for (k in 1:nrank) {
    lamSca <- log(ini$lamMax[k]) + log(c(0.9, control$minLambda))
    lambda <-
      c(exp(seq(lamSca[1], lamSca[2], length = nlambda - 1)), 0)

    UU <- matrix(ncol = nlambda, nrow = p)
    VV <- matrix(ncol = nlambda, nrow = q)
    DD <- rep(NA, nlambda)

    ##inital nonzero values
    Uini <- ini$Uini[, k]
    Vini <- ini$Vini[, k]
    Dini <- 1

    WU <- ini$Wuel[, k]
    WV <- ini$Wvel[, k]
    Yel <- ini$Yel[, , k]
    ## This line is wrong!!!
    ## Xtran <- X * WU ## WU is a vector now, extracted from that column
    Xtran <- X %*% diag(WU)

    ## cat("check")
    ## fitting
    for (i in 1:nlambda) {
      fit <- if (orthX) {
        surr.orthX.fit(Yel,
                       X,
                       lambda[i],
                       Uini,
                       Vini,
                       ini$Wuel[, k],
                       ini$Wvel[, k],
                       control = control)
      } else {
        ## if (!cpp) {
        ##      surr.fit(Yel, X, lambda[i], Uini, Vini, Dini,
        ##               WU, WV, Xtran, control)
        ## } else {
        surr_fit_Rcpp(Yel, X, lambda[i],
                      Uini, Vini, WU, WV,
                      Xtran, control)
        ## }
      }
      ## If converges

      ## if (fit$converged) {
      ## CONV[i] <- fit$converged
      ## BIC[i,k] <- fit$IC[1]
      ## BICP[i,k] <- fit$IC[2]
      ## AIC[i,k] <- fit$IC[3]
      IC_path[i, k, ] <- fit$ic

      ## SSE[i] <- fit$sse
      VV[, i] <- fit$V
      UU[, i] <- fit$U
      DD[i] <- fit$D
      if (fit$D != 0 & !is.nan(fit$D)) {
        Uini <- fit$U
        Vini <- fit$V
        Dini <- fit$D
      }
      ## U_star <- fit$U
    }

    logSY <- log(sum(Yel ^ 2))
    ## choose lambda
    minid <- apply(IC_path[, k, ], 2, which.min)

    ## Final solution based on IC
    if (logSY <= IC_path[minid[IC], k, IC]) {
      U[, k] <- 0
      V[, k] <- 0
      D[k] <- 0
      Lam[k] <- ini$lamMax[k]
    } else {
      U[, k] <- UU[, minid[IC]]
      V[, k] <- VV[, minid[IC]]
      D[k] <- DD[minid[IC]]
      Lam[k] <- lambda[minid[IC]]
    }

    U_path[, , k] <- UU
    V_path[, , k] <- VV
    D_path[, k] <- DD
  }

  if (screening) {
    U.reduced <- U
    V.reduced <- V
    U <- matrix(nrow = porg, ncol = nrank, 0)
    U[p.index, ] <- U.reduced
    V <- matrix(nrow = qorg, ncol = nrank, 0)
    V[q.index, ] <- V.reduced
  } else {
    p.index <- NULL
    q.index <- NULL
  }

  ret <- list(
    call = Call,
    p.index = p.index,
    q.index = q.index,
    U.path = U_path,
    V.path = V_path,
    D.path = D_path,
    ic.path = IC_path,
    U = U,
    V = V,
    D = D
  )
  class(ret) <- "rssvd"
  ret
}


### internal functions =========================================================

## Internal function for specifying computation parameters
##
## a list of internal computational parameters controlling optimization
##
## @param maxit maximum number of iterations
## @param epsilon convergence tolerance
## @param innerMaxit maximum number of iterations for inner steps
## @param innerEpsilon convergence tolerance for inner steps
## @param nlambda number of lambda values
## @param adaptive If Ture, use adaptive penalization
## @param gamma0 power parameter for constructing adaptive weights
## @param minLambda multiplicate factor to determine the minimum lambda
## @param niter.eea the number of iterations in the iterative exclusive
##     extraction algorithm
## @param df.tol tolerance
##
## @return a list of computational parameters.
rssvd.control <- function(maxit = 100L,
                          epsilon = 1e-4,
                          innerMaxit = 50L,
                          innerEpsilon = 1e-4,
                          nlambda = 100,
                          adaptive = TRUE,
                          gamma0 = 2,
                          minLambda = 1e-5,
                          niter.eea = 1,
                          df.tol = 1e-8) {
  maxit <- as.integer(maxit)
  innerMaxit <- as.integer(innerMaxit)
  nlambda <- as.integer(nlambda)
  list(
    maxit = maxit,
    epsilon = epsilon,
    innerMaxit = innerMaxit,
    innerEpsilon = innerEpsilon,
    nlambda = nlambda,
    adaptive = adaptive,
    gamma0 = gamma0,
    minLambda = minLambda,
    niter.eea = niter.eea,
    df.tol = df.tol
  )
}


## Initialization for EEA
rssvd.init <- function(Y,
                       X,
                       U0 = NULL,
                       V0 = NULL,
                       D0 = NULL,
                       Uw = U0,
                       Vw = V0,
                       nrank = 1,
                       control = list()) {
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (is.null(U0) || is.null(V0) || is.null(D0)) {
    ## Calculate LS estimator and RR estimator
    ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    C_ls <- ini$coef.ls
    C_rr <- ini$coef
    U0 <- ini$coefSVD$u
    V0 <- ini$coefSVD$v
    D0 <- ini$coefSVD$d
  } else {
    ## D0 <- diag(D0, nrow = nrank)
    ## C_rr <- U0 %*% D0 %*% t(V0)
    C_rr <- U0 %*% (D0 * t(V0))
    C_ls <- NULL
  }

  if (is.null(Uw))
    Uw <- U0
  if (is.null(Vw))
    Vw <- V0


  ## Compute exclusive layers
  Yel <- array(dim = c(n, q, nrank))
  Wuel <- matrix(1, p, nrank)
  Wvel <- matrix(1, q, nrank)

  ## Initial nonzero singular vectors
  Uini <- matrix(0, p, nrank)
  Vini <- matrix(0, q, nrank)

  control <- do.call("rssvd.control", control)
  ## Maximum lambda
  lamMax <- double(nrank)
  tol <- control$df.tol
  gamma0 <- control$gamma0
  for (k in 1:nrank) {
    Ck <- tcrossprod(U0[, k], V0[, k]) * D0[k]
    Cmk <- C_rr - Ck
    Yel[, , k] <- Y - X %*% Cmk

    if (control$adaptive) {
      if (all(Uw[, k] == 0)) {
        Wuel[, k] <- 1
      } else {
        Wuel[, k] <- abs(Uw[, k]) ^ gamma0
        Wuel[, k][Wuel[, k] == 0] <- tol
      }
      if (all(Vw[, k] == 0)) {
        Wvel[, k] <- 1
      } else {
        Wvel[, k] <- abs(Vw[, k]) ^ gamma0
        Wvel[, k][Wvel[, k] == 0] <- tol
      }
    }

    ww <- outer(Wuel[, k], Wvel[, k])
    xyk <-  crossprod(X, Yel[, , k])
    lam <- abs(ww * xyk)

    imax <- which.max(lam)
    b <- round(imax / p - 0.50001) + 1
    a <- imax - (b - 1) * p
    lamMax[k] <- lam[a, b]
    Uini[as.integer(a), k] <- 1
    Vini[as.integer(b), k] <- 1
  }

  list(
    Yel = Yel,
    C_ls = C_ls,
    C_rr = C_rr,
    Wuel = Wuel,
    Wvel = Wvel,
    U0 = U0,
    V0 = V0,
    D0 = D0,
    Uini = Uini,
    Vini = Vini,
    lamMax = lamMax,
    Mdim = c(p, q, n)
  )
}


## Initialization for EEA
rssvd.init.orth <- function(Y,
                            X,
                            U0 = NULL,
                            V0 = NULL,
                            D0 = NULL,
                            Uw = U0,
                            Vw = V0,
                            nrank = 1,
                            control = list()) {
  XY <- crossprod(X, Y)
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if (is.null(U0) || is.null(V0) || is.null(D0)) {
    ## Calculate LS estimator and RR estimator
    ##ini <- rrr.fit(Y, X, nrank=nrank, coefSVD = TRUE)
    ini <- svd(XY, nu = nrank, nv = nrank)
    ##C_ls <- ini$coef.ls
    ##C_rr <- ini$coef
    U0 <- ini$u
    V0 <- ini$v
    D0 <- ini$d[1:nrank]
    ## FIXME:
    ## C_rr <- U0 %*% diag(D0) %*% t(V0)
    C_rr <- U0 %*% (D0 * t(V0))
    C_ls <- XY

  } else {
    ## D0 <- diag(D0, nrow = nrank)
    ## C_rr <- U0 %*% D0 %*% t(V0)
    C_rr <- U0 %*% (D0 * t(V0))
    C_ls <- NULL
  }

  if (is.null(Uw))
    Uw <- U0
  if (is.null(Vw))
    Vw <- V0


  ## Compute exclusive layers
  Yel <- array(dim = c(n, q, nrank))
  Wuel <- matrix(1, p, nrank)
  Wvel <- matrix(1, q, nrank)

  ## Initial nonzero singular vectors
  Uini <- matrix(0, p, nrank)
  Vini <- matrix(0, q, nrank)

  control <- do.call("rssvd.control", control)
  ## Maximum lambda
  lamMax <- double(nrank)
  tol <- control$df.tol
  gamma0 <- control$gamma0
  for (k in 1:nrank) {
    Ck <- tcrossprod(U0[, k], V0[, k]) * D0[k]
    Cmk <- C_rr - Ck
    Yel[, , k] <- Y - X %*% Cmk

    if (control$adaptive) {
      if (all(Uw[, k] == 0)) {
        Wuel[, k] <- 1
      } else {
        Wuel[, k] <- abs(Uw[, k]) ^ gamma0
        Wuel[, k][Wuel[, k] == 0] <- tol
      }
      if (all(Vw[, k] == 0)) {
        Wvel[, k] <- 1
      } else {
        Wvel[, k] <- abs(Vw[, k]) ^ gamma0
        Wvel[, k][Wvel[, k] == 0] <- tol
      }
    }

    ww <- outer(Wuel[, k], Wvel[, k])
    xyk <-  crossprod(X, Yel[, , k])
    lam <- abs(ww * xyk)

    ##    lam <- matrix(0, p, q)
    ##     for (a in 1:p){
    ##       for (b in 1:q){
    ##         lam[a,b] <- abs(Wuel[a,k] * Wvel[b,k] *
    ##                         t(X[,a]) %*% Yel[,b,k])
    ##       }
    ##     }

    imax <- which.max(lam)
    b <- round(imax / p - 0.50001) + 1
    a <- imax - (b - 1) * p
    lamMax[k] <- lam[a, b]
    Uini[as.integer(a), k] <- 1
    Vini[as.integer(b), k] <- 1
  }

  list(
    Yel = Yel,
    C_ls = C_ls,
    C_rr = C_rr,
    Wuel = Wuel,
    Wvel = Wvel,
    U0 = U0,
    V0 = V0,
    D0 = D0,
    Uini = Uini,
    Vini = Vini,
    lamMax = lamMax,
    Mdim = c(p, q, n)
  )
}


## FIXME: is Wd needed?
surr.orthX.fit <- function(Y,
                           X,
                           lambda = 0,
                           U0 = NULL,
                           V0 = NULL,
                           Wu = NULL,
                           Wv = NULL,
                           Wd = 1.,
                           YX = NULL,
                           control = list()) {
  if (is.null(YX))
    YX <- crossprod(Y, X)
  control <- do.call("rssvd.control", control)
  innerControl <- penreg.control(epsilon = control$innerEpsilon,
                                 maxit   = control$innerMaxit)

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  if (is.null(U0) || is.null(V0) || is.null(Wu) || is.null(Wv)) {
    ini <- rssvd.init(Y, X, nrank = 1, control = control)
    U0 <- ini$Uini
    V0 <- ini$Vini
    Wu <- ini$Wuel
    Wv <- ini$Wvel
  }

  result <- .Call(
    "rssvd_orth",
    YX,
    1 / Wv,
    1 / Wu,
    Wd,
    V0,
    lambda,
    innerControl$epsilon,
    innerControl$maxit,
    PACKAGE = "rrpack"
  )
  V <- as.matrix(result[1:q])
  U <- as.matrix(result[(q + 1):(p + q)])
  D <- result[p + q + 1]

  ## calculate BIC
  residual <- Y - tcrossprod(X %*% U, V) * D

  sse <- sum(residual ^ 2)
  logsse <- log(sse)
  U <- round(U, 8)
  V <- round(V, 8)
  dfu0 <- sum(U != 0)
  dfv0 <- sum(V != 0)
  dfmodel <- dfu0 + dfv0 - 1
  BIC <- logsse + log(q * n) / q / n * dfmodel
  BICP <- logsse + 2 * log(p * q) / q / n * dfmodel
  AIC <- logsse + 2 / q / n * dfmodel

  ## FIXME: the output need to match that of surr.fit
  list(
    sse = sse,
    df = c(dfu0, dfv0),
    ic = c(BIC, BICP, AIC),
    U = U,
    V = V,
    D = D
  )
}
