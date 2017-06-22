##' Robust reduced-rank regression
##'
##' Perform robust reduced-rank regression.
##'
##' The model parameters can be controlled through argument \code{modstr}.
##' The available elements include
##' \itemize{
##'     \item{nlam}: parameter in the augmented Lagrangian function.
##'     \item{adaptive}: if TRUE, use leverage values for adaptive
##'         penalization.
##'     \item{weights}: user supplied weights for adaptive penalization.
##'     \item{minlam}: maximum proportion of outliers.
##'     \item{maxlam}: maximum proportion of good observations.
##'     \item{delid}: discarded observation indices for initial estimation.
##' }
##' The model fitting can be controlled through argument \code{control}.
##' The available elements include
##' \itemize{
##'     \item{epsilon}: convergence tolerance.
##'     \item{maxit}: maximum number of iterations.
##'     \item{qr.tol}: tolerance for qr decomposition.
##'     \item{tol}: tolerance.
##' }
##'
##' @usage
##' r4(Y, X, maxrank = min(dim(Y), dim(X)),
##'    method = c("rowl0", "rowl1", "entrywise"),
##'    Gamma = NULL, ic.type = c("AIC", "BIC", "PIC"),
##'    modstr = list(), control = list())
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param maxrank maximum rank for fitting
##' @param method  outlier detection method, either entrywise or rowwise
##' @param Gamma weighting matrix in the loss function
##' @param ic.type information criterion, AIC, BIC or PIC
##' @param modstr a list of model parameters controlling the model
##'     fitting
##' @param control a list of parameters for controlling the fitting process
##' @return a list consisting of
##'   \item{coef.path}{solutuon path of regression coefficients}
##'   \item{s.path}{solutuon path of sparse mean shifts}
##'   \item{s.norm.path}{solutuon path of the norms of sparse mean shifts}
##'   \item{ic.path}{paths of information criteria}
##'   \item{ic.smooth.path}{smoothed paths of information criteria}
##'   \item{lambda.path}{paths of the tuning parameter}
##'   \item{id.solution}{ids of the selected solutions on the path}
##'   \item{ic.best}{lowest values of the information criteria}
##'   \item{rank.best}{rank values of selected solutions}
##'   \item{coef}{estimated regression coefficients}
##'   \item{s}{estimated sparse mean shifts}
##'   \item{rank}{rank estimate}
##' @examples
##' \dontrun{
##' library(rrpack)
##' n <- 100; p <- 500; q <- 50
##' xrank <- 10; nrank <- 3; rmax <- min(n, p, q, xrank)
##' nlam <- 100; gamma <- 2
##' rho_E <- 0.3
##' rho_X <- 0.5
##' nlev <- 0
##' vlev <- 0
##' vout <- NULL
##' vlevsd <- NULL
##' nout <- 0.1 * n
##' s2n <- 1
##' voutsd <- 2
##' simdata <- rrr.sim5(n, p, q, nrank, rx = xrank, s2n = s2n,
##'                     rho_X = rho_X, rho_E = rho_E, nout = nout, vout = vout,
##'                     voutsd = voutsd,nlev = nlev,vlev=vlev,vlevsd=vlevsd)
##' Y <- simdata$Y
##' X <- simdata$X
##' fit <- r4(Y, X, maxrank = rmax,
##'                method = "rowl0", ic.type= "PIC")
##' summary(fit)
##' coef(fit)
##' which(apply(fit$s,1,function(a)sum(a^2))!=0)
##' }
##' @references
##'
##' She, Y. and Chen, K. (2017) Robust reduced-rank regression.
##' \emph{Biometrika}. In press.
##'
##' @importFrom MASS ginv
##' @importFrom stats quantile loess
##' @export
r4 <- function(Y,
               X,
               maxrank = min(dim(Y), dim(X)),
               method = c("rowl0", "rowl1", "entrywise"),
               Gamma = NULL,
               ic.type = c("AIC", "BIC", "PIC"),
               modstr = list(),
               control = list())
{
  Call <- match.call()

  n <- nrow(Y)
  q <- ncol(Y)
  p <- ncol(X)

  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

  last.min <- function(a) {
    l <- length(a)
    minid <- which.min(a)
    if (identical(minid, l)) {
      a <- a[l:1]
      minid <-
        l - which(c(a[-length(a)] - a[-1], 0) <= 0)[1] + 1
    }
    minid
  }

  control <- do.call("r4.control", control)
  modstr <- do.call("r4.modstr", modstr)

  ic <- match.arg(ic.type)
  ic <- switch(ic,
               "AIC" = 1,
               "BIC" = 2,
               "PIC" = 3)

  if (!is.null(Gamma)) {
    eigenGm <- eigen(Gamma)
    sqrtGm <- eigenGm$vectors %*%
      diag(sqrt(eigenGm$values)) %*%
      t(eigenGm$vectors)
    sqrtinvGm <- eigenGm$vectors %*%
      diag(sqrt(eigenGm$values ^ (-1))) %*%
      t(eigenGm$vectors)
    Y <- Y %*% sqrtGm
  }

  ## xrank <- sum(round(svd(X)$d,4)!=0)
  qrX <- qr(X, tol = control$qr.tol)
  ## C_ls <- qr.coef(qrX, Y)
  ## C_ls <- ifelse(is.na(C_ls), 0, C_ls) ## FIXME
  rX <- qrX$rank

  maxrank <- modstr$maxrank
  if (is.null(maxrank))
    maxrank <- min(rX, q)

  ## number of information criteria
  nIC <- 3

  ## Record the Best model for each rank and for each IC
  Barray <- array(dim = c(maxrank, nIC, p, q))
  Carray <- array(dim = c(maxrank, nIC, n, q))

  nlam <- modstr$nlam

  ## Record all the row l2 norm of the C matrix
  Cmatrix <- array(dim = c(maxrank, nlam, n))

  IC <- array(dim = c(maxrank, nlam, nIC))
  smoothIC <- array(dim = c(maxrank, nlam, nIC))
  ## Best one for each rank
  ## and for different information criteria
  Bestid <- matrix(nrow = maxrank, ncol = nIC)
  BestIC <- matrix(nrow = maxrank, ncol = nIC)

  SSE <- matrix(nrow = maxrank, ncol = nlam)
  lammatrix <- matrix(nrow = maxrank, ncol = nlam)

  ## Compute some matrices
  H0 <- ginv(crossprod(X)) %*% t(X)
  H <- X %*% H0

  ## adaptive penalty using leverage. Similar to Owen and She (2011)
  if (modstr$adaptive) {
    weights <- modstr$weights
    if (is.null(weights)) {
      ## Note that even for entrywise penalization,
      ## weights can be a n by 1 vector suppied to hardTH function
      weights <- sqrt(1 - diag(H)) + control$tol
    }
  } else {
    weights <- rep(1, n)
  }

  minlam <- modstr$minlam
  maxlam <- modstr$maxlam


  ## randomly generate some indecies
  delid <- modstr$delid
  if (is.null(delid))
    delid <- sample(seq_len(n), round(n * minlam))

  ## Fit model for each rank
  ## check lammax
  lamratio <- 1
  nr <- 1
  time <- vector()
  tj <- 1
  while (nr <= maxrank) {
    ## initial fitting
    ## determine lambda sequence
    ini <- rrr.fit(Y[-delid,], X[-delid,], nrank = nr)
    B0 <- ini$coef
    R <- Y - X %*% B0
    ## l2row  <- hardrowTH(R%*%Gamma,0)$l2

    if (method == "rowl0") {
      l2row <-  hardrowTH(R, 0)$l2
      lammax <-  max(l2row / weights) * maxlam * lamratio
      lammin <-  quantile(l2row / weights, 1 - minlam)
    } else if (method == "rowl1") {
      l2row <-  softrowTH(R, 0)$l2
      lammax <-  max(l2row / weights) * maxlam * lamratio
      lammin <-  quantile(l2row / weights, 1 - minlam)
    } else if (method == "entrywise") {
      ## Again, this division should still work for entrywise case
      lammax <-  max(R / weights) * maxlam * lamratio
      lammin <-  quantile(as.vector(R / weights), 1 - minlam)
    }
    lamseq <-  exp(seq(log(lammin), log(lammax), length = nlam))

    lammatrix[nr,] <- lamseq

    ##Store results for all the lambdas
    Clam <- array(dim = c(nlam, n, q))
    Blam <- array(dim = c(nlam, p, q))
    Clam[1, ,] <- 0
    Blam[1, ,] <- B0

    ptm <- proc.time()
    for (lamid in seq_len(nlam)) {
      C0 <- Clam[ifelse(lamid - 1 > 0, lamid - 1, 1), ,]
      B0 <- Blam[ifelse(lamid - 1 > 0, lamid - 1, 1), ,]
      iter <- 1
      diff <- 10 * control$epsilon

      while (diff >= control$epsilon & iter < control$maxit) {
        C1 <- C0
        B1 <- B0

        R <- Y - X %*% B1
        ## This is based on She(2009)
        ## Cth <- hardrowTH(C1 + (R-C1)%*%Gamma,lamseq[lamid]*weights)

        if (method == "rowl0") {
          Cth <- hardrowTH(R, lamseq[lamid] * weights)
          C0 <- Cth$C
        } else if (method == "rowl1") {
          Cth <- softrowTH(R, lamseq[lamid] * weights)
          C0 <- Cth$C
        } else if (method == "entrywise") {
          C0 <- hardThres(R, lamseq[lamid] * weights)
        }

        Y0 <- Y - C0
        ## XC <- H %*% Y0 %*% sqrtGm
        XC <- H %*% Y0

        ## SVD should be faster then eigen
        svdXC <- tryCatch(
          svd(XC, nu = nr, nv = nr),
          error = function(e)
            2
        )
        if (tryCatch(
          svdXC == 2,
          error = function(e)
            3
        ) == 3) {
          V <- svdXC$v[, seq_len(nr)]
        } else {
          eigenXC <- tryCatch(
            eigen(crossprod(XC)),
            error = function(e)
              2
          )
          if (tryCatch(
            eigenXC == 2,
            error = function(e)
              3
          ) == 3) {
            V <- eigenXC$vectors[, seq_len(nr)]
          }
        }

        ## B0 <- H0%*%Y0%*%sqrtGm%*%V%*%t(V)%*%sqrtinvGm
        B0 <- H0 %*% Y0 %*% tcrossprod(V)

        C1norm <- sqrt(sum((C1) ^ 2))
        B1norm <- sqrt(sum((B1) ^ 2))
        diff <- ifelse(C1norm == 0, 0,
                       sqrt(sum((C0 - C1) ^ 2)) / C1norm) +
          ifelse(B1norm == 0, 0,
                 sqrt(sum((B0 - B1) ^ 2)) / B1norm)
        iter <- iter + 1
      }

      Clam[lamid, ,] <- C0
      Blam[lamid, ,] <- B0

      Cmatrix[nr, lamid,] <- apply(C0 ^ 2, 1, sum)

      SSE[nr, lamid] <- sum((Y - C0 - X %*% B0) ^ 2)

      ## FIXME?
      if (method == "rowl0" | method == "rowl1") {
        nout <- sum(apply(C0 != 0, 1, sum) != 0)
        df <- (rX + q - nr) * nr + nout * q
      } else if (method == "entrywise") {
        nout <- sum(C0 != 0)
        df <- (rX + q - nr) * nr + sum(C0 != 0)
      }
      if (method == "rowl0" | method == "rowl1") {
        nout <- sum(apply(C0 != 0, 1, sum) != 0)
        IC[nr, lamid, 1] <-
          n * q * log(SSE[nr, lamid] / n / q) + 2 * df
        IC[nr, lamid, 2] <-
          n * q * log(SSE[nr, lamid] / n / q) + log(q * n) * df
        IC[nr, lamid, 3] <- n * q * log(SSE[nr, lamid] / n / q) +
          7 * df + ifelse(nout == 0, 0, 2.1 * nout * log(exp(1) *
                                                           n / nout))
      } else if (method == "entrywise") {
        nout <- sum(C0 != 0)
        IC[nr, lamid, 1] <-
          n * q * log(SSE[nr, lamid] / n / q) + 2 * df
        IC[nr, lamid, 2] <-
          n * q * log(SSE[nr, lamid] / n / q) + log(q * n) * df
        ## Need update
        IC[nr, lamid, 3] <- n * q * log(SSE[nr, lamid] / n / q) +
          7 * df + 2 * log(choose(n * q, nout))
      }
    }

    time[tj] <- (proc.time() - ptm)[3]
    tj <- tj + 1

    ## check whether lammax is too small
    if (sum(apply(Clam[nlam, , ], 1, sum) != 0) > n * 0.05) {
      lamratio <- lamratio + 1
    } else{
      ## lamratio <- 1

      ## smooth IC curves then select a model in favor of less outliers
      ## Select the best fit for the given rank
      for (i in seq_len(nIC)) {
        smoothIC[nr, , i] <- loess(IC[nr, , i] ~
                                     c(seq_len(nlam)))$fitted
      }
      Bestid[nr, ] <- apply(smoothIC[nr, , ], 2, last.min)
      ## Bestid[nr,] <- apply(smoothIC[nr,,],2,which.min)
      ## Bestid[nr,] <- apply(IC[nr,,],2,which.min)

      for (i in 1:nIC) {
        BestIC[nr, i] <- smoothIC[nr, Bestid[nr, i], i]
        ## BestIC[nr,2] <- smoothIC[nr,Bestid[nr,2],2]
        ## BestIC[nr,3] <- smoothIC[nr,Bestid[nr,3],3]
        if (is.null(Gamma)) {
          Carray[nr, i, , ] <- Clam[Bestid[nr, i], , ]
          Barray[nr, i, , ] <- Blam[Bestid[nr, i], , ]
        } else {
          Carray[nr, i, , ] <- Clam[Bestid[nr, i], , ] %*% sqrtinvGm
          Barray[nr, i, , ] <-
            Blam[Bestid[nr, i], , ] %*% sqrtinvGm
        }
        ## Barray[nr,2,,] <- Blam[Bestid[nr,2],,]
        ## Barray[nr,3,,] <- Blam[Bestid[nr,3],,]
      }
      nr <- nr + 1
    }
  }
  ## Select the best rank
  Bestrank <- apply(BestIC, 2, which.min)

  coef.b = Barray[Bestrank[ic], ic, , ]
  coef.c = Carray[Bestrank[ic], ic, , ]

  out <- list(
    call = Call,
    coef.path = Barray,
    s.path = Carray,
    s.norm.path = Cmatrix,
    ic.path = IC,
    ic.smooth.path = smoothIC,
    lambda.path = lammatrix,
    id.solution = Bestid,
    ic.best = BestIC,
    rank.best = Bestrank,
    ## time = time,
    ## sse = SSE,
    ## lamratio=lamratio
    coef = coef.b,
    s = coef.c,
    rank = Bestrank[ic]
  )
  class(out) <- "r4"
  out
}



##' Fitting reduced-rank ridge regression with given rank and shrinkage penalty
##'
##' Fitting reduced-rank ridge regression with given rank and shrinkage penalty
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param nrank an integer specifying the desired rank
##' @param lambda tunging parameter for the ridge penalty
##' @param coefSVD logical indicating the need for SVD for the
##'   coeffient matrix int the output
##' @return S3 \code{rrr} object, a list consisting of
##'   \item{coef}{coefficient of rrs}
##'   \item{coef.ls}{coefficient of least square}
##'   \item{fitted}{fitted value of rrs}
##'   \item{fitted.ls}{fitted value of least square}
##'   \item{A}{right singular matrix}
##'   \item{Ad}{sigular value vector}
##'   \item{nrank}{rank of the fitted rrr}
##' @examples
##' library(rrpack)
##' Y <- matrix(rnorm(400), 100, 4)
##' X <- matrix(rnorm(800), 100, 8)
##' rfit <- rrs.fit(Y, X)
##' @references
##'
##' Mukherjee, A. and Zhu, J. (2011) Reduced rank ridge regression and its
##' kernal extensions.
##'
##' Mukherjee, A., Chen, K., Wang, N. and Zhu, J. (2015) On the degrees of
##' freedom of reduced-rank estimators in multivariate
##' regression. \emph{Biometrika}, 102, 457--477.
##'
##' @importFrom MASS ginv
##' @export
rrs.fit <- function(Y,
                    X,
                    nrank = min(ncol(Y), ncol(X)),
                    lambda = 1,
                    coefSVD = FALSE)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  S_yx <- t(Y) %*% X
  ## This is a key difference
  S_xx <- t(X) %*% X + lambda * diag(p)

  ## S_xx_inv <- tryCatch(solve(S_xx+0.1*diag(p)),error=function(e)ginv(S_xx))
  ## S_xx_inv <- ginv(S_xx)
  ## Use the Woodbury matrix identity
  if (lambda != 0) {
    S_xx_inv <- 1 / lambda * diag(p) -
      lambda ^ (-2) * t(X) %*% ginv(diag(n) + lambda ^ (-1) * X %*%
                                      t(X)) %*% X
  } else{
    S_xx_inv <- ginv(S_xx)
    if (sum(is.na(S_xx_inv)) > 0) {
      S_xx_inv <- solve(S_xx + 0.1 * diag(p))
    }
  }

  C_ls <- S_xx_inv %*% t(S_yx)

  ypy.svd <- TRUE
  ##if(ypy.svd){
  ##This is another key difference
  XC <- rbind(X, sqrt(lambda) * diag(p)) %*% C_ls
  svdXC <- tryCatch(
    svd(XC, nu = nrank, nv = nrank),
    error = function(e)
      2)
  if (tryCatch(
    svdXC == 2,
    error = function(e)
      3) == 3) {
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
  } else{
    ypy.svd <- FALSE
  }
  #}
  if (!ypy.svd) {
    SS <- S_yx %*% C_ls
    SS <- (SS + t(SS)) / 2
    eigenSS <- eigen(SS, symmetric = TRUE)
    A <- as.matrix(eigenSS$vectors[, 1:nrank])
    Ad <- eigenSS$values[1:nrank]
  }

  AA <- A %*% t(A)
  C_rr <- C_ls %*% AA

  ##    if(c.svd){
  ##      svd_C <- svd(C_rr,nv=nrank,nu=nrank)
  ##      U <- as.matrix(svd_C$u[,1:nrank])
  ##      V <- as.matrix(svd_C$v[,1:nrank])
  ##      D <- diag(svd_C$d[1:nrank],nrow=nrank)
  ##
  ##      ####return ls estimator C_ls, reduced-rank estimator C_rr
  ##      ####return SVD of C_rr
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,U=U,V=V,D=D,C=C_rr,rank=nrank)
  ##    }else{
  ##      list(A=A,Ad=Ad,C_ls=C_ls,C_rr=C_rr,C=C_rr,rank=nrank)
  ##    }

  ret <- list(
    call = Call,
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    nrank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  class(ret) <- "rrs.fit"
  ret
}


## Internal function for specifying model parameters
##
## a list of internal model parameters controlling the model fitting
##
## @param nlam parameter in the augmented Lagrangian function
## @param adaptive if TRUE, use leverage values for adaptive penalization
## @param weights user supplied weights for adaptive penalization
## @param minlam maximum proportion of outliers
## @param maxlam maximum proportion of good observations
## @param delid discarded observation indices for initial estimation
##
## @return a list of model parameters.
r4.modstr <- function(nlam = 100,
                      adaptive = TRUE,
                      weights = NULL,
                      minlam = 0.3,
                      maxlam = 1,
                      delid = NULL)
{
  list(
    nlam = nlam,
    adaptive = adaptive,
    weights = weights,
    minlam = minlam,
    maxlam = maxlam,
    delid = delid
  )
}


## Internal function for specifying computation parameters
##
## a list of internal computational parameters controlling optimization
##
## @param epsilon convergence tolerance
## @param maxit maximum number of iterations
## @param qr.tol tolerance for qr decomposition
## @param tol tolerance
##
## @return a list of computational parameters.
r4.control <- function(epsilon = 1e-3,
                       maxit = 100L,
                       qr.tol = 1e-4,
                       tol = 1e-06)
{
  list(
    epsilon = epsilon,
    maxit = maxit,
    qr.tol = qr.tol,
    tol = tol
  )
}
