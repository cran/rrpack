##' Multivariate reduced-rank linear regression
##'
##' Produce solution paths of reduced-rank estimators and adaptive nuclear norm
##' penalized estimators; compute the degrees of freeom of the RRR estimators
##' and select a solution via certain information criterion.
##'
##' Model parameters can be specified through argument \code{modstr}.  The
##' available include
##' \itemize{
##'    \item{gamma}: A scalar power parameter of the adaptive weights in
##' \code{penalty == "ann"}.
##'
##'    \item{nlambda}: The number of lambda values; no effect if
##' \code{penalty == "count"}.
##'
##'    \item{lambda}: A vector of user-specified rank values if
##' \code{penalty == "count"} or a vector of penalty values if \code{penalty ==
##' "ann"}.
##' }
##'
##' The available elements for argument \code{control} include
##' \itemize{
##'    \item{sv.tol}: singular value tolerence.
##'    \item{qr.tol}: QR decomposition tolerence.
##' }
##'
##' @usage
##' rrr(Y, X, penaltySVD = c("rank", "ann"),
##'     ic.type = c("GIC","AIC","BIC","BICP","GCV"),
##'     df.type = c("exact","naive"), maxrank = min(dim(Y), dim(X)),
##'     modstr = list(), control = list())
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param penaltySVD `rank': rank-constrainted estimation; `ann': adaptive
##'     nuclear norm estimation.
##' @param maxrank an integer of maximum desired rank.
##' @param ic.type the information criterion to be used; currently supporting
##'     `AIC', `BIC', `BICP', `GCV', and `GIC'.
##' @param df.type `exact': the exact degrees of freedoms based on SURE theory;
##'     `naive': the naive degress of freedoms based on counting number of free
##'     parameters
##' @param modstr a list of model parameters controlling the model fitting
##' @param control a list of parameters for controlling the fitting process:
##'     `sv.tol' controls the tolerence of singular values; `qr.tol' controls
##'     the tolerence of QR decomposition for the LS fit
##'
##' @return
##'
##' S3 \code{rrr} object, a list consisting of
##' \item{call}{original function call}
##' \item{Y}{input matrix of response}
##' \item{X}{input matrix of covariate}
##' \item{A}{right singular matri x of the least square fitted matrix}
##' \item{Ad}{a vector of squared singular values of the least square
##'     fitted matrix}
##' \item{coef.ls}{coefficient estimate from LS}
##' \item{Spath}{a matrix, each column containing shrinkage factors of the
##'     singular values of a solution; the first four objects can be
##' used to recover all reduced-rank solutions}
##' \item{df.exact}{the exact degrees of freedom}
##' \item{df.naive}{the naive degrees of freedom}
##' \item{penaltySVD}{the method of low-rank estimation}
##' \item{sse}{a vecotr of sum of squard errors}
##' \item{ic}{a vector of information criterion}
##' \item{coef}{estimated coefficient matrix}
##' \item{U}{estimated left singular matrix such that XU/sqrtn is orthogonal}
##' \item{V}{estimated right singular matrix that is orthogonal}
##' \item{D}{estimated singular value matrix such that C = UDVt}
##' \item{rank}{estimated rank}
##'
##' @references
##'
##' Chen, K., Dong, H. and Chan, K.-S. (2013) Reduced rank regression via
##' adaptive nuclear norm penalization. \emph{Biometrika}, 100, 901--920.
##'
##' @examples
##' library(rrpack)
##' p <- 50; q <- 50; n <- 100; nrank <- 3
##' mydata <- rrr.sim1(n, p, q, nrank, s2n = 1, sigma = NULL,
##'                    rho_X = 0.5, rho_E = 0.3)
##' rfit <- with(mydata, rrr(Y, X, maxrank = 10))
##' summary(rfit)
##' coef(rfit)
##' plot(rfit)
##' @export
rrr <- function(Y,
                X,
                penaltySVD = c("rank", "ann"),
                ic.type = c("GIC", "AIC", "BIC", "BICP", "GCV"),
                df.type = c("exact", "naive"),
                maxrank = min(dim(Y), dim(X)),
                modstr = list(),
                control = list())
{
  ## record function call
  Call <- match.call()

  ## match arguments
  penaltySVD <- match.arg(penaltySVD)
  ic.type <- match.arg(ic.type)
  df.type <- match.arg(df.type)

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  if (n != nrow(X))
    stop("'nrow(X)' has to be equal to 'nrow(Y)'.")

  control <- do.call("rrr.control", control)
  modstr <- do.call("rrr.modstr", modstr)

  ## Obtain the LS estimate
  qrX <- qr(X, tol = control$qr.tol)
  C_ls <- qr.coef(qrX, Y)
  C_ls <- ifelse(is.na(C_ls), 0, C_ls) # FIXME
  rX <- qrX$rank

  nrank <-
    min(q, rX, maxrank)  # nrank is the user specified upper bound
  rmax <-
    min(rX, q)            # rmax is the maximum rank possible
  XC <- qr.fitted(qrX, Y)       # X %*% C_ls
  svdXC <- svd(XC, nu = rmax, nv = rmax)
  A <- svdXC$v[, 1:rmax]
  Ad <- (svdXC$d[1:rmax]) ^ 2

  ## FIXME: is the following necessary?
  Ad <- Ad[Ad > control$sv.tol]
  rmax <- length(Ad)            # modify rmax to be more practical

  ## f.d: weight function based on given thresholding rule
  ## f.d.derv: derivative of f.
  ## lambda: sequence of lambda values

  if (identical(penaltySVD, "ann")) {
    gamma <- modstr$gamma
    nlambda <- modstr$nlambda
    lambda <- modstr$lambda
    Adx <- Ad ^ ((1 + gamma) / 2)
    if (is.null(lambda)) {
      lambda <- exp(seq(log(max(Adx)), log(Adx[min(nrank + 1, rmax)]),
                        length = nlambda))
    }
    f.d <- function(lambda, Ad) {
      Adx <- (Ad) ^ ((1 + gamma) / 2)
      softTH(Adx, lambda) / Adx
    }
    f.d.derv <- function(lambda, Ad) {
      lambda * (gamma + 1) * Ad ^ ((-gamma - 2) / 2)
    }
  } else if (identical(penaltySVD, "rank")) {
    lambda <- modstr$lambda
    if (is.null(lambda))
      lambda <- seq_len(nrank)
    f.d <- function(lambda, Ad)
      rep(1, lambda)
    f.d.derv <- function(lambda, Ad)
      rep(0, lambda)
  }

  nlam <- length(lambda)

  ## Shrinkage factor
  Spath <- matrix(0, nrank, nlam)

  df.exact <- df.naive <- rep(0., nlam)
  for (i in seq_len(nlam)) {
    f <- f.d(lambda[i], Ad[seq_len(nrank)])
    Spath[seq_along(f), i] <- f
    r <- sum(f > control$sv.tol)
    if (r >= 1) {
      if (identical(df.type, "exact")) {
        # compute this only when requested
        f.derv <- f.d.derv(lambda[i], Ad[seq_len(r)])
        term1 <- max(rX, q) * sum(f)

        a <- vector()
        count = 1
        for (k in seq_len(r)) {
          for (s in (r + 1):rmax) {
            a[count] <- (Ad[k] + Ad[s]) * f[k] / (Ad[k] - Ad[s])
            count <- count + 1
          }
        }
        term2 <- sum(a)
        ## h <- length(a)
        ## if (r == maxrank) term2 <- sum(a[-c(h-1,h)])

        if (r == rmax)
          term2 <- 0
        ## FIXME: the following condition can be omitted?
        if (r == rmax & r == min(p, q))
          term2 <- 0

        b <- vector()
        count = 1
        for (k in seq_len(r)) {
          for (s in seq_len(r)) {
            if (s == k) {
              b[count] <- 0
            } else {
              b[count] <- Ad[k] * (f[k] - f[s]) / (Ad[k] - Ad[s])
            }
            count <- count + 1
          }
        }
        term3 <- sum(b)
        term4 <- sum(sqrt(Ad[seq_len(r)]) * f.derv)
        df.exact[i] <- term1 + term2 + term3 + term4
      }
      df.naive[i] <- r * (rX + q - r)
    }
  }

  ## code from rrr.select
  tempFit <- X %*% C_ls
  rankall <- sse <- rep(0., nlam)
  for (i in seq_len(nlam)) {
    dsth <- Spath[, i]
    rank <- sum(dsth != 0)
    rankall[i] <- rank
    if (rank != 0) {
      tempC <- A[, seq_len(rank)] %*%
        (dsth[seq_len(rank)] * t(A[, seq_len(rank)]))
      tempYhat <- tempFit %*% tempC
      sse[i] <- sum((Y - tempYhat) ^ 2)
    } else {
      sse[i] <- sum(Y ^ 2)
    }
  }

  logsse <- log(sse)
  df <- switch(df.type,
               "exact" = df.exact,
               "naive" = df.naive)
  ic <- switch(
    ic.type,
    "GCV" = n * q * sse / (n * q - df) ^ 2,
    "AIC" = n * q * log(sse / n / q) + 2 * df,
    "BIC" = n * q * log(sse / n / q) + log(q * n) * df,
    "BICP" = n * q * log(sse / n / q)  + 2 * log(p * q) * df,
    "GIC" = n * q * log(sse / n / q) +
      log(log(n * q)) * log(p * q) * df
  )

  min.id <- which.min(ic)
  rankest <- rankall[min.id]
  dsth <- Spath[, min.id]
  if (rankest != 0) {
    U <- C_ls %*% A[, seq_len(rankest)] %*%
      diag(1 / svdXC$d[seq_len(rankest)],
           nrow = rankest, ncol = rankest) * sqrt(n)
    D <-
      diag(svdXC$d[seq_len(rankest)] * dsth[seq_len(rankest)],
           nrow = rankest, ncol = rankest) / sqrt(n)
    V <- A[, seq_len(rankest)]
    C <- U %*% D %*% t(V)
  } else {
    C <- matrix(nrow = p, ncol = q, 0)
    U <- matrix(nrow = p, ncol = 1, 0)
    V <- matrix(nrow = q, ncol = 1, 0)
    D <- 0
  }

  rval <- list(
    call = Call,
    Y = Y,
    X = X,
    A = A,
    Ad = Ad,
    coef.ls = C_ls,
    ## singular value shrinkage after thresholding
    Spath = Spath,
    df.exact = df.exact,
    df.naive = df.naive,
    penaltySVD = penaltySVD,
    sse = sse,
    ic = ic,
    coef = C,
    U = U,
    V = V,
    D = D,
    rank = rankest
  )
  class(rval) <- "rrr"
  rval
}



##' Reduced-rank regression with rank selected by cross validation
##'
##' Reduced-rank regression with rank selected by cross validation
##'
##' @usage
##' cv.rrr(Y, X, nfold = 10, maxrank = min(dim(Y), dim(X)),
##'        norder = NULL, coefSVD = FALSE)
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nfold number of folds
##' @param maxrank maximum rank allowed
##' @param norder for constructing the folds
##' @param coefSVD If TRUE, svd of the coefficient is retuned
##'
##' @return a list containing rr estimates from cross validation
##'
##' @references
##' Chen, K., Dong, H. and Chan, K.-S. (2013) Reduced rank regression via
##' adaptive nuclear norm penalization. \emph{Biometrika}, 100, 901--920.
##'
##' @examples
##' library(rrpack)
##' p <- 50; q <- 50; n <- 100; nrank <- 3
##' mydata <- rrr.sim1(n, p, q, nrank, s2n = 1, sigma = NULL,
##'                    rho_X = 0.5, rho_E = 0.3)
##' rfit_cv <- with(mydata, cv.rrr(Y, X, nfold = 10, maxrank = 10))
##' summary(rfit_cv)
##' coef(rfit_cv)
##' @export
cv.rrr <- function (Y,
                    X,
                    nfold = 10,
                    maxrank = min(dim(Y), dim(X)),
                    norder = NULL,
                    coefSVD = FALSE)
{
  ## record function call
  Call <- match.call()

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  ndel <- round(n / nfold)
  if (is.null(norder))
    norder <- sample(seq_len(n), n)

  cr_path <- matrix(ncol = nfold, nrow = maxrank + 1, NA)
  for (f in seq_len(nfold)) {
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    Xf <- X[-iddel,]
    Xfdel <- X[iddel,]
    Yf <- Y[-iddel,]
    Yfdel <- Y[iddel,]
    ini <- rrr.fit(Yf, Xf, nrank = maxrank, coefSVD =  coefSVD)
    C_ls <- ini$coef.ls
    A <- ini$A
    tempFit <- Xfdel %*% C_ls
    tempC <- matrix(nrow = q, ncol = q, 0)
    for (i in seq_len(maxrank)) {
      tempC <- tempC + tcrossprod(A[, i])
      cr_path[i + 1, f] <-
        sum((Yfdel - tempFit %*% tempC) ^ 2)
    }
    cr_path[1, f] <- sum(Yfdel ^ 2)
  }
  index <- order(colSums(cr_path))
  crerr <- rowSums(cr_path[, index]) / length(index) * nfold
  minid <- which.min(crerr)
  rankest <- minid - 1
  ini <- rrr.fit(Y, X, nrank = maxrank)
  C_ls <- ini$coef.ls
  A <- ini$A

  out <- if (identical(rankest, 0)) {
    list(
      call = Call,
      cr.path = cr_path,
      cr.error = crerr,
      norder = norder,
      coef = matrix(0, nrow = p, ncol = q),
      rank = 0,
      coef.ls = C_ls
    )
  }
  else {
    list(
      call = Call,
      cr.path = cr_path,
      cr.error = crerr,
      norder = norder,
      coef.ls = C_ls,
      coef = C_ls %*% tcrossprod(A[, seq_len(rankest)]),
      rank = rankest
    )
  }
  class(out) <- "cv.rrr"
  out
}



##' Fitting reduced-rank regression with a specific rank
##'
##' Given a response matrix and a covariate matrix, this function fits reduced
##' rank regression for a specified rank. It reduces to singular value
##' decomposition if the covariate matrix is the identity matrix.
##'
##' @usage
##' rrr.fit(Y, X, nrank = 1, weight = NULL, coefSVD = FALSE)
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param nrank an integer specifying the desired rank
##' @param weight a square matrix of weight (q by q); The default is the
##'     identity matrix
##' @param coefSVD logical indicating the need for SVD for the coeffient matrix
##'     in the output; used in ssvd estimation
##' @return S3 \code{rrr} object, a list consisting of \item{coef}{coefficient
##'     of rrr} \item{coef.ls}{coefficient of least square} \item{fitted}{fitted
##'     value of rrr} \item{fitted.ls}{fitted value of least square}
##'     \item{A}{right singular matrix} \item{Ad}{a vector of sigular values}
##'     \item{rank}{rank of the fitted rrr}
##' @examples
##' Y <- matrix(rnorm(400), 100, 4)
##' X <- matrix(rnorm(800), 100, 8)
##' rfit <- rrr.fit(Y, X, nrank = 2)
##' coef(rfit)
##' @importFrom MASS ginv
##' @export
rrr.fit <- function(Y,
                    X,
                    nrank = 1,
                    weight = NULL,
                    coefSVD = FALSE)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  stopifnot(n == nrow(X))

  S_yx <- crossprod(Y, X)
  S_xx <- crossprod(X)

  ## FIXME: 0.1 is too arbitrary
  S_xx_inv <- tryCatch(
    ginv(S_xx),
    error = function(e)
      solve(S_xx + 0.1 * diag(p))
  )

  ## FIXME: if weighted, this needs to be weighted too
  C_ls <- tcrossprod(S_xx_inv, S_yx)

  if (!is.null(weight)) {
    stopifnot(nrow(weight) == q && ncol(weight) == q)
    eigenGm <- eigen(weight)
    ## FIXME: ensure eigen success?
    ## sqrtGm <- tcrossprod(eigenGm$vectors * sqrt(eigenGm$values),
    ##                      eigenGm$vectors)
    ## sqrtinvGm <- tcrossprod(eigenGm$vectors / sqrt(eigenGm$values),
    ##                         eigenGm$vectors)
    sqrtGm <- eigenGm$vectors %*% (sqrt(eigenGm$values) *
                                     t(eigenGm$vectors))
    sqrtinvGm <- eigenGm$vectors %*% (1 / sqrt(eigenGm$values) *
                                        t(eigenGm$vectors))

    XC <- X %*% C_ls %*% sqrtGm
    ## FIXME: SVD may not converge
    ## svdXC <- tryCatch(svd(XC,nu=nrank,nv=nrank),error=function(e)2)
    svdXC <- svd(XC, nrank, nrank)
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
    AA <- tcrossprod(A)
    C_rr <- C_ls %*% sqrtGm %*% AA %*% sqrtinvGm
  } else {
    ## unweighted
    XC <- X %*% C_ls
    svdXC <- svd(XC, nrank, nrank)
    A <- svdXC$v[, 1:nrank]
    Ad <- (svdXC$d[1:nrank]) ^ 2
    AA <- tcrossprod(A)
    C_rr <- C_ls %*% AA
  }

  ret <- list(
    call = Call,
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    rank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
    coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  class(ret) <- "rrr.fit"
  ret
}



##' Leverage scores and Cook's distance in reduced-rank regression
##' for model diagnostics
##'
##' Compute leverage scores and Cook's distance for model diagnostics
##' in \code{rrr} estimation.
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param nrank an integer specifying the desired rank
##' @param qr.tol tolerence to be passed to `qr'
##' @return `rrr.leverage' returns a list containing a vector of leverages
##'  and a scalar of the degrees of freedom (sum of leverages).
##' `rrr.cooks' returns a list containing
##'   \item{residuals}{resisuals matrix}
##'   \item{mse}{mean squared error}
##'   \item{leverage}{leverage}
##'   \item{cookD}{Cook's distance}
##'   \item{df}{degrees of freedom}
##' @references
##' Chen, K. Model diagnostics in reduced-rank estimation. \emph{Statistics and
##' Its interface}, 9, 469--484.
##' @importFrom MASS ginv
##' @export
rrr.leverage <- function(Y,
                         X = NULL,
                         nrank = 1,
                         qr.tol = 1e-7)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)

  if (!is.null(X)) {
    p <- ncol(X)

    S_xy <- crossprod(X, Y)
    S_xx <- crossprod(X)
    eigen_xx <- eigen(S_xx)
    eigen_xx$values
    qrX <- qr(X, tol = qr.tol)
    r <- qrX$rank

    Q <- eigen_xx$vectors[, 1:r]
    S <- diag(sqrt(eigen_xx$values[1:r]))
    Sinv <- diag(1 / sqrt(eigen_xx$values[1:r]))
    QSinv <- Q %*% Sinv
    H <- crossprod(QSinv, S_xy)
    XQSinv <- X %*% QSinv
  } else {
    H <- Y
    r <- n
    p <- n
  }

  svd_H <- svd(H)
  U <- as.matrix(svd_H$u[, 1:nrank])
  d <- svd_H$d[1:nrank]
  V <- as.matrix(svd_H$v[, 1:nrank])

  ## decompose df
  lev <- matrix(nrow = n, ncol = q, 0)

  ## compute the common quatities
  HHinv <- array(dim = c(q, q, nrank), 0)
  for (k in 1:nrank) {
    HHinv[, , k] <- ginv(t(H) %*% H - diag(q) * d[k] ^ 2)
  }
  HHinvH <- array(dim = c(q, r, nrank), 0)
  for (k in 1:nrank) {
    HHinvH[, , k] <- HHinv[, , k] %*% t(H)
  }

  VV <- array(dim = c(q, q, nrank), 0)
  for (k in 1:nrank) {
    VV[, , k] <- V[, k] %*% t(V[, k])
  }

  VH <- t(H %*% V)

  for (j in 1:q) {
    PART1 <- sum(V[j, ] ^ 2)
    aa <- HHinvH[, , 1] * V[j, 1] ^ 2
    k  <-  2
    while (k <= nrank) {
      aa <- aa + HHinvH[, , k] * V[j, k] ^ 2
      k <- k + 1
    }
    PART2 <-
      aa + HHinv[, j, ] %*% diag(V[j, ], nrow = nrank, ncol = nrank) %*% VH

    aa <- matrix(nrow = q, ncol = r, 0)
    for (k in 1:nrank) {
      bb <- H * V[j, k]
      bb[, j] <- bb[, j] + H %*% V[, k]
      aa <- aa + V[, k] %*% t(bb %*% HHinv[, j, k])
    }
    PART3 <- aa

    Hj <- -H %*% (PART2 + PART3)
    diag(Hj) <- diag(Hj) + PART1
    ##sum(abs(Hjnew-Hj))

    if (!is.null(X)) {
      lev[, j] <- diag(XQSinv %*% Hj %*% t(XQSinv))
    } else{
      lev[, j] <- diag(Hj)
    }
    ##cat("j = ",j,"\n")
  }

  ## lev is n x q matrix (same as Y); weight in the final solution
  ## df: degrees of freedom
  out <- list(call = Call,
              leverage = lev,
              df = sum(lev))
  class(out) <- "rrr.leverage"
  out
}



##' Cook's distance in reduced-rank regression for model diagnostics
##'
##' Compute Cook's distance for model diagnostics in \code{rrr} estimation.
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nrank model rank
##' @param qr.tol tolerance
##' @return a list containing diagnostics measures
##' @references
##' Chen, K. Model diagnostics in reduced-rank estimation. \emph{Statistics and
##' Its interface}, 9, 469--484.
##' @importFrom MASS ginv
##' @export
rrr.cookD <- function(Y,
                      X = NULL,
                      nrank = 1,
                      qr.tol = 1e-7)
{
  Call <- match.call()

  q <- ncol(Y)
  n <- nrow(Y)

  if (!is.null(X)) {
    p <- ncol(X)

    S_xy <- crossprod(X, Y)
    S_xx <- crossprod(X)
    eigen_xx <- eigen(S_xx)
    eigen_xx$values
    qrX <- qr(X, tol = qr.tol)
    r <- qrX$rank

    Q <- eigen_xx$vectors[, 1:r]
    S <- diag(sqrt(eigen_xx$values[1:r]))
    Sinv <- diag(1 / sqrt(eigen_xx$values[1:r]))
    QSinv <- Q %*% Sinv
    H <- crossprod(QSinv, S_xy)
    XQSinv <- X %*% QSinv
  } else {
    H <- Y
    r <- n
    p <- n
  }

  svd_H <- svd(H, nu = nrank, nv = nrank)
  U <- as.matrix(svd_H$u)
  d <- svd_H$d[1:nrank]
  V <- as.matrix(svd_H$v)

  Vprod <- tcrossprod(V)
  ## Residuals
  if (!is.null(X)) {
    Res <- Y - XQSinv %*% H %*% Vprod
  } else{
    Res <- Y - H %*% Vprod
  }

  mse <- sum(Res ^ 2) / (n - 1) / q

  ## compute the common quantities
  HHinv <- array(dim = c(q, q, nrank), 0)
  for (k in 1:nrank) {
    HHinv[, , k] <- ginv(t(H) %*% H - diag(q) * d[k] ^ 2)
  }
  VV <- array(dim = c(q, q, nrank), 0)
  for (k in 1:nrank) {
    VV[, , k] <- V[, k] %*% t(V[, k])
  }

  ## leverage scores
  lev <- matrix(nrow = n, ncol = q, 0)
  ## sum of squared terms in cook's distance
  Yss <- matrix(nrow = n, ncol = q, 0)
  ## par Hhat vs H[,j]
  Hj <- array(dim = c(r, r, q), 0)
  ## par Yhat vs Y[,j]
  Yj <- array(dim = c(n, n, q), 0)
  ## Ydiv <- array(dim=c(n,n,q,q),0)
  for (j in 1:q) {
    for (i in 1:r) {
      ## Comupte par H vs H[i,j]
      temp <- matrix(nrow = q, ncol = q, 0)
      HZij <- matrix(nrow = q, ncol = q, 0)
      HZij[, j] <- H[i, ]
      HZij[j, ] <- HZij[j, ] + H[i, ]
      for (k in 1:nrank) {
        aa <- HHinv[, , k] %*% HZij %*% VV[, , k]
        temp <- temp + aa
      }
      temp <- temp + t(temp)

      Hj[, i, ] <-  -H %*% temp
      Hj[i, i, ] <- Hj[i, i, ] + Vprod[, j]
    }

    if (!is.null(X)) {
      for (h in 1:q) {
        Yj[, , h] <- XQSinv %*% Hj[, , h] %*% t(XQSinv)
      }
      ##Is this correct? c(1) or c(2)
      ##Yss[,j] <- apply(Ydiv[,,,j],c(2),function(a)sum(a^2))
      Yss[, j] <- apply(Yj, c(2), function(a)
        sum(a ^ 2))
      ##lev[,j] <- diag(Ydiv[,,j,j])
      lev[, j] <- diag(Yj[, , j])
    } else {
      Yss[, j] <- apply(Hj, 2, function(a)
        sum(a ^ 2))
      lev[, j] <- diag(Hj[, , j])
    }
  }

  cook_entry <- Yss / (1 - lev) ^ 2 * Res ^ 2 / r / mse

  out <- list(
    call = Call,
    residuals = Res,
    mse = mse,
    leverage = lev,
    cookD = cook_entry,
    df = sum(lev)
  )
  class(out) <- "rrr.cookD"
  out
}


### internal functions =========================================================
rrr.control <- function(sv.tol = 1e-7, qr.tol = 1e-7)
{
  ## FIXME: add some simple checks?
  list(sv.tol = sv.tol,               # singular value tolerence
       qr.tol = qr.tol)               # QR decomposition tolerence
}

rrr.modstr <- function(gamma = 2,
                       nlambda = 100,
                       lambda = NULL)
{
  list(gamma = gamma,
       nlambda = nlambda,
       lambda = lambda)
}
