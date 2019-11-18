##' Sparse orthogonal factor regression
##'
##' Compute solution paths of sparse orthogonal factor regression
##'
##' The model parameters can be specified through argument \code{modstr}.
##' The available elements include
##' \itemize{
##'     \item{mu}: parameter in the augmented Lagrangian function.
##'     \item{mugamma}: increament of mu along iterations to speed up
##'         computation.
##'     \item{WA}: weight matrix for A.
##'     \item{WB}: weight matrix for B.
##'     \item{Wd}: weight matrix for d.
##'     \item{wgamma}: power parameter in constructing adaptive weights.
##' }
##' The model fitting can be controled through argument \code{control}.
##' The avilable elements include
##' \itemize{
##'    \item{nlam}: number of lambda triplets to be used.
##'    \item{lam.min.factor}: set the smallest lambda triplets as a fraction of the
##'     estimation lambda.max triplets.
##'    \item{lam.max.factor}: set the largest lambda triplets as a multiple of the
##'     estimation lambda.max triplets.
##'    \item{lam.AB.factor}: set the relative penalty level between A/B and D.
##'    \item{penA,penB,penD}: if TRUE, penalty is applied.
##'    \item{lamA}: sequence of tuning parameters for A.
##'    \item{lamB}: sequence of tuning parameters for B.
##'    \item{lamD}: sequence of tuning parameters for d.
##'    \item{methodA}: penalty for penalizing A.
##'    \item{methodB}: penalty for penalizing B.
##'    \item{epsilon}: convergence tolerance.
##'    \item{maxit}: maximum number of iterations.
##'    \item{innerEpsilon}: convergence tolerance for inner subroutines.
##'    \item{innerMaxit}: maximum number of iterations for inner subroutines.
##'    \item{sv.tol}: tolerance for singular values.
##' }
##'
##' @name sofar
##'
##' @usage
##' sofar(Y, X, nrank = 1, su = NULL, sv = NULL,
##'       ic.type = c("GIC", "BIC", "AIC", "GCV"),
##'       modstr = list(), control = list(), screening = FALSE)
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nrank an integer specifying the desired rank/number of factors
##' @param su a scaling vector for U such that \eqn{U^TU = diag(s_u)}.
##' @param sv a scaling vector for V such that \eqn{V^TV = diag(s_v)}.
##' @param ic.type select tuning method; the default is GIC
##' @param modstr a list of internal model parameters controlling the model
##'     fitting
##' @param control a list of internal computation parameters controlling
##'     optimization
##' @param screening If TRUE, marginal screening via lasso is performed before
##'     sofar fitting.
##'
##' @return
##' A \code{sofar} object containing
##'   \item{call}{original function call}
##'   \item{Y}{input response matrix}
##'   \item{X}{input predictor matrix}
##'   \item{Upath}{solution path of U}
##'   \item{Dpath}{solution path of D}
##'   \item{Vpath}{solution path of D}
##'   \item{Rpath}{path of estimated rank}
##'   \item{icpath}{path of information criteria}
##'   \item{lam.id}{ids of selected lambda for GIC, BIC, AIC and GCV}
##'   \item{p.index}{ids of predictors which passed screening}
##'   \item{q.index}{ids of responses which passed screening}
##'   \item{lamA}{tuning sequence for A}
##'   \item{lamB}{tuning sequence for B}
##'   \item{lamD}{tuning sequence for D}
##'   \item{U}{estimated left singular matrix that is orthogonal (factor weights)}
##'   \item{V}{estimated right singular matrix that is orthogonal (factor loadings)}
##'   \item{D}{estimated singular values}
##'   \item{rank}{estimated rank}
##'
##' @references
##'
##' Uematsu, Y., Fan, Y., Chen, K., Lv, J., & Lin, W. (2019). SOFAR: large-scale
##' association network learning. \emph{IEEE Transactions on Information
##' Theory}, 65(8), 4924--4939.
##'
##' @examples
##' \dontrun{
##' library(rrpack)
##' ## Simulate data from a sparse factor regression model
##' p <- 100; q <- 50; n <- 100; nrank <- 3
##' mydata <- rrr.sim1(n, p, q, nrank, s2n = 1,
##'                    sigma = NULL, rho_X = 0.5, rho_E = 0.3)
##' Y <- mydata$Y
##' X <- mydata$X
##'
##' fit1 <- sofar(Y, X, ic.type = "GIC", nrank = nrank + 2,
##'               control = list(methodA = "adlasso", methodB = "adlasso"))
##' summary(fit1)
##' plot(fit1)
##'
##' fit1$U
##' crossprod(fit1$U) #check orthogonality
##' fit1$V
##' crossprod(fit1$V) #check orthogonality
##' }
##'
##' @importFrom Rcpp evalCpp
##' @useDynLib rrpack
##' @export
sofar <- function(Y,
                  X,
                  nrank = 1,
                  su = NULL, sv = NULL,
                  ic.type = c("GIC", "BIC", "AIC", "GCV"),
                  modstr = list(),
                  control = list(),
                  screening = FALSE)
{
  Call <- match.call()

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if(is.null(su)) su <- rep(1,nrank)
  #   Su_sqrt <- diag(nrow=nrank, ncol=nrank)
  #   #Su_sqrt_inv <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Su_sqrt <- diag(sqrt(su),nrow=nrank, ncol=nrank)
  #   #Su_sqrt_inv <- diag(1/sqrt(su),nrow=nrank, ncol=nrank)
  # }
  if(is.null(sv))  sv <- rep(1,nrank)
  #   Sv_sqrt <- diag(nrow=nrank, ncol=nrank)
  #   #Sv_sqrt_inv <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Sv_sqrt <- diag(sqrt(sv),nrow=nrank, ncol=nrank)
  #   #Sv_sqrt_inv <- diag(1/sqrt(sv),nrow=nrank, ncol=nrank)
  # }

  ## model parameters
  modstr <- do.call("sofar.modstr", modstr)
  mu <- modstr$mu
  mugamma <- modstr$mugamma
  WA <- modstr$WA
  WB <- modstr$WB
  Wd <- modstr$Wd
  wgamma <- modstr$wgamma

  ## control
  control <- do.call("sofar.control", control)
  nlam <- control$nlam
  lam.min.factor <- control$lam.min.factor
  lam.max.factor <- control$lam.max.factor
  lam.AB.factor <- control$lam.AB.factor
  penA <- control$penA
  penB <- control$penB
  penD <- control$penD
  lamA <- control$lamA
  lamB <- control$lamB
  lamD <- control$lamD
  methodA <- control$methodA
  methodB <- control$methodB
  epsilon <- control$epsilon
  maxit <- control$maxit
  innerEpsilon <- control$innerEpsilon
  innerMaxit <- control$innerMaxit
  tol <- control$sv.tol

  ## Initial screening from lasso
  if(screening == TRUE){
    ini.screen <- sofar.init(Y, X, nrank)
    p.index <- ini.screen$p.index
    q.index <- ini.screen$q.index
  } else{
    p.index <- 1:p
    q.index <- 1:q
  }
  ## U <- init$coef.svd$U
  ## V <- init$coef.svd$V
  ## D <- init$coef.svd$D

  ## Reduce the problem and compute lambda sequence
  Y.reduced <- Y[, q.index]
  X.reduced <- X[, p.index]
  q.reduced <- length(q.index)
  p.reduced <- length(p.index)

  ## Construct weights and initial value.
  ini <- NULL
  if (is.null(WA) || is.null(WB) || is.null(Wd)) {
    if (is.null(ini))
      ini <-
        rrr.fit(Y.reduced, X.reduced, nrank, coefSVD = TRUE)

    WA <- switch(
      methodA,
      lasso = matrix(1., p.reduced, nrank),
      adlasso = (
          ini$coefSVD$u %*%
          diag(ini$coefSVD$d / sqrt(sv),
               nrow = length(ini$coefSVD$d)) + 1e-8
      ) ^ (-wgamma),
      glasso = rep(1., p.reduced),
      adglasso = sqrt(rowSums((
        ini$coefSVD$u %*%
          diag(ini$coefSVD$d/sqrt(sv),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )

    WB <- switch(
      methodB,
      lasso = matrix(1., q.reduced, nrank),
      adlasso = (ini$coefSVD$v %*%
                 diag(ini$coefSVD$d/sqrt(su),
                      nrow = length(ini$coefSVD$d)) + 1e-8
      ) ^ (-wgamma),
      glasso = rep(1., q.reduced),
      adglasso = sqrt(rowSums((
        ini$coefSVD$v %*%
          diag(ini$coefSVD$d/sqrt(su),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )

    if (substr(methodA, 1, 2) == "ad" ||
        substr(methodB, 1, 2) == "ad") {
      if (is.null(ini))
        ini <-
          rrr.fit(Y.reduced, X.reduced, nrank, coefSVD = TRUE)
      Wd <- (ini$coefSVD$d / sqrt(su) / sqrt(sv) + 1e-8) ^ (- wgamma)
    } else {
      Wd <- rep(1., nrank)
    }

    modstr$WA <- WA
    modstr$WB <- WB
    modstr$Wd <- Wd
  }

  ## if (is.null(init$U) || is.null(init$V) || is.null(init$D)) {
  U <- t(sqrt(su)*t(ini$coefSVD$u))
  V <- t(sqrt(sv)*t(ini$coefSVD$v))
  D <- ini$coefSVD$d/sqrt(su)/sqrt(sv)
  init <-
    list(
      U = U,
      V = V,
      D = D
    )
  ## } else {
  ##  U <- init$U
  ##  V <- init$V
  ##  D <- init$D
  ## }


  ## Compute lambda sequence
  if (is.null(lamA) || is.null(lamB) || is.null(lamD)) {
    lammax <- sofar.lammax(
      Y.reduced,
      X.reduced,
      U = U,
      V = V,
      D = diag(D, nrow = nrank),
      nrank = nrank,
      su = su,
      sv = sv,
      mu = mu,
      mugamma = mugamma,
      conv = 1e-4,
      maxit = maxit,
      methodA = methodA,
      methodB = methodB,
      WA = WA,
      WB = WB,
      Wd = Wd,
      wgamma = 2
    )
    lamA.max <- lammax$lamA.max
    lamB.max <- lammax$lamB.max
    lamD.max <- lammax$lamD.max
    if (penA) {
      lamA <- exp(seq(
        log(lamA.max * lam.min.factor),
        log(lamA.max * lam.max.factor),
        length = nlam
      ))
    } else {
      lamA <- rep(0, nlam)
    }
    if (penB) {
      lamB <- exp(seq(
        log(lamB.max * lam.min.factor * lam.AB.factor),
        log(lamB.max * lam.max.factor * lam.AB.factor),
        length = nlam
      ))
    } else {
      lamB <- rep(0, nlam)
    }
    if (penD) {
      lamD <- exp(seq(
        log(lamD.max * lam.min.factor),
        log(lamD.max * lam.max.factor),
        length = nlam
      ))
    } else {
      lamD <- rep(0, nlam)
    }

  } else {
    nlamA <- length(lamA)
    nlamB <- length(lamB)
    nlamD <- length(lamD)
    stopifnot(nlamA == nlamB & nlamB == nlamD)
  }

  XX.reduced <- crossprod(X.reduced)
  eigenXX.reduced <- eigen(XX.reduced)$values[1]
  XY.reduced <- crossprod(X.reduced, Y.reduced)
  Yvec.reduced <- as.vector(Y.reduced)

  fit <- sofar.path.reduced(
    Y.reduced,
    X.reduced,
    XX = XX.reduced,
    XY = XY.reduced,
    Yvec = Yvec.reduced,
    eigenXX = eigenXX.reduced,
    ic = ic.type,
    nrank = nrank,
    su = su,
    sv = sv,
    lamA = lamA,
    lamB = lamB,
    lamD = lamD,
    modstr = modstr,
    init = init,
    control = control
  )

  ## Results
  rank <- fit$rank
  if (rank > 0) {
    U <- matrix(nrow = p, ncol = rank, 0)
    U[p.index,] <- fit$U
    V <- matrix(nrow = q, ncol = rank, 0)
    V[q.index,] <- fit$V
    D <- fit$D
  } else{
    U <- matrix(nrow = p, ncol = 1, 0)
    V <- matrix(nrow = q, ncol = 1, 0)
    D <- fit$D
  }

  ICid <- switch(
    ic.type,
    GIC = fit$lam.id[1],
    BIC = fit$lam.id[2],
    AIC = fit$lam.id[3],
    GCV = fit$lam.id[4]
  )
  if(ICid==1) warning("Selected smallest lambda on the sequence. Decrease lambda.min.factor",call. = FALSE)
  if(ICid==control$nlam) warning("Selected largest lambda on the sequence. Increase lambda.max.factor",call. = FALSE)

    out <- list(
    call = Call,
    Y = Y,
    X = X,
    Upath = fit$Upath,
    Vpath = fit$Vpath,
    Dpath = fit$Dpath,
    Rpath = fit$Rpath,
    icpath = fit$ICpath,
    lam.id = fit$lam.id,
    p.index = p.index,
    q.index = q.index,
    lamA = lamA,
    lamB = lamB,
    lamD = lamD,
    U = U,
    V = V,
    D = D,
    rank = rank
  )
  class(out) <- "sofar"
  out
}


##' Sparse orthognal factor regression tuned by cross validation
##'
##' Sparse orthognal factor regression tuned by cross validation
##'
##' The model parameters can be specified through argument \code{modstr}.
##' The available elements include
##' \itemize{
##'     \item{mu}: parameter in the augmented Lagrangian function.
##'     \item{mugamma}: increament of mu along iterations to speed up
##'         computation.
##'     \item{WA}: weight matrix for A.
##'     \item{WB}: weight matrix for B.
##'     \item{Wd}: weight matrix for d.
##'     \item{wgamma}: power parameter in constructing adaptive weights.
##' }
##' The model fitting can be controled through argument \code{control}.
##' The avilable elements include
##' \itemize{
##'    \item{nlam}: number of lambda triplets to be used.
##'    \item{lam.min.factor}: set the smallest lambda triplets as a fraction of
##'        the estimation lambda.max triplets.
##'    \item{lam.max.factor}: set the largest lambda triplets as a multiple of
##'        the estimation lambda.max triplets.
##'    \item{lam.AB.factor}: set the relative penalty level between A/B and D.
##'    \item{penA,penB,penD}: if TRUE, penalty is applied.
##'    \item{lamA}: sequence of tuning parameters for A.
##'    \item{lamB}: sequence of tuning parameters for B.
##'    \item{lamD}: sequence of tuning parameters for d.
##'    \item{methodA}: penalty for penalizing A.
##'    \item{methodB}: penalty for penalizing B.
##'    \item{epsilon}: convergence tolerance.
##'    \item{maxit}: maximum number of iterations.
##'    \item{innerEpsilon}: convergence tolerance for inner subroutines.
##'    \item{innerMaxit}: maximum number of iterations for inner subroutines.
##'    \item{sv.tol}: tolerance for singular values.
##' }
##'
##' @usage
##' cv.sofar(Y, X, nrank = 1, su = NULL, sv = NULL, nfold = 5, norder = NULL, modstr = list(),
##'          control = list(), screening = FALSE)
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param nrank an integer specifying the desired rank/number of factors
##' @param su a scaling vector for U such that \eqn{U^{T}U = diag(s_{u})}
##' @param sv a scaling vector for V such that \eqn{V^{T}V = diag(s_{v})}
##' @param nfold number of fold; used for cv.sofar
##' @param norder observation orders to constrct data folds; used for cv.sofar
##' @param modstr a list of internal model parameters controlling the model
##'     fitting
##' @param control a list of internal computation parameters controlling
##'     optimization
##' @param screening If TRUE, marginal screening via lasso is performed before
##'     sofar fitting.
##'
##' @export
cv.sofar <- function(Y,
                     X,
                     nrank = 1,
                     su = NULL,
                     sv = NULL,
                     nfold = 5,
                     norder = NULL,
                     modstr = list(),
                     control = list(),
                     screening = FALSE)
{
  Call <- match.call()

  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  if(is.null(su)) su <- rep(1,nrank)
  #   Su_sqrt <- diag(nrow=nrank, ncol=nrank)
  #   #Su_sqrt_inv <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Su_sqrt <- diag(sqrt(su),nrow=nrank, ncol=nrank)
  #   #Su_sqrt_inv <- diag(1/sqrt(su),nrow=nrank, ncol=nrank)
  # }
  if(is.null(sv)) sv <- rep(1,nrank)
  #   Sv_sqrt <- diag(nrow=nrank, ncol=nrank)
  #   #Sv_sqrt_inv <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Sv_sqrt <- diag(sqrt(sv),nrow=nrank, ncol=nrank)
  #   #Sv_sqrt_inv <- diag(1/sqrt(sv),nrow=nrank, ncol=nrank)
  # }


  ## model parameters
  modstr <- do.call("sofar.modstr", modstr)
  mu <- modstr$mu
  mugamma <- modstr$mugamma
  WA <- modstr$WA
  WB <- modstr$WB
  Wd <- modstr$Wd
  wgamma <- modstr$wgamma

  ## control
  control <- do.call("sofar.control", control)
  nlam <- control$nlam
  lam.min.factor <- control$lam.min.factor
  lam.max.factor <- control$lam.max.factor
  lam.AB.factor <- control$lam.AB.factor
  penA <- control$penA
  penB <- control$penB
  penD <- control$penD
  lamA <- control$lamA
  lamB <- control$lamB
  lamD <- control$lamD
  methodA <- control$methodA
  methodB <- control$methodB
  epsilon <- control$epsilon
  maxit <- control$maxit
  innerEpsilon <- control$innerEpsilon
  innerMaxit <- control$innerMaxit
  tol <- control$sv.tol

  ## initial screening from lasso
  if(screening == TRUE){
    ini.screen <- sofar.init(Y, X, nrank)
    p.index <- ini.screen$p.index
    q.index <- ini.screen$q.index
  } else{
    p.index <- 1:p
    q.index <- 1:q
  }

  Y.reduced <- Y[, q.index]
  X.reduced <- X[, p.index]
  q.reduced <- length(q.index)
  p.reduced <- length(p.index)

  ## Construct weights and initial value.
  ini <- NULL
  if (is.null(WA) || is.null(WB) || is.null(Wd)) {
    if (is.null(ini))
      ini <-
        rrr.fit(Y.reduced, X.reduced, nrank, coefSVD = TRUE)

    WA <- switch(
      methodA,
      lasso = matrix(1., p.reduced, nrank),
      adlasso = (ini$coefSVD$u %*%
                   diag(ini$coefSVD$d/sqrt(sv),
                        nrow = length(ini$coefSVD$d))
      ) ^ (-wgamma),
      glasso = rep(1., p.reduced),
      adglasso = sqrt(rowSums((
        ini$coefSVD$u %*%
          diag(ini$coefSVD$d/sqrt(sv),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )

    WB <- switch(
      methodB,
      lasso = matrix(1., q.reduced, nrank),
      adlasso = (ini$coefSVD$v %*%
                   diag(ini$coefSVD$d/sqrt(su),
                        nrow = length(ini$coefSVD$d))
      ) ^ (-wgamma),
      glasso = rep(1., q.reduced),
      adglasso = sqrt(rowSums((
        ini$coefSVD$v %*%
          diag(ini$coefSVD$d/sqrt(su),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )

    if (substr(methodA, 1, 2) == "ad" ||
        substr(methodB, 1, 2) == "ad") {
      if (is.null(ini))
        ini <-
          rrr.fit(Y.reduced, X.reduced, nrank, coefSVD = TRUE)
      Wd <- (ini$coefSVD$d/sqrt(su)/sqrt(sv)) ^ (-wgamma)
    } else {
      Wd <- rep(1., nrank)
    }

    modstr$WA <- WA
    modstr$WB <- WB
    modstr$Wd <- Wd
  }

  ## if (is.null(init$U) || is.null(init$V) || is.null(init$D)) {
  U <- t(sqrt(su)*t(ini$coefSVD$u))
  V <- t(sqrt(sv)*t(ini$coefSVD$v))
  D <- ini$coefSVD$d/sqrt(su)/sqrt(sv)
  init <-
    list(
      U = U,
      V = V,
      D = D
    )
  ## } else {
  ##  U <- init$U
  ##  V <- init$V
  ##  D <- init$D
  ## }

  ## Compute lambda sequence
  if (is.null(lamA) || is.null(lamB) || is.null(lamD)) {
    lammax <- sofar.lammax(
      Y.reduced,
      X.reduced,
      U = U,
      V = V,
      D = diag(D, nrow = nrank),
      nrank = nrank,
      su = su,
      sv = sv,
      mu = mu,
      mugamma = mugamma,
      conv = 1e-4,
      maxit = maxit,
      methodA = methodA,
      methodB = methodB,
      WA = WA,
      WB = WB,
      Wd = Wd,
      wgamma = 2
    )
    lamA.max <- lammax$lamA.max
    lamB.max <- lammax$lamB.max
    lamD.max <- lammax$lamD.max
    if (penA) {
      lamA <- exp(seq(
        log(lamA.max * lam.min.factor),
        log(lamA.max * lam.max.factor),
        length = nlam
      ))
    } else {
      lamA <- rep(0, nlam)
    }
    if (penB) {
      lamB <- exp(seq(
        log(lamB.max * lam.min.factor * lam.AB.factor),
        log(lamB.max * lam.max.factor * lam.AB.factor),
        length = nlam
      ))
    } else {
      lamB <- rep(0, nlam)
    }
    if (penD) {
      lamD <- exp(seq(
        log(lamD.max * lam.min.factor),
        log(lamD.max * lam.max.factor),
        length = nlam
      ))
    } else {
      lamD <- rep(0, nlam)
    }

  } else {
    nlamA <- length(lamA)
    nlamB <- length(lamB)
    nlamD <- length(lamD)
    stopifnot(nlamA == nlamB & nlamB == nlamD)
  }

  XX.reduced <- crossprod(X.reduced)
  eigenXX.reduced <- eigen(XX.reduced)$values[1]
  XY.reduced <- crossprod(X.reduced, Y.reduced)
  Yvec.reduced <- as.vector(Y.reduced)

  fit <- sofar.cv(
    Y.reduced,
    X.reduced,
    XX = XX.reduced,
    XY = XY.reduced,
    Yvec = Yvec.reduced,
    eigenXX = eigenXX.reduced,
    nrank = nrank,
    nfold = nfold,
    su = su,
    sv = sv,
    norder = norder,
    lamA = lamA,
    lamB = lamB,
    lamD = lamD,
    modstr = modstr,
    init = init,
    control = control
  )

  ## Results
  rank <- fit$rank
  if (rank > 0) {
    U <- matrix(nrow = p, ncol = rank, 0)
    U[p.index, ] <- fit$U
    V <- matrix(nrow = q, ncol = rank, 0)
    V[q.index, ] <- fit$V
    D <- fit$D
  } else {
    U <- matrix(nrow = p, ncol = 1, 0)
    V <- matrix(nrow = q, ncol = 1, 0)
    D <- fit$D
  }

  if(fit$lam.id==1) warning("Selected smallest lambda on the sequence. Decrease lambda.min.factor",call. = FALSE)
  if(fit$lam.id==control$nlam) warning("Selected largest lambda on the sequence. Increase lambda.max.factor", call. = FALSE)

  out <- list(
    call = Call,
    Y = Y,
    X = X,
    crpath = fit$cr.path,
    crerror = fit$cr.error,
    norder = fit$norder,
    Upath = fit$Upath,
    Vpath = fit$Vpath,
    Dpath = fit$Dpath,
    Rpath = fit$Rpath,
    icpath = fit$ICpath,
    lam.id = fit$lam.id,
    #lammax = lammax,
    p.index = p.index,
    q.index = q.index,
    lamA = lamA,
    lamB = lamB,
    lamD = lamD,
    U = U,
    V = V,
    D = D,
    rank = rank
  )
  class(out) <- "cv.sofar"
  out
}



### internal functions =========================================================

## Internal function for specifying model parameters
##
## a list of internal model parameters controlling the model fitting
##
## @param mu parameter in the augmented Lagrangian function
## @param mugamma increament of mu along iterations to speed up computation
## @param WA weight matrix for A
## @param WB weight matrix for B
## @param Wd weight matrix for d
## @param wgamma power parameter in constructing adaptive weights
##
## @return a list of model parameters.
sofar.modstr <- function(mu = 1,
                         mugamma = 1.1,
                         WA = NULL,
                         WB = NULL,
                         Wd = NULL,
                         wgamma = 2) {
  list(
    mu = mu,
    mugamma = mugamma,
    WA = WA,
    WB = WB,
    Wd = Wd,
    wgamma = wgamma
  )
}


## Internal function for specifying computation parameters
##
## a list of internal computational parameters controlling optimization
## @param nlam number of lambda triplets to be used
## @param lam.min.factor set the smallest lambda triplets as a fraction of the
##     estimation lambda.max triplets
## @param lam.max.factor set the largest lambda triplets as a multiple of the
##     estimation lambda.max triplets
## @param lam.AB.factor set the relative penalty level between A/B and D
## @param penA,penB,penD if TRUE, penalty is applied
## @param lamA sequence of tuning parameters for A
## @param lamB sequence of tuning parameters for B
## @param lamD sequence of tuning parameters for d
## @param methodA penalty for penalizing A
## @param methodB penalty for penalizing B
## @param epsilon convergence tolerance
## @param maxit maximum number of iterations
## @param innerEpsilon convergence tolerance for inner subroutines
## @param innerMaxit maximum number of iterations for inner subroutines
## @param sv.tol tolerance for singular values
##
## @return a list of computational parameters.
sofar.control <- function(nlam = 50,
                          lam.min.factor = 1e-4,
                          lam.max.factor = 10,
                          lam.AB.factor = 1,
                          penA = TRUE,
                          penB = TRUE,
                          penD = TRUE,
                          lamA = NULL,
                          lamB = NULL,
                          lamD = NULL,
                          methodA = "adlasso",
                          methodB = "adlasso",
                          epsilon = 1e-3,
                          maxit = 200L,
                          innerEpsilon = 1e-3,
                          innerMaxit = 50L,
                          sv.tol = 1e-02)
{
  list(
    nlam = nlam,
    lam.min.factor = lam.min.factor,
    lam.max.factor = lam.max.factor,
    lam.AB.factor = lam.AB.factor,
    penA = penA,
    penB = penB,
    penD = penD,
    lamA = lamA,
    lamB = lamB,
    lamD = lamD,
    methodA = methodA,
    methodB = methodB,
    epsilon = epsilon,
    maxit = maxit,
    innerEpsilon = innerEpsilon,
    innerMaxit = innerMaxit,
    sv.tol = sv.tol
  )
}


### Preparation
### Generate initial values.
##' @importFrom stats coef
##' @importFrom glmnet cv.glmnet
sofar.init <- function(Y, X, nrank = 1)
{
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  ## 1. Use glmnet to get initial values.
  ## 2. Decompose to get initial values
  coef.lasso <- matrix(0, nrow = p, ncol = q)
  for (i in seq_len(q)) {
    fiti <- cv.glmnet(X, Y[, i], intercept = FALSE, standardize = FALSE)
    coef.lasso[, i] <- as.vector(coef(fiti))[-1]
  }
  p.index <- which(apply(coef.lasso, 1, lnorm) != 0)
  q.index <- which(apply(coef.lasso, 2, lnorm) != 0)

  ## Y.reduced <- Y[,q.index]
  ## X.reduced <- X[,p.index]
  ## coef.lasso.reduced <- coef.lasso[p.index,q.index]
  ## if(ini.type == "xc"){
  ##     fit.rrr.reduced <- rrr.fit(Y=Y.reduced,
  ##                                X=X.reduced,
  ##                                nrank = nrank)
  ##     fit.rrr.reduced <- rrr.fit(Y=X.reduced%*%coef.lasso.reduced,
  ##                                X=X.reduced,
  ##                                nrank = nrank)
  ##     coef.rrr.reduced <- fit.rrr.reduced$coef
  ##     coef.svd.reduced <- svd(coef.rrr.reduced,
  ##                             nu = nrank,
  ##                             nv = nrank)
  ## }
  ## if(ini.type == "c"){
  ##     coef.svd.reduced <- svd(coef.lasso.reduced,
  ##                             nu = nrank,
  ##                             nv = nrank)
  ##     coef.rrr.reduced <- coef.svd.reduced$u %*%
  ##         diag(coef.svd.reduced$d[1:nrank],nrow=nrank) %*%
  ##         t(coef.svd.reduced$v)
  ## }

  ## coef.svd$u
  ## coef.svd$d[1:nrank]
  ## coef.svd$v

  list(p.index = p.index,
       q.index = q.index)
}


### Solve the Procrustes problem in SoFAR
sofar.procrustes <- function(XY,
                             XX,
                             D,
                             rho2,
                             U = NULL,
                             identity = FALSE,
                             control = list(maxit = 50L, epsilon = 1e-3)) {
  control <- do.call("sofar.control", control)
  p <- ncol(XX)
  converged <- FALSE

  if (is.null(U) | identity == TRUE) {
    XYsvd <- svd(XY)
    U <- tcrossprod(XYsvd$u, XYsvd$v)
    ## U <- XYsvd$u%*%t(XYsvd$v)
    diff <- NULL
    niter = NULL
    converged <- TRUE
  }

  if (identity == FALSE) {
    epsilon <- control$epsilon
    maxit <- control$maxit
    niter <- 1
    ## rho2 <- eigen(XX)$values[1]
    diff <- vector()
    diff[1] <- 10 * epsilon

    ## D = diag(D,nrow=length(D))
    ## ZZ <- (rho2*diag(nrow=p)-XX)
    ZZ <- -XX
    diag(ZZ) <- rho2 + diag(ZZ)
    XYD <- XY %*% D
    Dsq <- D ^ 2
    while (niter <= maxit & diff[niter] > epsilon) {
      U0 <- U
      ## Z <- ((rho2*diag(nrow=p)-XX)%*%U0%*%D + XY)%*%D
      Z <- ZZ %*% U0 %*% Dsq + XYD
      Zsvd <- svd(Z, nu = nrow(D), nv = nrow(D))
      #U <- Zsvd$u%*%Zsvd$v
      U <- tcrossprod(Zsvd$u, Zsvd$v)
      niter <- niter + 1
      diff[niter] <- sqrt(sum((U - U0) ^ 2))
      ## diff[niter] <- sqrt(sum((Y-X%*%U%*%D)^2))
    }

  }

  list(
    U = U,
    diff = diff,
    niter = niter,
    converged = converged
  )

}


###Fit SoFAR with fixed tuning parameters
##' @importFrom lassoshooting lassoshooting
sofar.fit <- function(Y,
                      X,
                      XX = NULL,
                      XY = NULL,
                      Yvec = NULL,
                      eigenXX = NULL,
                      nrank = 1,
                      su = NULL,
                      sv = NULL,
                      modstr = list(),
                      init = list(U = NULL, V = NULL, D = NULL),
                      control = list()) {
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)

  if(is.null(su)) su <- rep(1,nrank)
    #Su_sqrt <- diag(nrow=nrank, ncol=nrank)
    #Su_sqrt_inv <- diag(nrow=nrank, ncol=nrank)
  # } else{
  #   Su_sqrt <- diag(sqrt(su),nrow=nrank, ncol=nrank)
  #   #Su_sqrt_inv <- diag(1/sqrt(su),nrow=nrank, ncol=nrank)
  # }
  if(is.null(sv)) sv <- rep(1,nrank)
    #Sv_sqrt <- diag(nrow=nrank, ncol=nrank)
    #Sv_sqrt_inv <- diag(nrow=nrank, ncol=nrank)
  # } else{
  #   Sv_sqrt <- diag(sqrt(sv),nrow=nrank, ncol=nrank)
  #   #Sv_sqrt_inv <- diag(1/sqrt(sv),nrow=nrank, ncol=nrank)
  # }


  modstr <- do.call("sofar.modstr", modstr)
  ##init <- do.call("sofar.0", init)
  control <- do.call("sofar.control", control)

  ## FIXME: xrank may become an input is the same X is used in multiple calls
  xrank <- sum(svd(X)$d > control$sv.tol)

  ## init
  U <- init$U
  V <- init$V
  D <- init$D
  ## modstr
  #nrank <- modstr$nrank;
  mu <- modstr$mu
  mugamma <- modstr$mugamma
  WA <- modstr$WA
  WB <- modstr$WB
  Wd <- modstr$Wd
  wgamma <- modstr$wgamma
  ## control
  lamA <-
    control$lamA
  lamB <- control$lamB
  lamD <- control$lamD

  methodA <- control$methodA
  methodB <- control$methodB

  epsilon <- control$epsilon
  maxit <- control$maxit

  innerEpsilon <-
    control$innerEpsilon
  innerMaxit <- control$innerMaxit

  tol <- control$sv.tol

  if (is.null(U) || is.null(V) || is.null(D)) {
    ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
    U <- t(sqrt(su)*t(ini$coefSVD$u))
    V <- t(sqrt(sv)*t(ini$coefSVD$v))
    D <- ini$coefSVD$d/sqrt(su)/sqrt(sv)
  } else {
    ini <- NULL
  }

  ## d <- diag(D,nrow=length(D),ncol=length(D))
  ## FIXME:
  ## diag(D,nrow=length(D),ncol=length(D)) need to have nrow specified for
  ## security diag(0.9) returns a matrix of 0 by 0 instead of a 1 by 1 matrix
  ## of 0.9
  A <- t(t(U) * D)  # U%*%D
  B <- t(t(V) * D)  # V%*%D
  GA <- matrix(0., p, nrank)
  GB <- matrix(0., q, nrank)
  nzid <- 1:nrank

  ## Adaptive weights
  if (is.null(WA) || is.null(WB) || is.null(Wd)) {
    if (is.null(ini))
      ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
    WA <- switch(
      methodA,
      lasso    = matrix(1., p, nrank),
      adlasso  = (ini$coefSVD$u %*%
                    diag(ini$coefSVD$d/sqrt(sv),
                         nrow = length(ini$coefSVD$d))
                  ) ^ (-wgamma),
      glasso   = rep(1., p),
      adglasso =
        sqrt(rowSums((
          ini$coefSVD$u %*%
            diag(ini$coefSVD$d/sqrt(sv),
                 nrow = length(ini$coefSVD$d))
        ) ^ 2)) ^ (-wgamma)
    )
    WB <- switch(
      methodB,
      lasso    = matrix(1., q, nrank),
      adlasso  = (ini$coefSVD$v %*%
                    diag(ini$coefSVD$d/sqrt(su),
                         nrow = length(ini$coefSVD$d))
                  ) ^ (-wgamma),
      glasso   = rep(1., q),
      adglasso =
        sqrt(rowSums((
          ini$coefSVD$v %*%
            diag(ini$coefSVD$d/sqrt(su),
                 nrow = length(ini$coefSVD$d))
        ) ^2)) ^ (-wgamma)
    )
    if (substr(methodA, 1, 2) == "ad" ||
        substr(methodB, 1, 2) == "ad") {
      if (is.null(ini))
        ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
      Wd <- (ini$coefSVD$d/sqrt(su)/sqrt(sv)) ^ (-wgamma)
    } else {
      Wd <- rep(1., nrank)
    }
  }

  ## Compute some products for later use
  if (is.null(XX))
    XX <- crossprod(X)    # t(X)%*%X
  if (is.null(eigenXX))
    eigenXX <- eigen(XX)$values[1]
  if (is.null(XY))
    XY <- crossprod(X, Y) # t(X)%*%Y
  if (is.null(Yvec))
    Yvec <- as.vector(Y)


  for (i in 1:maxit) {
    ##    while(diff[iter] > conv & iter < maxit){
    ## Store initial values before iterations
    ## update rank
    nzid <- which(D != 0)
    ## If zero solution occurs
    if (length(nzid) == 0)
      nzid <- 1  ## JY: WHY?

    nrank <- length(nzid)
    U <- U0 <- as.matrix(U[, nzid, drop = FALSE])
    V <- V0 <- as.matrix(V[, nzid, drop = FALSE])
    D <- D0 <- D[nzid] # , nzid])
    ## d <- d0 <- d[nzid]
    A <- A0 <- as.matrix(A[, nzid, drop = FALSE])
    B <- B0 <- as.matrix(B[, nzid, drop = FALSE])
    GA <- GA0 <- as.matrix(GA[, nzid, drop = FALSE])
    GB <- GB0 <- as.matrix(GB[, nzid, drop = FALSE])
    WA <- switch(
      methodA,
      lasso = as.matrix(WA[, nzid]),
      adlasso = as.matrix(WA[, nzid]),
      glasso = as.vector(WA),
      adglasso = as.vector(WA)
    )
    WB <-  switch(
      methodB,
      lasso = as.matrix(WB[, nzid]),
      adlasso = as.matrix(WB[, nzid]),
      glasso = as.vector(WB),
      adglasso = as.vector(WB)
    )
    Wd <- Wd[nzid]
    ee <- diag(nrow = nrank)
    su <- su[nzid]
    sv <- sv[nzid]

    ## L2 step -------------------------------------------------------------
    ## 1. U step
    ## Note that only XY and XX are needed
    XXmu <- XX + mu * diag(nrow = p)
    XYu <- XY %*% V + mu * (A - GA)
    XXu <- XXmu
    Du <- D
    ## Ustep <- sofar.procrustes(XYu, XXu,
    ##                           diag(Du,nrow=length(Du),
    ##                                ncol=length(Du)),
    ## rho2 = eigenXX + mu, U,
    ## control = list(maxit=innerMaxit,
    ##                epsilon = innerEpsilon))
    Ustep <- procrustes_RCpp(
      XYu,
      XXu,
      diag(Du*sqrt(su), nrow = length(Du), ncol = length(Du)),
      rho2 = eigenXX + mu,
      t((1/sqrt(su))*t(U)),
      control = list(maxit = innerMaxit,
                     epsilon = innerEpsilon)
    )
    U <- t(sqrt(su)*t(Ustep$U))

    ## 2. V step
    tmp <- crossprod(XY, U) + mu * (B - GB)
    tmp <- t(D*sqrt(sv)* t(tmp))
    XYvsvd <- svd(tmp)
    V <- t(sqrt(sv)*tcrossprod(XYvsvd$v,XYvsvd$u))

    ## 3. D step
    Xd1 <- matrix(0., n * q, nrank)
    Xd2 <- matrix(0., p * nrank, nrank)
    Xd3 <- matrix(0., q * nrank, nrank)

    XU <- X %*% U
    for (k in 1:nrank) {
      Xd1[, k] <-
        kron_RcppArma(V[, k, drop = FALSE], XU[, k, drop = FALSE])
      ## ek <- as.vector(rep(0., nrank))
      ## ek[k] <- 1
      Xd2[, k] <-
        kron_RcppArma(U[, k, drop = FALSE], ee[, k, drop = FALSE]) # ek)
      Xd3[, k] <-
        kron_RcppArma(V[, k, drop = FALSE], ee[, k, drop = FALSE]) # ek)
    }


    Xd <-
      crossprod(Xd1) + mu * (crossprod(Xd2) + crossprod(Xd3))
    Yd <-
      crossprod(Xd1, Yvec) + mu * (t(Xd2) %*% as.vector(t(A - GA)) +
                                     t(Xd3) %*% as.vector(t(B - GB)))
    if (lamD == 0) {
      D <- solve(Xd, Yd)
    } else {
      ## Wdminv <- solve(diag(Wd, nrow=nrank,ncol=nrank))
      Wdminv <- 1. / Wd  ## a vector
      Xdw <- Wdminv * Xd %*% diag(Wdminv, nrow = length(Wdminv))
      Ydw <- Wdminv * Yd
      D <- Wdminv * as.vector(lassoshooting(
        XtX = Xdw,
        Xty = Ydw,
        lambda = lamD
      )$coef)
    } ## use penreg in place lassoshooting

    ## FIXME: should this be a tolerence?
    D[D < 0] <- 0
    D <- as.vector(D)
    ## D <- as.vector(d) ## diag(nrow=nrank,ncol=nrank, as.vector(d))
    ## If some d becomes zero, update U and V
    U[, which(D == 0)] <- 0
    V[, which(D == 0)] <- 0


    ## Thresholding step ---------------------------------------------------
    ##There might be some errors in the memo
    ##A <- softTH(U%*%D+GA,lamA/mu*WA)
    ##B <- softTH(V%*%D+GB,lamB/mu*WB)

    A <- switch(
      methodA,
      lasso = softTH(U %*% diag(D, nrow = length(D)) +
                       GA, lamA / mu * WA, tol = ifelse(lamA / mu * WA <=tol, 0, tol)),
      adlasso = softTH(U %*% diag(D, nrow = length(D)) +
                         GA, lamA / mu * WA, tol = ifelse(lamA / mu * WA <=tol, 0, tol)),
      glasso = softrowTH(U %*% diag(D, nrow = length(D)) +
                           GA, lamA / mu * WA, tol = ifelse(lamA / mu * WA <=tol, 0, tol))$C,
      adglasso = softrowTH(U %*% diag(D, nrow = length(D)) +
                             GA, lamA / mu * WA, tol = ifelse(lamA / mu * WA <=tol, 0, tol))$C
    )
    B <- switch(
      methodB,
      lasso = softTH(V %*% diag(D, nrow = length(D)) +
                       GB, lamB / mu * WB, tol = ifelse(lamB / mu * WB <=tol, 0, tol)),
      adlasso = softTH(V %*% diag(D, nrow = length(D)) +
                         GB, lamB / mu * WB, tol = ifelse(lamB / mu * WB <=tol, 0, tol)),
      glasso = softrowTH(V %*% diag(D, nrow = length(D)) +
                           GB, lamB / mu * WB, tol = ifelse(lamB / mu * WB <=tol, 0, tol))$C,
      adglasso = softrowTH(V %*% diag(D, nrow = length(D)) +
                             GB, lamB / mu * WB, tol = ifelse(lamB / mu * WB <=tol, 0, tol))$C
    )

    ## Dual step
    GA <-
      (GA + (t(D * t(U)) - A)) / mugamma # (GA + (U %*% D - A)) / mugamma
    GB <-
      (GB + (t(D * t(V)) - B)) / mugamma # (GB + (V %*% D - B)) / mugamma

    mu <- mu * mugamma

    del <- sqrt(sum((A - A0) ^ 2) + sum((B - B0) ^ 2))
    if (del < epsilon){
      #cat("del = ", del , "\n")
      break
    }

  }

  ## Final results
  invD <- D
  invD[D != 0] <- 1 / D[D != 0]
  U <- A %*% diag(invD, nrow = length(invD)) # ginv(D)
  V <- B %*% diag(invD, nrow = length(invD)) # ginv(D)
  d0id <- which(colSums(abs(U)) == 0)
  D[d0id] <- 0
  d0id <- which(colSums(abs(V)) == 0)
  D[d0id] <- 0
  ## some rounding##
  D[abs(D) / sum(D) < tol] <- 0
  #dorder <- order(D, decreasing = TRUE)
  dorder <- 1:length(invD)

  ## D <- as.vector(d) # diag(nrow=nrank, ncol=nrank, as.vector(d[dorder]))

  ## Update estimated rank
  nrank <- sum(D != 0)

  ## Reordering
  U <- as.matrix(A[, dorder] %*% diag(invD, nrow = length(invD)))
  V <- as.matrix(B[, dorder] %*% diag(invD, nrow = length(invD)))

  drank <- max(nrank, 1)
  U <- as.matrix(U[, 1:drank, drop = FALSE])
  V <- as.matrix(V[, 1:drank, drop = FALSE])
  D <- D[1:drank]

  ## Compute infomation criteria
  ## residual <- (Y-X%*%U%*%D%*%t(V))
  residual <- (Y - X %*% U %*% diag(D, nrow = drank) %*% t(V))
  sse <- sum(residual ^ 2)

  dfu0 <-  switch(
    methodA,
    lasso = sum(U != 0),
    adlasso = sum(U != 0),
    glasso = sum(apply(U, 1, function(a)
      sum(a ^ 2)) != 0) * nrank,
    adglasso = sum(apply(U, 1, function(a)
      sum(a ^ 2)) != 0) * nrank
  )
  dfv0 <-  switch(
    methodB,
    lasso = sum(V != 0),
    adlasso = sum(V != 0),
    glasso = sum(apply(V, 1, function(a)
      sum(a ^ 2)) != 0) * nrank,
    adglasso = sum(apply(V, 1, function(a)
      sum(a ^ 2)) != 0) * nrank
  )

  df <- dfu0 * xrank / p + dfv0 - nrank ^ 2

  BIC <-  log(sse) + log(q * n) / q / n * df
  GIC <- log(sse) + log(log(n * q)) * log(p * q) / q / n * df
  AIC <-  log(sse) + 2 / q / n * df
  GCV <- sse / q / n / (1 - df / q / n) ^ 2

  converged <- del <= epsilon

  list(
    # diff=diff,iter=iter,
    ic = c(
      GIC = GIC,
      BIC = BIC,
      AIC = AIC,
      GCV = GCV
    ),
    sse = sse,
    df = df,
    converged = converged,
    U = U,
    V = V,
    D = D,
    rank = nrank
  )
}



sofar.lammax <- function(Y,
                         X,
                         U = NULL,
                         V = NULL,
                         D = NULL,
                         XX = NULL,
                         XY = NULL,
                         Yvec = NULL,
                         eigenXX = NULL,
                         nrank = 3,
                         su = NULL,
                         sv = NULL,
                         ## lamA=0,lamB=0,lamD=0,
                         mu = 1,
                         mugamma = 1.2,
                         conv = 1e-3,
                         maxit = 100,
                         ## inner.conv=1e-3,inner.iter=50,
                         methodA = c("lasso", "adlasso", "glasso", "adglasso")[2],
                         methodB = c("lasso", "adlasso", "glasso", "adglasso")[2],
                         WA = NULL,
                         WB = NULL,
                         Wd = NULL,
                         wgamma = 2)
{
  tol <- 1e-2
  lamA = 0
  lamB = 0
  lamD = 0

  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Y)
  nr <- nrank
  xrank <- sum(svd(X)$d > 1e-4)


  if(is.null(su)) su <- rep(1,nrank)
  #   Su_inv <- diag(nrow=nrank, ncol=nrank)
  #   Su_inv_sqrt <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Su_inv <- diag(sqrt(su),nrow=nrank, ncol=nrank)
  #   Su_inv_sqrt <- diag(1/sqrt(su),nrow=nrank, ncol=nrank)
  # }
  if(is.null(sv)) sv <- rep(1,nrank)
  #   Sv_inv <- diag(nrow=nrank, ncol=nrank)
  #   Sv_inv_sqrt <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Sv_inv <- diag(sqrt(sv),nrow=nrank, ncol=nrank)
  #   Sv_inv_sqrt <- diag(1/sqrt(sv),nrow=nrank, ncol=nrank)
  # }


  if (is.null(U) | is.null(V) | is.null(D)) {
    ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    U <- t(sqrt(su)*t(ini$coefSVD$u))
    V <- t(sqrt(sv)*t(ini$coefSVD$v))
    D <- diag(ini$coefSVD$d/sqrt(su)/sqrt(sv), nrow = length(ini$coefSVD$d))
  } else{
    ini <- NULL
  }

  d <- diag(D)
  A <- U %*% D
  B <- V %*% D
  O <- D
  o <- diag(O)
  GA <- matrix(nrow = p, ncol = nr, 0)
  GB <- matrix(nrow = q, ncol = nr, 0)
  GO <- matrix(nrow = nr, ncol = nr, 0)

  ## Adaptive weights
  if (is.null(WA) | is.null(WB) | is.null(Wd)) {
    if (is.null(ini))
      ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
    WA <- switch(
      methodA,
      lasso = matrix(nrow = p, ncol = nrank, 1),
      adlasso = ((
          ini$coefSVD$u %*%
          diag(ini$coefSVD$d/sqrt(sv), nrow = length(ini$coefSVD$d))
      ) + 1e-8) ^ (- wgamma),
      glasso = rep(1, p),
      adglasso = apply(ini$coefSVD$u
                       %*% diag(ini$coefSVD$d/sqrt(sv),
                                nrow = length(ini$coefSVD$d)),
                       1, function(a)
                         sqrt(sum(a ^ 2))) ^ {
                           -wgamma
                         }
    )
    WB <- switch(
      methodB,
      lasso = matrix(nrow = q, ncol = nrank, 1),
      adlasso = ((
          ini$coefSVD$v %*%
          diag(ini$coefSVD$d/sqrt(su), nrow = length(ini$coefSVD$d))
      ) + 1e-8) ^ (- wgamma),
      glasso = rep(1, q),
      adglasso = apply(ini$coefSVD$v %*%
                         diag(ini$coefSVD$d/sqrt(su),
                              nrow = length(ini$coefSVD$d)),
                       1, function(a)
                         sqrt(sum(a ^ 2))) ^ {
                           -wgamma
                         }
    )
    if (substr(methodA, 1, 2) == "ad" |
        substr(methodB, 1, 2) == "ad") {
      if (is.null(ini))
        ini <- rrr.fit(Y, X, nrank = nrank, coefSVD = TRUE)
      Wd <- (ini$coefSVD$d / sqrt(su) / sqrt(sv) + 1e-8) ^ {
        -wgamma
      }
    } else{
      Wd <- rep(1, nrank)
    }
  }


  ##Compute some products for later use
  if (is.null(XX))
    XX <- crossprod(X)    # t(X)%*%X
  if (is.null(eigenXX))
    eigenXX <- eigen(XX)$values[1]
  if (is.null(XY))
    XY <- crossprod(X, Y) # t(X)%*%Y
  if (is.null(Yvec))
    Yvec <- as.vector(Y)

  lamA.max <- 0
  lamB.max <- 0
  lamD.max <- 0

  iter <- 1
  diff <- vector()
  diff[1] <- conv * 2

  while (diff[iter] > conv && iter < maxit) {
    ## Store initial values before iterations
    U0 <- U
    V0 <- V
    D0 <- D
    d0 <- d
    A0 <- A
    B0 <- B
    O0 <- O
    o0 <- o
    GA0 <- GA
    GB0 <- GB
    GO0 <- GO

    ## L2 step
    ## 1. U step
    ## Note that only XY and XX are needed
    XXmu <- XX + mu * diag(nrow = p)
    ## Yu <- XXsqrtm%*%(XY%*%V+mu*(A-GA))
    XYu <- XY %*% V + mu * (A - GA)
    ##Xu <- XXsqrt
    XXu <- XXmu
    Du <- D

    Ustep <- procrustes_RCpp(XYu,
                             XXu,
                             sqrt(su)*D,
                             rho2 = eigenXX + mu,
                             t((1/sqrt(su))*t(U)),
                             control = list(maxit = 50,
                                            epsilon = 1e-3))
    U <- t(sqrt(su)*t(Ustep$U))

    ## 2. V step
    XYvsvd <- svd((t(XY) %*% U + mu * (B - GB)) %*% (sqrt(sv)*D))
    V <- t(sqrt(sv)*tcrossprod(XYvsvd$v, XYvsvd$u))

    ## 3. D step
    Xd1 <- matrix(nrow = n * q, ncol = nr)
    Xd2 <- matrix(nrow = p * nr, ncol = nr)
    Xd3 <- matrix(nrow = q * nr, ncol = nr)
    Xd4 <- vector()

    XU <- X %*% U
    for (k in 1:nr) {
      Xd1[, k] <- kronecker(V[, k], XU[, k])
      ek <- as.vector(rep(0, nr))
      ek[k] <- 1
      Xd2[, k] <- kronecker(U[, k], ek)
      Xd3[, k] <- kronecker(V[, k], ek)
    }
    Xd4 <- as.vector(diag(O - GO))

    Xd <-
      t(Xd1) %*% Xd1 + mu * (t(Xd2) %*% Xd2 + t(Xd3) %*% Xd3 + diag(nrow = nrank))
    Yd <- t(Xd1) %*% Yvec +
      mu * (t(Xd2) %*% as.vector(t(A - GA)) +
              t(Xd3) %*% as.vector(t(B - GB)) + as.vector(diag(O - GO)))
    d <- solve(Xd) %*% Yd
    d[d < tol] <- 0
    D <- diag(nrow = nr, ncol = nr, as.vector(d))

    ## If some d becomes zero, update U and V
    U[, which(d == 0)] <- 0
    V[, which(d == 0)] <- 0

    ## Thresholding step
    A <- switch(
      methodA,
      lasso = softTH(U %*% D + GA, lamA / mu * WA),
      adlasso = softTH(U %*% D + GA, lamA / mu * WA),
      glasso = softrowTH(U %*% D + GA, lamA / mu * WA)$C,
      adglasso = softrowTH(U %*% D + GA, lamA / mu * WA)$C
    )
    ##Estimate lamA.max
    lamA.max <- switch(
      methodA,
      lasso = max(abs((U %*% D + GA) * mu / WA)),
      adlasso = max(abs((U %*% D + GA) * mu / WA)),
      glasso = max(abs(lnorm(U %*% D + GA, 2) * mu /
                         WA)),
      adglasso = max(abs(lnorm(U %*% D + GA, 2) * mu /
                           WA))
    )
    B <- switch(
      methodB,
      lasso = softTH(V %*% D + GB, lamB / mu * WB),
      adlasso = softTH(V %*% D + GB, lamB / mu * WB),
      glasso = softrowTH(V %*% D + GB, lamB / mu * WB)$C,
      adglasso = softrowTH(V %*% D + GB, lamB / mu * WB)$C
    )
    lamB.max <- switch(
      methodB,
      lasso = max(abs((V %*% D + GB) * mu / WB)),
      adlasso = max(abs((V %*% D + GB) * mu / WB)),
      glasso = max(abs(lnorm(V %*% D + GB, 2) * mu /
                         WB)),
      adglasso = max(abs(lnorm(V %*% D + GB, 2) * mu /
                           WB))
    )

    o <- softTH(diag(D + GO), lamD / mu * Wd)
    lamD.max <- max(abs(diag(D + GO) * mu / Wd))

    o[o < 0] <- 0
    O <- diag(nrow = nr, ncol = nr, o)

    ## Dual step
    GA <- (GA + (U %*% D - A)) / mugamma
    GB <- (GB + (V %*% D - B)) / mugamma
    GO <- (GO + (D - O)) / mugamma

    mu <- mu * mugamma

    iter <- iter + 1
    ## diff <- sqrt(sum((A%*%t(B)-A0%*%t(B0))^2))
    diff[iter] <-
        sqrt(sum((A - A0) ^ 2) + sum((B - B0) ^ 2) + sum((O - O0) ^ 2))
  }

  list(lamA.max = lamA.max,
       lamB.max = lamB.max,
       lamD.max = lamD.max)

}


sofar.path.reduced <- function(Y,
                               X,
                               XX = NULL,
                               XY = NULL,
                               Yvec = NULL,
                               eigenXX = NULL,
                               nrank = 1,
                               su = NULL,
                               sv = NULL,
                               lamA = 1,
                               lamB = 1,
                               lamD = 1,
                               ic = c("GIC", "BIC", "AIC", "GCV"),
                               modstr = list(),
                               init = list(U = NULL, V = NULL, D = NULL),
                               control = list()) {
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  nlamA <- length(lamA)
  nlamB <- length(lamB)
  nlamD <- length(lamD)
  stopifnot(nlamA == nlamB & nlamB == nlamD)

  if(is.null(su)) su <- rep(1,nrank)
  #   Su_inv <- diag(nrow=nrank, ncol=nrank)
  #   Su_inv_sqrt <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Su_inv <- diag(sqrt(su),nrow=nrank, ncol=nrank)
  #   Su_inv_sqrt <- diag(1/sqrt(su),nrow=nrank, ncol=nrank)
  # }
  if(is.null(sv)) sv <- rep(1,nrank)
  #   Sv_inv <- diag(nrow=nrank, ncol=nrank)
  #   Sv_inv_sqrt <- diag(nrow=nrank, ncol=nrank)
  # }else{
  #   Sv_inv <- diag(sqrt(sv),nrow=nrank, ncol=nrank)
  #   Sv_inv_sqrt <- diag(1/sqrt(sv),nrow=nrank, ncol=nrank)
  # }

  modstr <- do.call("sofar.modstr", modstr)
  control <- do.call("sofar.control", control)

  ## nrank <- modstr$nrank

  nlam <- nlamA

  ICpath <- array(dim = c(nlam, 4), NA)
  Upath <- array(dim = c(nlam, p, nrank), 0)
  Vpath <- array(dim = c(nlam, q, nrank), 0)
  Dpath <- array(dim = c(nlam, nrank), 0)
  Rpath <- rep(0, nlam)

  ## initial value
  ##init <- do.call("sofar.0", init)
  U <- init$U
  V <- init$V
  D <- init$D

  ini <- NULL
  ## Adaptive weights
  WA <-
    modstr$WA
  WB <- modstr$WB
  Wd <- modstr$Wd
  wgamma <- modstr$wgamma
  methodA <- control$methodA
  methodB <- control$methodB

  if (is.null(WA) || is.null(WB) || is.null(Wd)) {
    if (is.null(ini))
      ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
    WA <- switch(
      methodA,
      lasso    = matrix(1., p, nrank),
      adlasso  = (ini$coefSVD$u %*%
                    diag(ini$coefSVD$d/sqrt(sv),
                         nrow = length(ini$coefSVD$d))) ^ (-wgamma),
      glasso   = rep(1., p),
      adglasso = sqrt(rowSums((
        ini$coefSVD$u %*%
          diag(ini$coefSVD$d/sqrt(sv),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )
    WB <- switch(
      methodB,
      lasso    = matrix(1., q, nrank),
      adlasso  = (ini$coefSVD$v %*%
                    diag(ini$coefSVD$d/sqrt(su),
                         nrow = length(ini$coefSVD$d))) ^ (-wgamma),
      glasso   = rep(1., q),
      adglasso = sqrt(rowSums((
        ini$coefSVD$v %*%
          diag(ini$coefSVD$d/sqrt(su),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )
    if (substr(methodA, 1, 2) == "ad" ||
        substr(methodB, 1, 2) == "ad") {
      if (is.null(ini))
        ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
      Wd <- (ini$coefSVD$d/sqrt(su)/sqrt(sv)) ^ (-wgamma)
    } else {
      Wd <- rep(1., nrank)
    }
  }

  ## initial
  Upath[1, ,] <- U
  Vpath[1, ,] <- V
  Dpath[1,] <- D
  Rpath[1] <- nrank

  for (h in 1:nlam) {
    hidx <- max(h - 1, 1)
    U <- as.matrix(Upath[hidx, ,])
    V <- as.matrix(Vpath[hidx, ,])
    D <- Dpath[hidx,]
    rfit <- Rpath[hidx]

    U <- as.matrix(U[, 1:rfit])
    V <- as.matrix(V[, 1:rfit])
    D <- D[1:rfit]

    ## weights according to the rank
    WA1 <- switch(
      methodA,
      lasso = as.matrix(WA[, 1:rfit]),
      adlasso = as.matrix(WA[, 1:rfit]),
      glasso = as.vector(WA),
      adglasso = as.vector(WA)
    )
    WB1 <-  switch(
      methodB,
      lasso = as.matrix(WB[, 1:rfit]),
      adlasso = as.matrix(WB[, 1:rfit]),
      glasso = as.vector(WB),
      adglasso = as.vector(WB)
    )
    Wd1 <- Wd[1:rfit]

    ## modstr$nrank <- rfit
    modstr$WA <- WA1
    modstr$WB <- WB1
    modstr$Wd <- Wd1
    control$lamA <- lamA[h]
    control$lamB <- lamB[h]
    control$lamD <- lamD[h]

    ## Use warm start
    fit <- sofar.fit(
      Y,
      X,
      XX = XX,
      XY = XY,
      Yvec = Yvec,
      eigenXX = eigenXX,
      nrank  = rfit,
      su = su,
      sv = sv,
      modstr = modstr,
      init = list(U = U, V = V, D = D),
      control = control
    )
    U <- fit$U
    V <- fit$V
    D <- fit$D
    rfit <- fit$rank

    if (rfit == 0) {
      ICpath[h:nlam, ] <- matrix(fit$ic,
                                 nrow = nlam - h + 1,
                                 ncol = 4,
                                 byrow = TRUE)
      Rpath[h:nlam] <- 0
      break
    }
    ICpath[h, ] <- fit$ic
    Upath[h, , seq_len(rfit)] <- U
    Vpath[h, , seq_len(rfit)] <- V
    Dpath[h, seq_len(rfit)] <- D
    Rpath[h] <- rfit
  }

  minidx <- apply(ICpath, 2, which.min)

  IC <- match.arg(ic)

  ## Which IC to use?
  ICid <- switch(
    IC,
    GIC = minidx[1],
    BIC = minidx[2],
    AIC = minidx[3],
    GCV = minidx[4]
  )
  ## Note that minid are reported based on the original tuning sequences.

  list(
    Upath = Upath,
    Vpath = Vpath,
    Dpath = Dpath,
    Rpath = Rpath,
    ICpath = ICpath,
    lam.id = minidx,
    U = Upath[ICid, , 1:Rpath[ICid]],
    V = Vpath[ICid, , 1:Rpath[ICid]],
    D = Dpath[ICid, 1:Rpath[ICid]],
    rank = Rpath[ICid]
  )
}



sofar.cv <- function(Y,
                     X,
                     XX = NULL,
                     XY = NULL,
                     Yvec = NULL,
                     eigenXX = eigenXX,
                     nrank = 1,
                     su = NULL,
                     sv = NULL,
                     lamA = 1,
                     lamB = 1,
                     lamD = 1,
                     nfold = 5,
                     norder = NULL,
                     #fold.drop = 0,
                     modstr = list(),
                     init = list(U = NULL, V = NULL, D = NULL),
                     control = list())
{
  p <- ncol(X)
  q <- ncol(Y)
  n <- nrow(Y)

  nlamA <- length(lamA)
  nlamB <- length(lamB)
  nlamD <- length(lamD)
  stopifnot(nlamA == nlamB & nlamB == nlamD)
  nlam <- nlamA

  if(is.null(sv)) su <- rep(1,nrank)
  if(is.null(sv)) sv <- rep(1,nrank)

  modstr <- do.call("sofar.modstr", modstr)
  ## init <- do.call("sofar.0", init)
  control <- do.call("sofar.control", control)

  ## Adaptive weights are the same for different fold
  WA <- modstr$WA
  WB <- modstr$WB
  Wd <- modstr$Wd
  wgamma <- modstr$wgamma
  methodA <- control$methodA
  methodB <- control$methodB
  ini <- NULL
  if (is.null(WA) || is.null(WB) || is.null(Wd)) {
    if (is.null(ini))
      ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
    WA <- switch(
      methodA,
      lasso = matrix(1., p, nrank),
      adlasso = (ini$coefSVD$u %*%
                   diag(ini$coefSVD$d/sqrt(sv),
                        nrow = length(ini$coefSVD$d))) ^ (-wgamma),
      glasso = rep(1., p),
      adglasso = sqrt(rowSums((
        ini$coefSVD$u %*%
          diag(ini$coefSVD$d/sqrt(sv),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )
    WB <- switch(
      methodB,
      lasso = matrix(1., q, nrank),
      adlasso = (ini$coefSVD$v %*%
                   diag(ini$coefSVD$d/sqrt(su),
                        nrow = length(ini$coefSVD$d))) ^ (-wgamma),
      glasso = rep(1., q),
      adglasso = sqrt(rowSums((
        ini$coefSVD$v %*%
          diag(ini$coefSVD$d/sqrt(su),
               nrow = length(ini$coefSVD$d))
      ) ^ 2)) ^ (-wgamma)
    )
    if (substr(methodA, 1, 2) == "ad" ||
        substr(methodB, 1, 2) == "ad") {
      if (is.null(ini))
        ini <- rrr.fit(Y, X, nrank, coefSVD = TRUE)
      Wd <- (ini$coefSVD$d/sqrt(su)/sqrt(sv)) ^ (-wgamma)
    } else {
      Wd <- rep(1., nrank)
    }
    modstr$WA <- WA
    modstr$WB <- WB
    modstr$Wd <- Wd
  }
  if (is.null(init$U) || is.null(init$V) || is.null(init$D)) {
    init <- list(
      U = t(sqrt(su)*t(ini$coefSVD$u)),
      V = t(sqrt(sv)*t(ini$coefSVD$v)),
      D = ini$coefSVD$d/sqrt(su)/sqrt(sv)
    )
  }


  ## Construct folds
  ndel <- round(n / nfold)
  if (is.null(norder))
    norder <- sample(n)

  nlam <- length(lamA)
  cr_path <- array(NA, dim = c(nlam, nfold))

  for (f in seq_len(nfold)) {
    ## determine cr sample
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    } else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    ## idkeep <- (1:n)[-iddel]

    Xf <- X[-iddel, ]
    Xfdel <- X[iddel, ]
    Yf <- Y[-iddel, ]
    Yfdel <- Y[iddel, ]

    XXf <- crossprod(Xf)    # t(X)%*%X
    XYf <- crossprod(Xf, Yf) # t(X)%*%Y
    Yvecf <- as.vector(Yf)
    eigenXXf <- eigen(XXf)$values[1]

    fitf <- sofar.path.reduced(
      Yf,
      Xf,
      XX = XXf,
      XY = XYf,
      Yvec = Yvecf,
      eigenXX = eigenXXf,
      nrank = nrank,
      su = su,
      sv = sv,
      lamA = lamA,
      lamB = lamB,
      lamD = lamD,
      ic = "GIC",
      modstr = modstr,
      control = control,
      init = init
    )

    for (ll in 1:nlam) {
      cr_path[ll, f] <- sum((
        Yfdel - Xfdel %*% fitf$Upath[ll, ,] %*%
          diag(fitf$Dpath[ll,], nrow = nrank) %*% t(fitf$Vpath[ll, ,])
      ) ^ 2)
    }

  }

  index <- order(apply(cr_path, 2, sum))
  crerr <- apply(cr_path[, index], 1, sum) / length(index) * nfold
  minid <- which.min(crerr)

  ## Refit model with all data.
  ## modstr$WA <- WA; modstr$WB <- WB; modstr$Wd <- Wd
  control$lamA <- lamA[minid]
  control$lamB <- lamB[minid]
  control$lamD <- lamD[minid]

  ## Initial values
  ## init <- list(U=fitf$Upath[minid,,],
  ##              V=fitf$Vpath[minid,,],
  ##              D=fitf$Dpath[minid,])
  fitall <- sofar.path.reduced(
    Y,
    X,
    XX = XX,
    XY = XY,
    Yvec = Yvec,
    eigenXX = eigenXX,
    ic = "GIC",
    nrank = nrank,
    su = su,
    sv = sv,
    lamA = lamA,
    lamB = lamB,
    lamD = lamD,
    modstr = modstr,
    init = init,
    control = control
  )
  ## fitall <- sofar.fit(Y, X, XX = XX, XY = XY,
  ##                     Yvec = Yvec, eigenXX = eigenXX,
  ##                     nrank = nrank, modstr = modstr,
  ##                     control = control, init = init)

  list(
    cr.path = cr_path,
    cr.error = crerr,
    lam.id = c(minid),
    norder = norder,
    ## W = list(WA=WA, WB=WB, Wd=Wd),
    Upath = fitall$Upath,
    Vpath = fitall$Vpath,
    Dpath = fitall$Dpath,
    Rpath = fitall$Rpath,
    U = fitall$Upath[minid, , 1:fitall$rank],
    V = fitall$Vpath[minid, , 1:fitall$rank],
    D = fitall$Dpath[minid, 1:fitall$rank],
    ## U = fitall$Upath[minid,,],
    ## V = fitall$Vpath[minid,,],
    ## D = fitall$Dpath[minid,],
    ## U = fitall$U,
    ## V = fitall$V,
    ## D = fitall$D,
    rank = fitall$rank
  )
}
