##' Generalized or mixed-response reduced-rank regression
##'
##' Peforms either rank constrained maximum likelihood estimation or
##' singular value penalized estimation.
##'
##' The model fitting process can be fine tuned through argument \code{control}.
##' The available elements for \code{control} include
##' \itemize{
##'    \item{epsilon}: positive convergence tolerance epsilon; the
##'        iterations converge when |new - old | / (old + 0.1) < epsilon.
##'        treated as zero.
##'    \item{sv.tol}: tolerance for singular values.
##'    \item{maxit}: integer giving the maximal number of iterations.
##'    \item{trace:}{logical indicating if tracing the objective is needed.}
##'    \item{conv.obj:}{if TRUE, track objective function.}
##'    \item{equal.phi:}{if TRUE, use a single dispersion parameter for Gaussian
##'        responses.}
##'    \item{plot.obj:}{if TRUE, plot obj values along iterations; for checking
##'        only}
##'    \item{plot.cv:}{if TRUE, plot cross validation error.}
##'    \item{gammaC0:}{adaptive scaling to speed up computation.}
##' }
##' Similarly, the available elements for arguments \code{penstr} specifying
##' penalty structure of SVD include
##' \itemize{
##'    \item{penaltySVD}: penalty for reducing rank
##'    \item{lambdaSVD}: tuning parameter. For penaltySVD = rankCon, this is the
##'     specified rank.
##' }
##'
##' @usage
##' mrrr(Y, X, is.pca = NULL, offset = NULL, ctrl.id = c(),
##'      family = list(gaussian(), binomial()),
##'      familygroup = NULL, maxrank = min(ncol(Y), ncol(X)),
##'      penstr = list(), init = list(), control = list())
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param is.pca If TRUE, mixed principal component analysis with X=I
##' @param offset matrix of the same dimension as Y for offset
##' @param ctrl.id indices of unpenalized predictors
##' @param family a list of family functions as used in \code{glm}
##' @param familygroup a list of family indices of the responses
##' @param maxrank integer giving the maximum rank allowed. Usually this can be
##'     set to min(n,p,q)
##' @param penstr a list of penalty structure of SVD, contains penstr$penaltySVD
##'     is the penalty of SVD, penstr$lambdaSVD is the regularization parameter
##' @param init a list of initial values of kappaC0, kappaS0, C0, and S0
##' @param control a list of controling parameters for the fitting
##'
##' @return
##'
##' S3 \code{mrrr} object, a list containing
##' \item{obj}{the objective function tracking}
##' \item{converged}{TRUE/FALSE for convergence}
##' \item{coef}{the estimated coefficient matrix}
##' \item{outlier}{the estimated outlier matrix}
##' \item{nrank}{the rank of the fitted model}
##'
##' @examples
##' library(rrpack)
##' simdata <- rrr.sim3(n = 100, p = 30, q.mix = c(5, 20, 5),
##'                     nrank = 2, mis.prop = 0.2)
##' Y <- simdata$Y
##' Y_mis <- simdata$Y.mis
##' X <- simdata$X
##' X0 <- cbind(1, X)
##' C <- simdata$C
##' family <- simdata$family
##' familygroup <- simdata$familygroup
##' svdX0d1 <- svd(X0)$d[1]
##' init1 = list(kappaC0 = svdX0d1 * 5)
##' offset = NULL
##' control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 2000,
##'                trace = FALSE, gammaC0 = 1.1, plot.cv = TRUE,
##'                conv.obj = TRUE)
##' fit.mrrr <- mrrr(Y_mis, X, family = family, familygroup = familygroup,
##'                  penstr = list(penaltySVD = "rankCon", lambdaSVD = 2),
##'                  control = control, init = init1)
##' summary(fit.mrrr)
##' coef(fit.mrrr)
##' par(mfrow = c(1, 2))
##' plot(fit.mrrr$obj)
##' plot(C ~ fit.mrrr$coef[- 1 ,])
##' abline(a = 0, b = 1)
##' @importFrom stats gaussian binomial
##' @importFrom graphics plot
##' @importFrom MASS ginv
##' @export
mrrr <- function(Y,
                 X,
                 is.pca = NULL,
                 offset = NULL,
                 ctrl.id = c(),
                 family = list(gaussian(), binomial()),
                 familygroup = NULL,
                 maxrank = min(ncol(Y), ncol(X)),
                 penstr = list(),
                 init = list(),
                 control = list())
{
  Call <- match.call()

  control <- do.call("mrrr.control", control)
  penstr <- do.call("mrrr.penstr", penstr)
  init <- do.call("mrrr.init", init)

  conv.obj <- control$conv.obj
  equal.phi <- control$equal.phi
  plot.obj <- control$plot.obj
  gammaC0 <- control$gammaC0

  if (is.null(is.pca)) {
    if (is.null(X)) {
      cat("Doing generalized PCA...\n")
      X = diag(nrow(Y))
      is.pca <- TRUE
    } else if (nrow(X) == ncol(X) &&
               sum(abs(X - diag(nrow(Y)))) < 0.01) {
      is.pca <- TRUE
      cat("Doing generalized PCA...\n")
    } else if (sum(is.na(X)) > 0 ||
               is.null(X) || nrow(X) != nrow(Y)) {
      stop(
        paste(
          "X need to have the same rows as Y,",
          "no missing value!",
          "OR set is.pca=TRUE for generalized PCA!"
        )
      )
    } else
      is.pca <- FALSE
  } else if (is.pca) {
    X = diag(nrow(Y))
    cat("Doing generalized PCA...\n")
  } else if (!is.pca) {
    if (sum(is.na(X)) > 0 || is.null(X) || nrow(X) != nrow(Y))
      stop("X need to have the same rows as Y, no missing value!\n")
  } else {
    stop("is.pca ERROR\n")
  }
  ## get the dimensions
  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  X0 <- cbind(1, X)
  ## put intercept and control variables' coef together, need not to be
  ## penalized
  id <- if (is.null(ctrl.id))
    c(1)
  else
    c(1, ctrl.id + 1)

  ## should also check family and family group
  nfamily <- length(family)
  if (nfamily == 1 &
      is.null(familygroup))
    familygroup <- rep(1, q)
  ## characters of families
  cfamily <- unique(familygroup)

  kappaC0 <- init$kappaC0

  C0 <- init$C0

  Phi0 <- init$Phi0
  if (is.null(kappaC0)) {
    temp <- vector()
    svdX0d1 <- svd(X0)$d[1]
    for (j in 1:nfamily) {
      temp[j] <- switch(
        family[[j]]$family,
        'gaussian' = svdX0d1,
        'binomial' = svdX0d1 / 2,
        'poisson' = svdX0d1 * 10
      )
    }
    kappaC0 <- max(temp)
  }
  kappaC0_ini <- kappaC0

  ## initial value and scaling
  C <-
    if (is.null(C0))
      matrix(0, nrow = p + 1, ncol = q)
  else
    C0 * kappaC0

  ##initialize dispersion parameters
  if (is.null(Phi0))
    Phi0 <- as.vector(rep(1, q))
  Phi <- Phi0

  objFun <- function(Y, mu, Phi) {
    temp <- vector()
    for (j in 1:nfamily) {
      Sig <- if (family[[j]]$family == "gaussian")
        t(matrix(Phi[familygroup == cfamily[j]],
                 sum(familygroup == cfamily[j]), n))
      else
        1
      temp[j] <-
        -sum(logLikehood(Y[, familygroup == cfamily[j]],
                         mu[, familygroup == cfamily[j]],
                         sqrt(Sig),
                         family[[j]]$family),
             na.rm = TRUE)
    }
    sum(temp)
  }

  penaltySVD <- penstr$penaltySVD
  lambdaSVD <- penstr$lambdaSVD
  lambdaSVD <- switch(
    penaltySVD,
    'nn' = lambdaSVD / kappaC0,
    'rankCon' = lambdaSVD,
    'rankPen' = lambdaSVD / kappaC0
  )
  thresholdSVD <- switch(
    penaltySVD,
    'nn' = softThres,
    'rankCon' = countThres,
    'rankPen' = hardThres
  )
  penaltySVD <- switch(
    penaltySVD,
    'nn' = softPen,
    'rankCon' = countPen,
    'rankPen' = hardPen
  )

  diff <- vector()
  obj <- vector()
  obj.inc <- 0
  Phitrace <- matrix(NA, control$maxit, q)
  Phitrace [1, ] <- Phi
  kappaC0Trace <- vector()
  kappaC0Trace[1] <- kappaC0
  converged <- FALSE
  if (is.null(offset))
    offset <- 0
  MU <- matrix(nrow = n, ncol = q, NA)

  if (is.pca) {
    MU0 <- offset + t(t(C[-1, ]) + C[1, ]) / kappaC0
    for (j in 1:nfamily)
      MU[, familygroup == cfamily[j]] <-
        family[[j]]$linkinv(MU0[, familygroup == cfamily[j]])

    if (conv.obj)
      obj[1] <- objFun(Y, MU, Phi0) +
        svdPen(C / kappaC0, maxrank, penaltySVD, lambdaSVD)
    flag <- FALSE
    iter <- 1
    iter_actural <- 1

    ##while (iter_actural <= (control$maxit - 1) ) {
    while (iter < control$maxit &
           iter_actural < control$maxit) {
      C_temp <- C
      MU_temp <- MU
      Phi_temp <- Phi
      kappaC0_temp <- kappaC0

      ## Update Intercept and nonpenalized coefficients
      Res <- Y - MU
      Res[is.na(Res)] <- 0
      C[id,] <- C[id,] + t(t(rbind(apply(Res, 2, sum),
                                   Res[ctrl.id, ])) / Phi) /
        kappaC0
      MU0 <- offset + t(t(C[-1, ]) + C[1, ]) / kappaC0
      for (j in 1:nfamily)
        MU[, familygroup == cfamily[j]] <-
        family[[j]]$linkinv(MU0[, familygroup == cfamily[j]])

      ## Update penalized cofficients
      Res <- Y - MU
      Res[is.na(Res)] <- 0

      C[-id,] <- C[-id,] +
        t(t(if (is.null(ctrl.id))
          Res
          else
            Res[-ctrl.id, ]) / Phi) / kappaC0
      lamt <- if (penstr$penaltySVD == "rankCon")
        penstr$lambdaSVD
      else
        penstr$lambdaSVD / kappaC0
      C[-id,] <-
        svdThres(C[-id,], maxrank, thresholdSVD, lamt)$C

      ## update Phi
      MU0 <- offset + t(t(C[-1, ]) + C[1, ]) / kappaC0
      for (j in 1:nfamily)
        MU[, familygroup == cfamily[j]] <-
        family[[j]]$linkinv(MU0[, familygroup == cfamily[j]])
      Res <- Y - MU
      Res[is.na(Res)] <- 0

      for (j in 1:nfamily)
        Phi[familygroup == cfamily[j]] <-
        if (family[[j]]$family == "gaussian")
          mean(Res[, familygroup == cfamily[j]] ^ 2, na.rm =
                 TRUE)
      else
        1


      Phitrace[iter + 1,] <- Phi

      if (conv.obj)
        obj[iter + 1] <- objFun(Y, MU, Phi) +
        svdPen(C, maxrank, penaltySVD, lambdaSVD / kappaC0)

      kappaC0Trace[iter + 1] <- kappaC0
      diff[iter] <- if (conv.obj)
        abs(obj[iter] - obj[iter + 1]) /
        (abs(obj[iter]) + 0.1)
      else
        sum((C - C_temp) ^ 2) / sum(C_temp ^ 2)

      if (control$trace)
        cat("iter =", iter, "obj/Cnorm_diff =",
            if (conv.obj)
              obj[iter]
            else
              diff[iter], "\n")
      if (diff[iter] < control$epsilon) {
        converged <- TRUE
        break
      }
      ## Decrease kappaC0 if obj is decreasing.
      ## This should also be done separately
      ## Check two iterations apart.
      if (iter > 1 && conv.obj) {
        ## if((obj[iter + 1] < obj[iter])) {
        if ((obj[iter + 1] < obj[ifelse(iter == 1, 1, iter - 1)])) {
          C <- C / gammaC0
          kappaC0 <- kappaC0 / gammaC0
          obj.inc <- 0
        } else {
          ## Return to the previous step and change C
          ## C <- C/kappaC0_temp*kappaC0_ini
          ## cat(kappaC0,"\n")
          C <-
            C_temp / kappaC0_temp * kappaC0_ini #*gammaC0^obj.inc
          kappaC0 <- kappaC0_ini #*gammaC0^obj.inc
          MU <- MU_temp
          Phi <- Phi_temp
          obj.inc <- obj.inc + 1
          iter <- iter - 1
        }
      }

      if (conv.obj && obj.inc > 5) {
        ## if (control$trace)
        ##     cat("Obj stops decreasing. Use larger kappaC0!\n")
        converged <- TRUE
        break
      }
      iter <- iter + 1
      iter_actural <- iter_actural + 1
    }
  } else {
    MU0 <- offset + X0 %*% C / kappaC0
    for (j in 1:nfamily)
      MU[, familygroup == cfamily[j]] <-
        family[[j]]$linkinv(MU0[, familygroup == cfamily[j]])

    if (conv.obj)
      obj[1] <- objFun(Y, MU, Phi0) +
        svdPen(C / kappaC0, maxrank, penaltySVD, lambdaSVD)
    ## flag <- FALSE
    iter <- 1
    iter_actural <- 1
    ## cat("initialized!", " \n" )

    while (iter < control$maxit &
           iter_actural < control$maxit) {
      ## while (iter_actural <= (control$maxit - 1)) {
      C_temp <- C
      MU_temp <- MU
      Phi_temp <- Phi
      kappaC0_temp <- kappaC0
      ## cat("kappaC0 = ", kappaC0, " \n")

      ## Update Intercept and nonpenalized cofficients
      Res <- Y - MU
      Res[is.na(Res)] <- 0

      C[id,] <-
        C[id,] + t(crossprod(Res, X0[, id]) / Phi) / kappaC0
      ## crossprod(X0[, id], Res)%*%diag(1/Phi) / kappaC0
      MU0 <- offset + X0 %*% C / kappaC0
      for (j in 1:nfamily)
        MU[, familygroup == cfamily[j]] <-
        family[[j]]$linkinv(MU0[, familygroup == cfamily[j]])

      ## Update penalized cofficients
      Res <- Y - MU
      Res[is.na(Res)] <- 0

      C[-id,] <-
        C[-id,] + t(crossprod(Res, X0[,-id]) / Phi) / kappaC0
      ## crossprod(X0[, -id], Res)%*%diag(1/Phi) / kappaC0
      lamt <- if (penstr$penaltySVD == "rankCon")
        penstr$lambdaSVD
      else
        penstr$lambdaSVD / kappaC0
      C[-id,] <-
        svdThres(C[-id,], maxrank, thresholdSVD, lamt)$C

      ## update Phi
      MU0 <- offset + X0 %*% C / kappaC0
      for (j in 1:nfamily)
        MU[, familygroup == cfamily[j]] <-
        family[[j]]$linkinv(MU0[, familygroup == cfamily[j]])

      Res <- Y - MU
      Res[is.na(Res)] <- 0

      if (equal.phi) {
        for (j in 1:nfamily)
          Phi[familygroup == cfamily[j]] <-
            if (family[[j]]$family == "gaussian")
              mean(Res[, familygroup == cfamily[j]] ^ 2, na.rm =
                     TRUE)
        else
          1
      } else {
        for (j in 1:nfamily)
          Phi[familygroup == cfamily[j]] <-
            if (family[[j]]$family == "gaussian") {
              apply(Res[, familygroup == j], 2,
                    function(a)
                      mean(a ^ 2, na.rm = TRUE))
            }
        else
          1
      }

      Phitrace[iter + 1,] <- Phi

      if (conv.obj)
        obj[iter + 1] <- objFun(Y, MU, Phi) +
        svdPen(C, maxrank, penaltySVD, lambdaSVD / kappaC0)

      ## Decrease kappaC0 if obj is decreasing.
      ## This should also be done separately
      ## Check two iterations apart.
      if (iter > 1 && conv.obj) {
        ## if((obj[iter + 1] < obj[iter] )) {
        if ((obj[iter + 1] < obj[ifelse(iter == 1, 1, iter - 1)])) {
          C <- C / gammaC0
          kappaC0 <- kappaC0 / gammaC0
          obj.inc <- 0
        } else {
          ## Return to the previous step and change C
          ## C <- C/kappaC0_temp*kappaC0_ini
          C <- C_temp / kappaC0_temp * kappaC0_ini
          kappaC0 <- kappaC0_ini    #*2^obj.inc
          MU <- MU_temp
          Phi <- Phi_temp
          obj.inc <- obj.inc + 1
          iter <- iter - 1
        }
      }

      kappaC0Trace[iter + 1] <- kappaC0
      diff[iter] <- if (conv.obj)
        abs(obj[iter] - obj[iter + 1]) /
        (abs(obj[iter]) + 0.1)
      else
        sum((C - C_temp) ^ 2) / sum(C_temp ^ 2)

      if (control$trace)
        cat("iter = ", iter, " obj/Cnorm_diff = ",
            if (conv.obj)
              obj[iter]
            else
              diff[iter], " \n")

      if (diff[iter] < control$epsilon) {
        converged <- TRUE
        break
      }
      if (conv.obj && obj.inc > 5) {
        ## if (control$trace)
        ##     cat("Obj stops decreasing. Use larger kappaC0!\n")
        converged <- TRUE
        break
      }
      iter <- iter + 1
      iter_actural <- iter_actural + 1
    }
  }
  if (!conv.obj)
    obj[iter] <- objFun(Y, MU, Phi) +
    svdPen(C, maxrank, penaltySVD, lambdaSVD / kappaC0)

  nrank <- sum(svd(C[-id, ])$d / kappaC0 > control$sv.tol)

  ## if(obj[length(obj)] != min(obj)) flag <- TRUE
  ## flag == TRUE &
  if (plot.obj) {
    ## warning("Objective not monotone decreasing; increase kappa!")
    ## warning("Check convergence")
    if (conv.obj)
      plot(obj, main = paste0("iter = ", iter, ", log10(diff) = ",
                              round(log10(diff[iter]), 3)))
    else
      plot(log10(diff), main = paste0("iter = ", iter, ", log10(diff) = ",
                                      round(log10(diff[iter]), 3)))
  }

  rval <- list(
    call = Call,
    iter = iter,
    iter_actural = iter_actural,
    obj = obj,
    converged = converged,
    diff = diff,
    ctrl.id = id,
    family = family,
    Phitrace = Phitrace,
    kappaC0Trace = kappaC0Trace,
    penaltySVD = penaltySVD,
    lambdaSVD = lambdaSVD,
    mu = MU,
    coef = C / kappaC0,
    ## contains the intercept
    dispersion = Phi,
    rank = nrank
  )
  class(rval) <- "mrrr"
  rval
}



##' Mixed-response reduced-rank regression with rank selected by
##' cross validation
##'
##' @usage
##' cv.mrrr(Y, X, is.pca = NULL, offset = NULL, ctrl.id = c(),
##'         family = list(gaussian(), binomial(), poisson()),
##'         familygroup = NULL, maxrank = min(ncol(Y), ncol(X)),
##'         penstr = list(), init = list(), control = list(), nfold = 5,
##'         foldid = NULL, nlam = 20, warm = FALSE)
##'
##' @param Y response matrix
##' @param X covariate matrix
##' @param is.pca If TRUE, mixed principal component analysis with X=I
##' @param offset matrix of the same dimension as Y for offset
##' @param ctrl.id indices of unpenalized predictors
##' @param family a list of family functions as used in \code{glm}
##' @param familygroup a list of family indices of the responses
##' @param maxrank integer giving the maximum rank allowed.
##' @param penstr a list of penalty structure of SVD.
##' @param init a list of initial values of kappaC0, kappaS0, C0, and S0
##' @param control a list of controling parameters for the fitting
##' @param nfold number of folds in cross validation
##' @param foldid to specify the folds if desired
##' @param nlam number of tuning parameters; not effective when using rank
##'     constrained estimation
##' @param warm if TRUE, use warm start in fitting the solution paths
##' @return S3 \code{mrrr} object, a list containing \item{fit}{the output from
##'     the selected model} \item{dev}{deviance measures}
##' @examples
##' \dontrun{
##' library(rrpack)
##' simdata <- rrr.sim3(n = 100, p = 30, q.mix = c(5, 20, 5),
##'                     nrank = 2, mis.prop = 0.2)
##' Y <- simdata$Y
##' Y_mis <- simdata$Y.mis
##' X <- simdata$X
##' X0 <- cbind(1,X)
##' C <- simdata$C
##' family <- simdata$family
##' familygroup <- simdata$familygroup
##' svdX0d1 <- svd(X0)$d[1]
##' init1 = list(kappaC0 = svdX0d1 * 5)
##' offset = NULL
##' control = list(epsilon = 1e-4, sv.tol = 1e-2, maxit = 2000,
##'                trace = FALSE, gammaC0 = 1.1, plot.cv = TRUE,
##'                conv.obj = TRUE)
##' fit.cv.mrrr <- cv.mrrr(Y_mis, X, family = family,
##'                        familygroup = familygroup,
##'                        maxrank = 20,
##'                        penstr = list(penaltySVD = "rankCon",
##'                                      lambdaSVD = c(1 : 6)),
##'                        control = control, init = init1,
##'                        nfold = 10, nlam = 50)
##' summary(fit.cv.mrrr)
##' coef(fit.cv.mrrr)
##' fit.mrrr <- fit.cv.mrrr$fit
##'
##' ## plot(svd(fit.mrrr$coef[- 1,])$d)
##' plot(C ~ fit.mrrr$coef[- 1, ])
##' abline(a = 0, b = 1)
##' }
##' @importFrom stats gaussian binomial poisson
##' @importFrom graphics abline plot
##' @export
cv.mrrr <-
  function(Y,
           X,
           is.pca = NULL,
           offset = NULL,
           ctrl.id = c(),
           family = list(gaussian(), binomial(), poisson()),
           familygroup = NULL,
           maxrank = min(ncol(Y), ncol(X)),
           penstr = list(),
           init = list(),
           control = list(),
           nfold = 5,
           foldid = NULL,
           nlam = 20,
           warm = FALSE)
  {
    Call <- match.call()

    control <- do.call("mrrr.control", control)
    penstr <- do.call("mrrr.penstr", penstr)
    init <- do.call("mrrr.init", init)

    ## conv.obj <- control$conv.obj
    ## equal.phi <- control$equal.phi
    plot.cv <- control$plot.cv
    ## gammaC0 <- control$gammaC0

    ## get the dimensions
    q <- ncol(Y)
    n <- nrow(Y)
    p <- ncol(X)
    X0 <- cbind(1, X)

    ## distribution families ##
    nfamily <- length(family)
    if (nfamily == 1 &
        is.null(familygroup))
      familygroup <- rep(1, q)
    ## characters of families
    cfamily <- unique(familygroup)

    if (length(penstr$lambdaSVD) < 2) {
      cat("No lambdaSVD sequence provided, use default!\n")
      if (penstr$penaltySVD == "rankCon") {
        penstr$lambdaSVD <- c(1:maxrank)
      } else {
        yt <- rep(0, n)
        for (i in 1:nfamily) {
          yt <- yt +
            switch(
              family[[i]]$family,
              'gaussian' = apply(abs(2 * Y[, familygroup == i]),
                                 1, function(a)
                                   sum(a, na.rm = TRUE)),
              'binomial' = apply(abs(2 - Y[, familygroup == i]),
                                 1, function(a)
                                   sum(a, na.rm = TRUE)),
              'poisson' = apply(abs(2 * Y[, familygroup == i] -
                                      2),
                                1, function(a)
                                  sum(a, na.rm = TRUE))
            )
        }
        lambdaSVDmax <- sum(yt * apply(abs(X), 1, max))
        penstr$lambdaSVD <- 10 ^ (seq(log10(lambdaSVDmax * 0.0001),
                                      log10(lambdaSVDmax), len = nlam))
      }
    }

    nlambdaSVD <- length(penstr$lambdaSVD)
    if (is.null(offset))
      offset <- 0

    objFun <- function(Y, mu, Phi) {
      temp <- vector()
      for (j in 1:nfamily) {
        Sig = if (family[[j]]$family == "gaussian")
          t(matrix(Phi[familygroup == cfamily[j]],
                   sum(familygroup == cfamily[j]), n))
        else
          1
        temp[j] <-
          -sum(logLikehood(Y[, familygroup == cfamily[j]],
                           mu[, familygroup == cfamily[j]],
                           sqrt(Sig),
                           family[[j]]$family),
               na.rm = TRUE)
      }
      sum(temp)
    }


    N.nna <- sum(!is.na(Y))
    ind.nna <- which(!is.na(Y))
    ## store the deviance of the test data
    dev <- matrix(NA, nfold, nlambdaSVD)
    conv <- matrix(NA, nfold, nlambdaSVD)
    N.iter <- matrix(NA, nfold, nlambdaSVD)

    if (is.null(foldid)) {
      size <- floor(N.nna / nfold)
      ID <- rep(1:nfold, len = N.nna)
      ID <- sample(ID, N.nna, replace = FALSE)
    } else {
      ID <- foldid
    }

    lambdaSVD1n <- if (penstr$penaltySVD == "rankCon")
      nlambdaSVD:1
    else
      1:nlambdaSVD
    penstr_i = penstr
    penstr_i$lambdaSVD = penstr$lambdaSVD[lambdaSVD1n[1]]

    fit.ini <- mrrr(
      Y,
      X,
      is.pca = is.pca,
      offset = offset,
      ctrl.id = ctrl.id,
      family = family,
      familygroup = familygroup,
      maxrank = maxrank,
      penstr = penstr_i,
      init = init,
      control = control
    )

    C0.tt <- fit.ini$coef

    for (ifold in 1:nfold) {
      ind.test <- ind.nna[which(ID == ifold)]
      Y_test <- Y
      Y_test[-ind.test] <- NA
      Y_train <- Y
      Y_train[ind.test] <- NA

      C0.t <- C0.tt
      for (lambdaSVDnum in lambdaSVD1n) {
        lambdaSVD <- penstr$lambdaSVD[lambdaSVDnum]
        init_i = init
        if (warm)
          init_i$C0 = C0.t
        penstr_i = penstr
        penstr_i$lambdaSVD = lambdaSVD

        fit1 <- mrrr(
          Y_train,
          X,
          is.pca = is.pca,
          offset = offset,
          ctrl.id = ctrl.id,
          family = family,
          familygroup = familygroup,
          maxrank = maxrank,
          penstr = penstr_i,
          init = init_i,
          control = control
        )
        mu.test <- fit1$mu
        mu.test[-ind.test] <- NA
        C0.t <- fit1$coef
        conv[ifold, lambdaSVDnum] <- fit1$conv
        N.iter[ifold, lambdaSVDnum] <- fit1$iter
        dev[ifold, lambdaSVDnum] <- 2 * objFun(Y_test, mu.test,
                                               fit1$dispersion)
      }
    }

    dev.mean <- apply(dev, 2, mean)
    lambdaSVDnum <- which.min(dev.mean)

    if (plot.cv) {
      plot(
        dev.mean,
        xlab = "lambdaSVDnum",
        ylab = "Deviance",
        main = paste(
          "CV, lam.opt=",
          round(penstr$lambdaSVD[lambdaSVDnum], 3),
          ", conv=",
          mean(conv)
        )
      )
      abline(v = lambdaSVDnum)
    }

    penstr.opt = penstr
    penstr.opt$lambdaSVD = penstr$lambdaSVD[lambdaSVDnum]

    if (warm)
      init$C0 <- C0.tt
    fit1 <-
      mrrr(
        Y,
        X,
        is.pca = is.pca,
        offset = offset,
        ctrl.id = ctrl.id,
        family = family,
        familygroup = familygroup,
        maxrank = maxrank,
        penstr = penstr.opt,
        init = init,
        control = control
      )

    out <- list(
      call = Call,
      foldid = ID,
      fit.ini = fit.ini,
      fit = fit1,
      dev = dev,
      conv = conv,
      N.iter = N.iter,
      dev.mean = dev.mean,
      lambdaSVD.opt = penstr$lambdaSVD[lambdaSVDnum],
      lambdaSVD = penstr$lambdaSVD
    )
    class(out) <- "cv.mrrr"
    out
  }



### internal functions =========================================================

## Auxiliary for controlling mrrr
##
## Auxiliary function for \code{mrrr}. Typically used by default
## but may be constructed to fine tune the fitting.
##
## @param epsilon positive convergence tolerance epsilon; the iterations
##     converge when |new - old | / (old + 0.1) < epsilon.
##     treated as zero
## @param sv.tol tolerance for singular values
## @param maxit integer giving the maximal number of iterations.
## @param trace logical indicating if tracing the objective is needed.
## @param conv.obj if TRUE, track objective function.
## @param equal.phi if TRUE, use a single dispersion parameter for Gaussian
##     responses.
## @param plot.obj if TRUE, plot obj values along iterations; for checking only
## @param plot.cv if TRUE, plot cross validation error.
## @param gammaC0 adaptive scaling to speed up computation.
## @return a list with components names as the arguments.
mrrr.control <- function(epsilon = 1e-6,
                         sv.tol = 1e-6,
                         maxit = 3000,
                         trace = FALSE,
                         conv.obj = FALSE,
                         equal.phi = FALSE,
                         plot.obj = FALSE,
                         plot.cv = TRUE,
                         gammaC0 = 1.05)
{
  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(sv.tol) || sv.tol <= 0)
    stop("value of 'sv.tol' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  if (!is.numeric(gammaC0) || gammaC0 < 1)
    stop("value of 'gammaC0' must be > 1")

  list(
    epsilon = epsilon,
    sv.tol = sv.tol,
    maxit = maxit,
    trace = trace,
    conv.obj = conv.obj,
    equal.phi = equal.phi,
    plot.obj = plot.obj,
    plot.cv = plot.cv,
    gammaC0 = gammaC0
  )
}


## Auxiliary for controlling mrrr
##
## Auxiliary function for \code{mrrr}. Typically used by default
## but may be constructed to fine tune the fitting.
##
## @param penaltySVD penalty for reducing rank
## @param lambdaSVD tuning parameter. For penaltySVD = rankCon, this is the
##     specified rank.
mrrr.penstr <- function(penaltySVD = c('rankCon', 'nn', 'rankPen'),
                        lambdaSVD  = 0)
{
  ## lambdaSVD = 0,
  ## For 'rankCon', it is the specified rank.
  ## lambdaS = 0,
  ## For 'rankCon', it is the specified number of outliers.
  penaltySVD <- match.arg(penaltySVD)
  #penaltyS <- match.arg(penaltyS)
  stopifnot(lambdaSVD >= 0 || lambdaSVD >= 0)
  stopifnot(lambdaSVD >= 0)
  list(penaltySVD = penaltySVD,
       ## penaltyS   = penaltyS,
       lambdaSVD  = lambdaSVD)
  ## lambdaS    = lambdaS)
}


mrrr.init <- function(kappaC0 = NULL,
                      C0 = NULL,
                      kappaS0 = NULL,
                      S0 = NULL,
                      Phi0 = NULL)
{
  list(
    kappaC0 = kappaC0,
    C0 = C0,
    kappaS0 = kappaS0,
    S0 = S0,
    Phi0 = Phi0
  )
}


##' @importFrom stats dnorm dpois dbinom
logLikehood <- function(Y, MU, Sigma = 1, family)
{
  switch(
    family,
    "gaussian" =   dnorm(Y, MU, Sigma, log = TRUE),
    "poisson"  =   dpois(Y, MU, log = TRUE),
    "binomial" =   dbinom(Y, 1, MU, log = TRUE)
  )
}
