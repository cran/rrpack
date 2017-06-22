##Soft thresholding
##lambda can be a matrix with same dim as C
softTH <- function(C, Lam, tol = 1e-4) {
  Ct <- abs(C) - Lam
  Ct[Ct < tol] <- 0
  sign(C) * Ct
}



##Hard thresholding for rows
##lambda can be a vector
softrowTH <- function(C, Lam, tol = 1e-4) {
  l2row <-
    sqrt(rowSums(C ^ 2)) # apply(C, 1, function(x) sqrt(sum(x^2)))
  under <- l2row - Lam < tol
  C[under,] <- 0
  C[!under,] <- C[!under,] * (1 - Lam / l2row)[!under]
  list(C = C, l2 = l2row)
}

## soft thresholding function
softThres <- function(x, lambda) {
  xt <- abs(x) - lambda
  xt[xt <= 0] <- 0
  sign(x) * xt
}


## penalty function for soft thresholding
softPen <- function(x, lambda) {
  lambda * sum(abs(x))
}


## hard thresholding function
hardThres <- function(x, lambda) {
  ##  ifelse(abs(x) > lambda, x, 0) ## nice but too slow
  xt <- abs(x) - lambda
  x[xt <= 0] <- 0
  x
}

## penalty for hard thresholding
hardPen <- function(x, lambda) {
  absx <- abs(x)
  0.5 * sum(lambda ^ 2 - (absx - lambda) ^ 2 * (absx < lambda))
}

## Cadinality constraint
countThres <- function(x, lambda) {
  ## lambda is the cardinality constraint
  if (lambda == 0) {
    x[TRUE] <- 0
  } else{
    xt <- abs(x)
    cut <- xt[order(x, decreasing = TRUE)[lambda]]
    x[xt < cut] <- 0
  }
  x
}

####Since it is constrained estiamtion,
####it does not contribute to the objective function
countPen <- function(x, lambda) {
  0
}

####Need to include SCAD and MCP
# scadTH<-function(C, lambda, a=3.7)
# {
#   absC <- abs(C)
#   signC <- sign(C)
#   alambda <- a*lambda
#   C*(absC>alambda) +
#     (((a-1)*C-(signC*alambda))/(a-2))*(absC<=alambda & absC>2*lambda) +
#     (signC*max(0,(absC-lambda)))*(absC<=2*lambda)
#
# }


## singular value thresholding
svdThres <-
  function(C,
           maxrank = min(dim(C)),
           threshold = hardThres,
           lambda = 0) {
    svdC <- svd(C, nu = maxrank, nv = maxrank)
    d <- threshold(svdC$d, lambda)
    rank <- sum(d != 0)
    list(
      C = svdC$u %*% diag(d, nrow = maxrank, ncol = maxrank) %*% t(svdC$v),
      rank = rank
    )
  }

## singular value penalty
svdPen <- function(C,
                   maxrank = min(dim(C)),
                   penalty = softPen,
                   lambda = 0) {
  svdC <- svd(C, nu = maxrank, nv = maxrank)
  penalty(svdC$d, lambda)
}



## ##Hard thresholding for rows
## ##lambda can be a vector
hardrowTH <- function(C, lambda, eta = 0) {
  l2row <- apply(C, 1, function(x)
    sqrt(sum(x ^ 2)))
  C <- C / (1 + eta)
  C[l2row <= lambda, ] <- 0
  list(C = C, l2 = l2row)
}

## ###Rowwise penalties and thresholding rules
## hardrowTH <- function(C, lambda) {
##   l2row <- apply(C,1,function(x)sqrt(sum(x^2)))
##   ##C <- C/(1+eta)
##   C[l2row<=lambda,] <- 0
##   list(C=C,l2=l2row)
## }

## hardrowPEN <- function(C,lambda) {
##   l2row <- apply(C,1,function(x)sqrt(sum(x^2)))
##   lambda^2/2*sum(l2row>0)
## }


## softrowTH <- function(C,lambda) {
##   l2row <- apply(C,1,function(x)sqrt(sum(x^2)))
##   ##C <- C/(1+eta)
##   l2row <- softTH(l2row,lambda)
##   C <- diag(l2row,nrow=nrow(C),ncol=nrow(C))%*%C
##   list(C=C,l2=l2row)
## }

## softrowPEN <- function(C,lambda) {
##   l2row <- apply(C,1,function(x)sqrt(sum(x^2)))
##   lambda*sum(l2row)
## }

## countrowTH <- function(C,count) {
##   ##l2row <- apply(C,1,function(x)sqrt(sum(x^2)))
##   C <- C/(1+eta)
##   cut <- sort(l2row,decreasing = TRUE)[count+1]
##   C[l2row<=cut,] <- 0
##   l2row[l2row<=cut] <- 0
##   list(C=C,l2=l2row)
## }

## countrowPEN <- function(C,count) {
##   0
## }







## ##Mean function for glm
## ##kappa values are scaling factors
## ##X0: design matrix
## ##offset: offset matrix
## ##C0: coefficient matrix
## ##S0: mean-shift matrix
## ##kappa_C0, kappa_S0: scaling factor to ensure convergence
## mu.gaussian <- function(X0, offset,
##                        kappa_C0, C0,
##                        kappa_S0 = 1, S0 = matrix(nrow=nrow(X0),ncol=ncol(C0),0)){

##   offset + X0%*%C0/kappa_C0 + S0/kappa_S0

## }


## ##negative log-likelihood
## obj.gaussian <- function(Y, X0, offest,
##                         kappa_C0 = 1, C0,
##                         kappa_S0 = 1, S0 = matrix(nrow=nrow(Y),ncol=ncol(Y),0)){

##   MU <- mu.gaussian(X0,offset,kappa_C0,C0,kappa_S0,S0)
##   -sum(dnorm(Y,MU,log=TRUE))

## }

## ##Mean function for glm-poisson
## mu.poisson <- function(X0, offset,
##                        kappa_C0, C0,
##                        kappa_S0 = 1, S0 = matrix(nrow=nrow(X0),ncol=ncol(C0),0)){

##   exp(offset + X0%*%C0/kappa_C0 + S0/kappa_S0)

## }

## obj.poisson <- function(Y, X0, offest,
##                         kappa_C0 = 1, C0,
##                         kappa_S0 = 1, S0 = matrix(nrow=nrow(Y),ncol=ncol(Y),0)){

##   MU <- mu.poisson(X0,offset,kappa_C0,C0,kappa_S0,S0)
##   -sum(dpois(Y,MU,log=TRUE))

## }

## ##Use kappa values to scale
## mu.logistic <- function(X0, offset,
##                         kappa_C0 = 1, C0,
##                         kappa_S0 = 1, S0=matrix(nrow=nrow(X0),ncol=ncol(C0),0)){
##   1/(1+exp(-offset - X0%*%C0/kappa_C0 - S0/kappa_S0))
##   ##1/(1+exp(- X0%*%C0/kappa_C0 - S0/kappa_S0))
## }


## obj.logistic <- function(Y,X0, offset,
##                          kappa_C0=1, C0,
##                          kappa_S0=1, S0 = matrix(nrow=nrow(Y),ncol=ncol(Y),0)){

##   MU <- mu.logistic(X0,offset,kappa_C0,C0,kappa_S0,S0)
##   -sum(dbinom(Y,1,MU,log=TRUE))

## }
