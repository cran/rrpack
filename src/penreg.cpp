#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

double softThres(double x, double lambda) {
    return((x > lambda) ? x - lambda :
           (x < -lambda) ? x + lambda : 0.);
}

// Penalized regression with lasso penalty
//
// @param Y response vector
// @param X covariate matrix
// @param lambda scalar penalty parameter
// @param beta0 initial value of regression coefficient vector
// @param control a list of parameters controling the fitting process
// @return estimated regression coefficient vector
// [[Rcpp::export]]
arma::vec penreg_Rcpp(arma::vec Y,
                      arma::mat X,
                      double lambda,
                      arma::vec beta0,
                      List control) {
    int p = X.n_cols, n = X.n_rows, i, j;
    arma::vec vlambda(p), ss(p), bchg(p), r(n);

    int maxit = control["maxit"];
    double epsilon = control["epsilon"], bj;
    //bool trace = control["trace"];
    ss = arma::conv_to< arma::vec >::from(sqrt(sum(X % X)));

    for (i = 0; i < p; i++) X.col(i) /= ss(i);
    vlambda = lambda * pow(ss, -1.);

    beta0 = beta0 % ss;
    r = Y - X * beta0;
    arma::vec beta(beta0);

    for (i = 0; i < maxit; i++) {
        for (j = 0; j < p; j++) {
            bj = sum(X.col(j) % r) + beta(j);
            bj = softThres(bj, vlambda(j));
           // if (control["trace"]) {
                // trace here
           // }
            beta(j) = bj;
            bchg(j) = beta(j) - beta0(j);
            if (bchg(j) != 0) r -= bchg(j) * X.col(j);
        }
        if (sqrt(sum(bchg % bchg) / sum(beta0 % beta0)) < epsilon) break;
        beta0 = beta;
    }
    beta /=  ss;
    return(beta);
}



// Penalized regression with lasso penalty
//
// @param XY cross product of X and Y
// @param XX Gram matrix
// @param lambda scalar penalty parameter
// @param beta0 initial value of regression coefficient vector
// @param control a list of parameters controling the fitting process
// @return estimated regression coefficient vector
// [[Rcpp::export]]
arma::vec penreg_Rcpp_XY(arma::vec XY,
              arma::mat XX,
                      double lambda,
                      arma::vec beta0,
                      List control) {
    int p = XX.n_cols, i, j;
    arma::vec vlambda(p), bchg(p);
    vlambda = lambda / diagvec(XX);

    int maxit = control["maxit"];
    double epsilon = control["epsilon"], bj;

    //bool trace = control["trace"];
    //ss = arma::conv_to< arma::vec >::from(sqrt(sum(X % X)));
    //for (i = 0; i < p; i++) X.col(i) /= ss(i);
    //vlambda = lambda * pow(ss, -1.);

    //beta0 = beta0 % ss;
    //r = Y - X * beta0;
    arma::vec beta(beta0);

  for (i = 0; i < maxit; i++) {
        for (j = 0; j < p; j++) {

            //bj = sum(X.col(j) % r) + beta(j);
      bj = (XY(j)-sum(XX.col(j) % beta))/XX(j,j) + beta(j);
            //bj = softThres(bj, vlambda(j));
            bj = softThres(bj,vlambda(j));
     // if (control["trace"]) {
                // trace here
           // }
            beta(j) = bj;
            bchg(j) = beta(j) - beta0(j);
            //if (bchg(j) != 0) r -= bchg(j) * X.col(j);
        }
        if (sqrt(sum(bchg % bchg) / sum(beta0 % beta0)) < epsilon) break;
        beta0 = beta;
    }
    //beta /=  ss;
    return(beta);
}
