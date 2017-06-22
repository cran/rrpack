#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Kronecker product, Rcpp version
// @param A first matrix
// @param B second matrix
// [[Rcpp::export]]
arma::mat kron_RcppArma(arma::mat A, arma::mat B) {
    return(arma::kron(A, B));
}


