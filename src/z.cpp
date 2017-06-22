#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


// [[Rcpp::export()]]
arma::mat a5 (arma::mat x, arma::vec y) {
int R = x.n_rows ;
int C = x.n_cols ;
 int p = y.n_elem; 
Rcout << "Rows: " << R << std::endl ;
Rcout << "Cols: " << C << std::endl ;
 Rcout <<"p:    " << p << std::endl ;
return(x) ;
}

// [[Rcpp::export()]]
arma::vec a6 (arma::mat x, arma::vec y) {
int R = x.n_rows ;
int C = x.n_cols ;
 int p = y.n_elem; 
Rcout << "Rows: " << R << std::endl ;
Rcout << "Cols: " << C << std::endl ;
 Rcout <<"p:    " << p << std::endl ;
 arma::vec v = arma::conv_to< arma::vec >::from(sum(x));
 return(v / y) ;
}

// [[Rcpp::export()]]
arma::mat a17(arma::mat x) {
return( exp(x) ) ;
}

// [[Rcpp::export()]]
arma::mat a18(arma::mat x) {
    return( sum(x % x));
}

// [[Rcpp::export()]]
arma::mat a19(arma::mat x) {
    x /= x;
    return(x);
}

// [[Rcpp::export()]]
int a20(List x) {
    return(x["maxit"]);
}
