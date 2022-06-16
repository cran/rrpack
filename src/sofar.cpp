#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace sugar;

// [[Rcpp::export]]
Rcpp::List procrustes_RCpp (arma::mat XY, arma::mat XX, arma::mat D, double rho2,
                             arma::mat U,
                             Rcpp::List control){

  int maxit = control["maxit"];
  double epsilon = control["epsilon"];
  int niter = 1;
  double diff = 10*epsilon;
  arma::mat Z;
  arma::mat u;
  arma::mat v;
  arma::vec s;
  Rcpp::List out;

  arma::mat XYD = XY*D;
  arma::mat Dsq = D*D;
  arma::mat U0 = U;

  XX = -XX;
  XX.diag() += rho2;

  while ((diff > epsilon) & (niter < maxit)) {
    U0 = U;
    Z = XX*U0*Dsq + XYD;
    svd(u, s, v, Z);  //qXq qXr rXr
    u = u.cols(0, v.n_cols-1);
    U = u*v.t();
    niter +=1;
    //diff = norm(U-U0)/norm(U0);
    diff = norm(U-U0);
  }


  out["U"]  = U;
  out["diff"]  = diff;
  out["niter"]  = niter;
  return(out);

}
