#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;
using namespace sugar;

// [[Rcpp::export]]
Rcpp::List MGlasso_C (arma::mat Y, arma::mat X, arma::vec lam, arma::mat B0, double conv, int maxiter) {
  // min |Y-XB|^2 + lam*|B|
  int  p=X.n_cols, iter=0, j;  // n=Y.n_rows, q=Y.n_cols,
  double  diff=10*conv, l2B1, sse;
  arma::rowvec sh;
  //arma::mat mat1=eye(p,p);
  arma::mat B1;
  arma::mat res1;
  arma::mat res1j;
  arma::mat XRj;
  Rcpp::List out;
  if (lam.size() == 1) {lam = as_scalar(lam)*ones(p);}
  sh = sum(square(X), 0);
  //if (B0.is_finite()) {
    B1 = B0;
  //} else {
  //  Rcpp::List ini = rrr_cpp(Y, X, 1, FALSE, mat1, TRUE, TRUE);
  //  arma::mat iniC = ini["C_ls"];
  //  B1 = iniC;
  //}
  res1 = Y - X * B1;
  while ((diff > conv) & (iter < maxiter)) {
    B0 = B1;
    for (j = 0; j < p; j++) {
      res1j = res1 +  X.col(j)* B1.row(j); //n q
      XRj =   trans(X.col(j)) * res1j;    //1 q
      arma::rowvec t1=XRj/as_scalar(sh(j))*max(0.0,1-lam(j)/sqrt(accu(square(XRj))));
      B1.row(j) = t1;
      res1 = res1j - X.col(j)* B1.row(j);
    }
    l2B1 = accu(square(B1));
    if (l2B1 == 0) {
      iter = maxiter;
    } else {
      diff = sqrt(accu(square(B0 - B1))/l2B1);
      iter = iter + 1;
    }
  }
  sse = accu(square(Y - X * B0));
  out["B"] = B1;
  out["sse"] = sse;
  out["iter"] = iter;
  return(out);
}

// Row-sparse Reduced-Rank Regression
//
// @param Y response matrix
// @param X covariate matrix
// @param method method
// @param A0 initial value
// @param V0 initial value
// @param nrank rank
// @param lambda tuning parameter
// @param conv conv
// @param maxiter maxiter
// @param inner_conv inner conv
// @param inner_iter inner iter
// @param WA weights
// @return estimation results
// [[Rcpp::export]]
Rcpp::List srrr_Rcpp (arma::mat Y, arma::mat X, String method, arma::mat A0, arma::mat V0,
                   int nrank, double lambda,
                   double conv, int maxiter,
                   double inner_conv, int inner_iter,
                   arma::vec WA) {
// min |Y-XAV|^2 + lamA*|A|, s.t. V'V = I_nrank

  int n = Y.n_rows, p = X.n_cols, q = Y.n_cols;
  arma::mat YX = Y.t() * X;
  //vec lamA = as_scalar(lambda)*ones(p);

  int xrank = accu(svd(X) > 0.01);
  int iter=0, dfu0, dfv0;
  bool conv_flag;
  double l2C1, sse, df, BIC, BICP, AIC, GCV, GIC;
  arma::vec s;
  arma::vec diff(maxiter+1);
  diff.fill(2*conv);
  arma::mat mat1 = eye(q, q);   ///, iniU, iniD, iniV, iniC;
  arma::mat V1;
  arma::mat A1;
  arma::mat C1;
  arma::mat C0;
  arma::mat W;
  arma::mat u;
  arma::mat v;
  arma::mat residual;

  Rcpp::List ini;
  Rcpp::List out_MGlasso;
  Rcpp::List out;

  //ini = rrr_cpp(Y, X, nrank, FALSE, mat1, TRUE, TRUE);
  //arma::mat iniU = ini["U"];
  //arma::vec inid = ini["D"];
  //arma::mat iniV = ini["V"];
  //arma::mat iniC = ini["C"];
  //arma::mat iniD = diagmat(inid);

  //if (V0.is_empty() || A0.is_empty()) {
  //  V1 = iniV;  //q r
  //  A1 = iniU * iniD;  //p r
  //  C1 = iniC; //p q
  //} else {
    V1 = V0;
    A1 = A0;
    C1 = A0*V0.t();
  //}

//  if (WA.is_empty()) {
//    if((method=="glasso")) {
//        WA = ones(p);
//    } else if((method=="adglasso")) {
//      A1 = iniU * iniD;  //p r
//      //vec A1norm(p);
//      for(int i=0;i<p;i++){
//        WA(i) = pow(sqrt(accu(square(A1.row(i)))),-wgamma);
//      }
//      //WA = pow(A1norm, -wgamma);
//    }
//  }

  while ((iter < maxiter) & (diff(iter) > conv)) {
  //while (iter < 10) {
    V0 = V1;
    A0 = A1;
    C0 = C1;

    arma::mat YV0 = Y*V0;
    out_MGlasso = MGlasso_C(YV0, X, lambda*WA, A0, inner_conv, inner_iter);
    arma::mat MGlassoB = out_MGlasso["B"];
    A1 = MGlassoB;   //p r
    W = YX * A1;      //q r
    svd(u, s, v, W);  //qXq qXr rXr
    u = u.cols(0, nrank-1);
    V1 = u * v.t();
    C1 = A1*V1.t();
    l2C1 = accu(square(C1));
    if (l2C1 == 0) {
      diff(iter) = 0;
    } else {
      iter = iter + 1;
      diff(iter) = sqrt(accu(square(C0 - C1))/l2C1);
    }

  }

  diff = diff.subvec(0, iter);
  residual = Y - X * C1;
  sse = accu(square(residual));
  dfu0 = accu(A1 != 0);
  dfv0 = accu(V1 != 0);
  df = dfu0 * xrank/p + dfv0 - nrank*nrank;
  double logqn = log(q * n);
  double logsse = log(sse);
  BIC = logsse + df*logqn/(q*n);
  BICP = logsse + 2*df*logqn/(q*n);
  AIC = logsse + 2*df/(q*n);
  GIC = logsse + df*log(logqn)*log(p*q)/(q*n);
  double dfqn2 = pow((1 - df/(q*n)), 2);
  GCV = sse/(q*n*dfqn2);
  if (diff(iter) <= conv) {
    conv_flag = TRUE;
  } else {
    conv_flag = FALSE;
  }

  out["diff"] = diff;
  out["iter"] = iter;
  out["BIC"]  = BIC;
  out["BICP"] = BICP;
  out["AIC"]  = AIC;
  out["GCV"]  = GCV;
  out["GIC"]  = GIC;
  out["sse"]  = sse;
  out["df"]   = df;
  out["conv_flag"] = conv_flag;
  out["A"]  = A1;
  out["V"]  = V1;
  out["C"]  = C1;
  return(out);
}
