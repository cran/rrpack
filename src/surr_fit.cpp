#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
using namespace arma;
using namespace Rcpp;
using namespace std;

arma::vec penreg_Rcpp(arma::vec Y,
          arma::mat X,
                      double lambda,
                      arma::vec beta0,
                      List control);

arma::vec penreg_Rcpp_XY(arma::vec XY,
          arma::mat XX,
                      double lambda,
                      arma::vec beta0,
                      List control);


//using namespace std;
//using namespace sugar;

//
//double softThres(double x, double lambda) {
//    return((x > lambda) ? x - lambda :
//     (x < -lambda) ? x + lambda : 0.);
//}

//
//// Penalized regression with lasso penalty
////
//// @param Y response vector
//// @param X covariate matrix
//// @param lambda scalar penalty parameter
//// @param beta0 initial value of regression coefficient vector
//// @param control a list of parameters controling the fitting process
//// @return estimated regression coefficient vector
//// [[Rcpp::export]]
//arma::vec penreg_Rcpp(arma::vec Y,
//            arma::mat X,
//                    double lambda,
//                    arma::vec beta0,
//                    List control) {
//    int p = X.n_cols, n = X.n_rows, i, j;
//    arma::vec vlambda(p), ss(p), bchg(p), r(n);
//
//    int maxit = control["maxit"];
//    double epsilon = control["epsilon"], bj;
//    //bool trace = control["trace"];
//    ss = arma::conv_to< arma::vec >::from(sqrt(sum(X % X)));
//
//    for (i = 0; i < p; i++) X.col(i) /= ss(i);
//    vlambda = lambda * pow(ss, -1.);
//
//    beta0 = beta0 % ss;
//    r = Y - X * beta0;
//    arma::vec beta(beta0);
//
//    for (i = 0; i < maxit; i++) {
//	for (j = 0; j < p; j++) {
//          bj = sum(X.col(j) % r) + beta(j);
//          bj = softThres(bj, vlambda(j));
//         // if (control["trace"]) {
//		// trace here
//         // }
//          beta(j) = bj;
//          bchg(j) = beta(j) - beta0(j);
//          if (bchg(j) != 0) r -= bchg(j) * X.col(j);
//	}
//	if (sqrt(sum(bchg % bchg) / sum(beta0 % beta0)) < epsilon) break;
//	beta0 = beta;
//    }
//    beta /=  ss;
//    return(beta);
//}

//
// // Sparse unit-rank regression, Rcpp version
// // @param Y response matrix
// // @param X covariate matrix
// // @param lambda scalar penalty parameter
// // @param U0 initial value of U
// // @param V0 initial value of V
// // @param WU weight vector for U
// // @param WV weight vector for V
// // @param Xtran weight transformed X; if NULL X\%*\% WU is used
// // @param control a list of parameters controlling the fitting process
// // @return S3 \code{surr} object as `surr.fit'
// // [[Rcpp::export]]
// Rcpp::List surr_fit_Rcpp_old(arma::mat Y, arma::mat X, double lambda,
//         arma::vec U0, arma::vec V0, //double D0,
//         arma::vec WU, arma::vec WV,
//         arma::mat Xtran, Rcpp::List control){
//
//   int p = X.n_cols;
//   int q = Y.n_cols;
//   int n = Y.n_rows;
//   arma::vec y = vectorise(Y);
//   Rcpp::List out;
//
//   int maxit = control["maxit"];
//   double epsilon = control["epsilon"];
//   Rcpp::List innerControl;
//   innerControl["epsilon"] = control["innerEpsilon"];
//   innerControl["maxit"] = control["innerMaxit"];
//   //innerControl["trace"] = TRUE;
//
//   //matrices used in iterations
//   arma::vec U = U0;
//   arma::vec V = V0;
//   double D = 0.0;
//   //double D = D0;
//
//   arma::vec Up;
//   arma::vec Vp;
//   //double Dp;
//
//   arma::mat Xu;
//   arma::mat Xv;
//   //quantites used in iterations
//   double norm_UV;
//   double conv = 0.0;
//   double sse, BIC, BICP, AIC;
//   arma::vec res;
//   bool converged;
//   int dfu, dfv;
//
//   //Initial values should be passed in already.
//   //U, V, WU, WV
//   //if (Xtran.is_empty()) {
//   //  Xtran = X*diagmat(WU);
//   //}
//   U = U/WU;
//   V = V/WV;
//   U = U/norm(U,2);
//   V = V/norm(V,2);
//
//
//   bool flag = FALSE;
//   for(int j = 1; j < maxit; j++){
//
//     Up = U;
//     Vp = V;
//     //Dp = D;
//     norm_UV = norm(U*V.t()*as_scalar(D),2);
//
//     Xu = kron(WV%V,Xtran);
//     U = penreg_Rcpp(y,Xu,lambda,Up,innerControl);
//     D = norm(U,1);
//     if (D <  std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
//       D = 0.0;
//       flag = TRUE;
//       conv = epsilon;
//       conv *= 2;
//       break;
//     }
//     U = U/as_scalar(D);
//
//     Xv = kron(diagmat(WV),Xtran*U);
//     V = penreg_Rcpp(y,Xv,lambda,Vp,innerControl);
//     D = norm(V,1);
//     if (D <  std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
//       D = 0.0;
//       flag = TRUE;
//       conv = epsilon;
//       conv *= 2;
//       break;
//     }
//     V = V/as_scalar(D);
//
//     conv = norm(U*V.t()-Up*Vp.t(),2)/norm_UV;
//     if(conv <= epsilon) break;
//   }
//
//   if(flag == TRUE){
//     res = y;
//     sse = pow(norm(y,2),2);
//     V = V.zeros();
//     U = U.zeros();
//     D = 0.0;
//     dfu = 0.0;
//     dfv = 0.0;
//     BIC = log(sse);
//     BICP = log(sse);
//     AIC = log(sse);
//   }else{
//     res = y - Xv*V*as_scalar(D);
//     sse = pow(norm(res,2),2);
//     dfu =  accu(U != 0);
//     dfv =  accu(V != 0);
//     V = WV%V*as_scalar(D);
//     U = WU%U;
//     double dv = norm(V,2);
//     double du = norm(U,2);
//     D = dv*du;
//     U = U/du;
//     V = V/dv;
//     BIC = log(sse) + log(q*n)/q/n*(dfu+dfv-1);
//     BICP = log(sse) + 2*log(p*q)/q/n*(dfu+dfv-1);
//     AIC = log(sse) + 2/q/n*(dfu+dfv-1);
//   }
//
//   if(conv <= epsilon){
//     converged = TRUE;
//   }else{
//     converged = FALSE;
//   }
//   arma::vec df(2);
//   df(0) = dfu;
//   df(1) = dfv;
//   arma::vec ic(3);
//   ic(0) = BIC;
//   ic(1) = BICP;
//   ic(2) = AIC;
//
//   out["sse"] = sse;
//   out["df"] = df;
//   out["ic"] = ic;
//   out["U"] = U;
//   out["V"] = V;
//   out["D"] = D;
//   out["converged"] = converged;
//   return(out);
// }
//
//



// Sparse unit-rank regression, Rcpp version
// @param Y response matrix
// @param X covariate matrix
// @param lambda scalar penalty parameter
// @param U0 initial value of U
// @param V0 initial value of V
// @param WU weight vector for U
// @param WV weight vector for V
// @param Xtran weight transformed X; if NULL X\%*\% WU is used
// @param control a list of parameters controlling the fitting process
// @return S3 \code{surr} object as `surr.fit'
// [[Rcpp::export]]
Rcpp::List surr_fit_Rcpp(arma::mat Y, arma::mat X, double lambda,
        arma::vec U0, arma::vec V0, //double D0,
        arma::vec WU, arma::vec WV,
        arma::mat Xtran, Rcpp::List control){

  int p = X.n_cols;
  int q = Y.n_cols;
  int n = Y.n_rows;
  arma::vec y = vectorise(Y);
  Rcpp::List out;

  arma::mat XXtran = Xtran.t()*Xtran;
  arma::mat XYtran = Xtran.t()*Y;
  //arma::mat YXtran = Y.t()*Xtran;

  int maxit = control["maxit"];
  double epsilon = control["epsilon"];
  Rcpp::List innerControl;
  innerControl["epsilon"] = control["innerEpsilon"];
  innerControl["maxit"] = control["innerMaxit"];
  //innerControl["trace"] = TRUE;

  //matrices used in iterations
  arma::vec U = U0;
  arma::vec V = V0;
  double D = 0.0;
  //double D = D0;

  arma::vec Up;
  arma::vec Vp;
  //double Dp;

  arma::mat Xu;
  arma::mat Xv;
  //quantites used in iterations
  double norm_UV;
  double conv = 0.0;
  double sse, BIC, BICP, AIC;
  arma::vec res;
  bool converged;
  double dfu, dfv;

  //Initial values should be passed in already.
  //U, V, WU, WV
  //if (Xtran.is_empty()) {
  //  Xtran = X*diagmat(WU);
  //}

  U = U/WU;
  V = V/WV;
  U = U/norm(U,1);
  V = V/norm(V,1);

  bool flag = FALSE;
  for(int j = 1; j < maxit; j++){

    Up = U;
    Vp = V;
    //Dp = D;
    norm_UV = norm(U*V.t()*as_scalar(D),2);

    //Xu = kron(WV%V,Xtran);
    //U = penreg_Rcpp(y,Xu,lambda,Up,innerControl);
    //XXu = XXtran*norm(WV%V)^2;
    //XYu = XYtran%WV%V;
    U = penreg_Rcpp_XY(XYtran*(WV%V),XXtran*pow(norm(WV%V,2),2),lambda,Up,innerControl);

    D = norm(U,1);
    if (D <  std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
      D = 0.0;
      flag = TRUE;
      conv = epsilon;
      conv *= 2;
      break;
    }
    U = U/as_scalar(D);

    Xv = kron(diagmat(WV),Xtran*U);
    //V = penreg_Rcpp(y,Xv,lambda,Vp,innerControl);
    //XXv = diagmat(WV%WV)*sum(pow(Xtran*U,2));
    //XYv = diagmat(WV)*YXtran*U;
    V = penreg_Rcpp_XY(WV%(XYtran.t()*U),diagmat(WV%WV*as_scalar(pow(norm(Xtran*U,2),2))),lambda,Vp,innerControl);
    //V = penreg_Rcpp_XY(WV%(XYtran.t()*U),Xv.t()*Xv,lambda,Vp,innerControl);

    D = norm(V,1);
    if (D <  std::numeric_limits<double>::epsilon() && D > -std::numeric_limits<double>::epsilon()){
      D = 0.0;
      flag = TRUE;
      conv = epsilon;
      conv *= 2;
      break;
    }
    V = V/as_scalar(D);

    conv = norm(U*V.t()-Up*Vp.t(),2)/norm_UV;
    if(conv <= epsilon) break;
  }

  if(flag == TRUE){
    res = y;
    sse = pow(norm(y,2),2);
    V = V.zeros();
    U = U.zeros();
    D = 0.0;
    dfu = 0.0;
    dfv = 0.0;
    BIC = std::log(sse);
    BICP = std::log(sse);
    AIC = std::log(sse);
  }else{
    res = y - Xv*V*as_scalar(D);
    sse = pow(norm(res,2),2);
    dfu =  accu(U != 0);
    dfv =  accu(V != 0);
    V = WV%V*as_scalar(D);
    U = WU%U;
    double dv = norm(V,2);
    double du = norm(U,2);
    D = dv*du;
    U = U/du;
    V = V/dv;
    BIC = std::log(sse) + std::log(static_cast <double>(q*n))/q/n*(dfu+dfv-1);
    BICP = std::log(sse) + 2*std::log(static_cast <double>(p*q))/q/n*(dfu+dfv-1);
    AIC = std::log(sse) + 2/q/n*(dfu+dfv-1);
  }

  if(conv <= epsilon){
    converged = TRUE;
  }else{
    converged = FALSE;
  }
  arma::vec df(2);
  df(0) = dfu;
  df(1) = dfv;
  arma::vec ic(3);
  ic(0) = BIC;
  ic(1) = BICP;
  ic(2) = AIC;

  out["sse"] = sse;
  out["df"] = df;
  out["ic"] = ic;
  out["U"] = U;
  out["V"] = V;
  out["D"] = D;
  out["converged"] = converged;
  return(out);
}
