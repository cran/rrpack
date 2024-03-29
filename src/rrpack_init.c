#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _rrpack_kron_RcppArma(SEXP, SEXP);
extern SEXP _rrpack_lasso_shooting(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rrpack_MGlasso_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rrpack_penreg_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rrpack_penreg_Rcpp_XY(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rrpack_procrustes_RCpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rrpack_srrr_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _rrpack_surr_fit_Rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rssvd_orth(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_rrpack_kron_RcppArma",   (DL_FUNC) &_rrpack_kron_RcppArma,    2},
    {"_rrpack_lasso_shooting",  (DL_FUNC) &_rrpack_lasso_shooting,   5},
    {"_rrpack_MGlasso_Rcpp",    (DL_FUNC) &_rrpack_MGlasso_Rcpp,     6},
    {"_rrpack_penreg_Rcpp",     (DL_FUNC) &_rrpack_penreg_Rcpp,      5},
    {"_rrpack_penreg_Rcpp_XY",  (DL_FUNC) &_rrpack_penreg_Rcpp_XY,   5},
    {"_rrpack_procrustes_RCpp", (DL_FUNC) &_rrpack_procrustes_RCpp,  6},
    {"_rrpack_srrr_Rcpp",       (DL_FUNC) &_rrpack_srrr_Rcpp,       12},
    {"_rrpack_surr_fit_Rcpp",   (DL_FUNC) &_rrpack_surr_fit_Rcpp,    9},
    {"rssvd_orth",              (DL_FUNC) &rssvd_orth,               8},
    {NULL, NULL, 0}
};

void R_init_rrpack(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
