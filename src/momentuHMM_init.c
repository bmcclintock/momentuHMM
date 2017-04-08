#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton (R 3.4)

/* .Call calls */
extern SEXP momentuHMM_dbeta_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dexp_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dgamma_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dlnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dpois_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dweibull_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_dwrpcauchy_rcpp(SEXP, SEXP, SEXP);
extern SEXP momentuHMM_nLogLike_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP momentuHMM_trMatrix_rcpp(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"momentuHMM_dbeta_rcpp",      (DL_FUNC) &momentuHMM_dbeta_rcpp,       3},
  {"momentuHMM_dexp_rcpp",       (DL_FUNC) &momentuHMM_dexp_rcpp,        3},
  {"momentuHMM_dgamma_rcpp",     (DL_FUNC) &momentuHMM_dgamma_rcpp,      3},
  {"momentuHMM_dlnorm_rcpp",     (DL_FUNC) &momentuHMM_dlnorm_rcpp,      3},
  {"momentuHMM_dpois_rcpp",      (DL_FUNC) &momentuHMM_dpois_rcpp,       3},
  {"momentuHMM_dvm_rcpp",        (DL_FUNC) &momentuHMM_dvm_rcpp,         3},
  {"momentuHMM_dweibull_rcpp",   (DL_FUNC) &momentuHMM_dweibull_rcpp,    3},
  {"momentuHMM_dwrpcauchy_rcpp", (DL_FUNC) &momentuHMM_dwrpcauchy_rcpp,  3},
  {"momentuHMM_nLogLike_rcpp",   (DL_FUNC) &momentuHMM_nLogLike_rcpp,   11},
  {"momentuHMM_trMatrix_rcpp",   (DL_FUNC) &momentuHMM_trMatrix_rcpp,    3},
  {NULL, NULL, 0}
};

void R_init_momentuHMM(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
