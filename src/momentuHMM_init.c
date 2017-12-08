#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton (R 3.4)

/* .Call calls */
extern SEXP _momentuHMM_dbern_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dbeta_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dexp_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dgamma_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dlnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dpois_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dweibull_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dwrpcauchy_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_getDM_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_nLogLike_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_trMatrix_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_XBloop_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_momentuHMM_dbern_rcpp",      (DL_FUNC) &_momentuHMM_dbern_rcpp,       3},
    {"_momentuHMM_dbeta_rcpp",      (DL_FUNC) &_momentuHMM_dbeta_rcpp,       3},
    {"_momentuHMM_dexp_rcpp",       (DL_FUNC) &_momentuHMM_dexp_rcpp,        3},
    {"_momentuHMM_dgamma_rcpp",     (DL_FUNC) &_momentuHMM_dgamma_rcpp,      3},
    {"_momentuHMM_dlnorm_rcpp",     (DL_FUNC) &_momentuHMM_dlnorm_rcpp,      3},
    {"_momentuHMM_dnorm_rcpp",      (DL_FUNC) &_momentuHMM_dnorm_rcpp,       3},
    {"_momentuHMM_dpois_rcpp",      (DL_FUNC) &_momentuHMM_dpois_rcpp,       3},
    {"_momentuHMM_dvm_rcpp",        (DL_FUNC) &_momentuHMM_dvm_rcpp,         3},
    {"_momentuHMM_dweibull_rcpp",   (DL_FUNC) &_momentuHMM_dweibull_rcpp,    3},
    {"_momentuHMM_dwrpcauchy_rcpp", (DL_FUNC) &_momentuHMM_dwrpcauchy_rcpp,  3},
    {"_momentuHMM_getDM_rcpp",      (DL_FUNC) &_momentuHMM_getDM_rcpp,       7},
    {"_momentuHMM_nLogLike_rcpp",   (DL_FUNC) &_momentuHMM_nLogLike_rcpp,   11},
    {"_momentuHMM_trMatrix_rcpp",   (DL_FUNC) &_momentuHMM_trMatrix_rcpp,    3},
    {"_momentuHMM_XBloop_rcpp",     (DL_FUNC) &_momentuHMM_XBloop_rcpp,      8},
    {NULL, NULL, 0}
};

void R_init_momentuHMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}