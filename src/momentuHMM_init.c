#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "expm.h"

// Mostly generated with tools::package_native_routine_registration_skeleton(getwd(),,,FALSE) (R 4.0)
// additionally requires: 
// #include "expm.h" in header
// expm = (void (*) (double*, int, double*, precond_type)) R_GetCCallable("expm", "expm"); in R_init_momentuHMM

/* .Call calls */
extern SEXP _momentuHMM_cbindmean2(SEXP, SEXP);
extern SEXP _momentuHMM_cbindmean3(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_cbindsigma2(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_cbindsigma3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_combine(SEXP, SEXP);
extern SEXP _momentuHMM_dbern_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dbeta_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dcat_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dcrwrice_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dcrwvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dexp_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dgamma_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dlnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dlogis_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dmvnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dnbinom_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dnorm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dpois_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dt_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dvm_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dweibull_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dwrpcauchy_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_expmatrix_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_getDM_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_nLogLike_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_stationary_rcpp(SEXP);
extern SEXP _momentuHMM_trMatrix_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_XBloop_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_momentuHMM_cbindmean2",      (DL_FUNC) &_momentuHMM_cbindmean2,       2},
    {"_momentuHMM_cbindmean3",      (DL_FUNC) &_momentuHMM_cbindmean3,       3},
    {"_momentuHMM_cbindsigma2",     (DL_FUNC) &_momentuHMM_cbindsigma2,      3},
    {"_momentuHMM_cbindsigma3",     (DL_FUNC) &_momentuHMM_cbindsigma3,      6},
    {"_momentuHMM_combine",         (DL_FUNC) &_momentuHMM_combine,          2},
    {"_momentuHMM_dbern_rcpp",      (DL_FUNC) &_momentuHMM_dbern_rcpp,       3},
    {"_momentuHMM_dbeta_rcpp",      (DL_FUNC) &_momentuHMM_dbeta_rcpp,       3},
    {"_momentuHMM_dcat_rcpp",       (DL_FUNC) &_momentuHMM_dcat_rcpp,        3},
    {"_momentuHMM_dcrwrice_rcpp",   (DL_FUNC) &_momentuHMM_dcrwrice_rcpp,    3},
    {"_momentuHMM_dcrwvm_rcpp",     (DL_FUNC) &_momentuHMM_dcrwvm_rcpp,      3},
    {"_momentuHMM_dexp_rcpp",       (DL_FUNC) &_momentuHMM_dexp_rcpp,        3},
    {"_momentuHMM_dgamma_rcpp",     (DL_FUNC) &_momentuHMM_dgamma_rcpp,      3},
    {"_momentuHMM_dlnorm_rcpp",     (DL_FUNC) &_momentuHMM_dlnorm_rcpp,      3},
    {"_momentuHMM_dlogis_rcpp",     (DL_FUNC) &_momentuHMM_dlogis_rcpp,      3},
    {"_momentuHMM_dmvnorm_rcpp",    (DL_FUNC) &_momentuHMM_dmvnorm_rcpp,     3},
    {"_momentuHMM_dnbinom_rcpp",    (DL_FUNC) &_momentuHMM_dnbinom_rcpp,     3},
    {"_momentuHMM_dnorm_rcpp",      (DL_FUNC) &_momentuHMM_dnorm_rcpp,       3},
    {"_momentuHMM_dpois_rcpp",      (DL_FUNC) &_momentuHMM_dpois_rcpp,       3},
    {"_momentuHMM_dt_rcpp",         (DL_FUNC) &_momentuHMM_dt_rcpp,          3},
    {"_momentuHMM_dvm_rcpp",        (DL_FUNC) &_momentuHMM_dvm_rcpp,         3},
    {"_momentuHMM_dweibull_rcpp",   (DL_FUNC) &_momentuHMM_dweibull_rcpp,    3},
    {"_momentuHMM_dwrpcauchy_rcpp", (DL_FUNC) &_momentuHMM_dwrpcauchy_rcpp,  3},
    {"_momentuHMM_expmatrix_rcpp",  (DL_FUNC) &_momentuHMM_expmatrix_rcpp,   3},
    {"_momentuHMM_getDM_rcpp",      (DL_FUNC) &_momentuHMM_getDM_rcpp,       7},
    {"_momentuHMM_nLogLike_rcpp",   (DL_FUNC) &_momentuHMM_nLogLike_rcpp,   16},
    {"_momentuHMM_stationary_rcpp", (DL_FUNC) &_momentuHMM_stationary_rcpp,  1},
    {"_momentuHMM_trMatrix_rcpp",   (DL_FUNC) &_momentuHMM_trMatrix_rcpp,   10},
    {"_momentuHMM_XBloop_rcpp",     (DL_FUNC) &_momentuHMM_XBloop_rcpp,     11},
    {NULL, NULL, 0}
};

void R_init_momentuHMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    expm = (void (*) (double*, int, double*, precond_type)) R_GetCCallable("expm", "expm");
}
