#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton (R 3.6)

/* .Call calls */
extern SEXP _momentuHMM_cbindmean2(SEXP, SEXP);
extern SEXP _momentuHMM_cbindmean3(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_cbindsigma2(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_cbindsigma3(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_combine(SEXP);
extern SEXP _momentuHMM_dbern_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dbeta_rcpp(SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_dcat_rcpp(SEXP, SEXP, SEXP);
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
extern SEXP _momentuHMM_getDM_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_nLogLike_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_trMatrix_rcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _momentuHMM_XBloop_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_momentuHMM_cbindmean2",      (DL_FUNC) &_momentuHMM_cbindmean2,       2},
    {"_momentuHMM_cbindmean3",      (DL_FUNC) &_momentuHMM_cbindmean3,       3},
    {"_momentuHMM_cbindsigma2",     (DL_FUNC) &_momentuHMM_cbindsigma2,      3},
    {"_momentuHMM_cbindsigma3",     (DL_FUNC) &_momentuHMM_cbindsigma3,      6},
    {"_momentuHMM_combine",         (DL_FUNC) &_momentuHMM_combine,          1},
    {"_momentuHMM_dbern_rcpp",      (DL_FUNC) &_momentuHMM_dbern_rcpp,       3},
    {"_momentuHMM_dbeta_rcpp",      (DL_FUNC) &_momentuHMM_dbeta_rcpp,       3},
    {"_momentuHMM_dcat_rcpp",       (DL_FUNC) &_momentuHMM_dcat_rcpp,        3},
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
    {"_momentuHMM_getDM_rcpp",      (DL_FUNC) &_momentuHMM_getDM_rcpp,       7},
    {"_momentuHMM_nLogLike_rcpp",   (DL_FUNC) &_momentuHMM_nLogLike_rcpp,   13},
    {"_momentuHMM_trMatrix_rcpp",   (DL_FUNC) &_momentuHMM_trMatrix_rcpp,    4},
    {"_momentuHMM_XBloop_rcpp",     (DL_FUNC) &_momentuHMM_XBloop_rcpp,     11},
    {NULL, NULL, 0}
};

void R_init_momentuHMM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
