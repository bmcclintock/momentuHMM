// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// dgamma_rcpp
arma::colvec dgamma_rcpp(NumericVector x, double mu, double sigma);
RcppExport SEXP momentuHMM_dgamma_rcpp(SEXP xSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(dgamma_rcpp(x, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// dweibull_rcpp
arma::colvec dweibull_rcpp(NumericVector x, double shape, double scale);
RcppExport SEXP momentuHMM_dweibull_rcpp(SEXP xSEXP, SEXP shapeSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(dweibull_rcpp(x, shape, scale));
    return rcpp_result_gen;
END_RCPP
}
// dlnorm_rcpp
arma::colvec dlnorm_rcpp(NumericVector x, double meanlog, double sdlog);
RcppExport SEXP momentuHMM_dlnorm_rcpp(SEXP xSEXP, SEXP meanlogSEXP, SEXP sdlogSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type meanlog(meanlogSEXP);
    Rcpp::traits::input_parameter< double >::type sdlog(sdlogSEXP);
    rcpp_result_gen = Rcpp::wrap(dlnorm_rcpp(x, meanlog, sdlog));
    return rcpp_result_gen;
END_RCPP
}
// dexp_rcpp
arma::colvec dexp_rcpp(NumericVector x, double rate, double foo);
RcppExport SEXP momentuHMM_dexp_rcpp(SEXP xSEXP, SEXP rateSEXP, SEXP fooSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type foo(fooSEXP);
    rcpp_result_gen = Rcpp::wrap(dexp_rcpp(x, rate, foo));
    return rcpp_result_gen;
END_RCPP
}
// dvm_rcpp
arma::colvec dvm_rcpp(NumericVector x, double mu, double kappa);
RcppExport SEXP momentuHMM_dvm_rcpp(SEXP xSEXP, SEXP muSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(dvm_rcpp(x, mu, kappa));
    return rcpp_result_gen;
END_RCPP
}
// dwrpcauchy_rcpp
arma::colvec dwrpcauchy_rcpp(NumericVector x, double mu, double rho);
RcppExport SEXP momentuHMM_dwrpcauchy_rcpp(SEXP xSEXP, SEXP muSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(dwrpcauchy_rcpp(x, mu, rho));
    return rcpp_result_gen;
END_RCPP
}
// dbeta_rcpp
arma::colvec dbeta_rcpp(NumericVector x, double shape1, double shape2);
RcppExport SEXP momentuHMM_dbeta_rcpp(SEXP xSEXP, SEXP shape1SEXP, SEXP shape2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type shape1(shape1SEXP);
    Rcpp::traits::input_parameter< double >::type shape2(shape2SEXP);
    rcpp_result_gen = Rcpp::wrap(dbeta_rcpp(x, shape1, shape2));
    return rcpp_result_gen;
END_RCPP
}
// dpois_rcpp
arma::colvec dpois_rcpp(NumericVector x, double rate, double foo);
RcppExport SEXP momentuHMM_dpois_rcpp(SEXP xSEXP, SEXP rateSEXP, SEXP fooSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< double >::type foo(fooSEXP);
    rcpp_result_gen = Rcpp::wrap(dpois_rcpp(x, rate, foo));
    return rcpp_result_gen;
END_RCPP
}
// nLogLike_rcpp
double nLogLike_rcpp(int nbStates, arma::mat beta, arma::mat covs, DataFrame data, std::string stepDist, std::string angleDist, std::string omegaDist, std::string dryDist, std::string diveDist, std::string iceDist, std::string landDist, arma::mat stepPar, arma::mat anglePar, arma::mat omegaPar, arma::mat dryPar, arma::mat divePar, arma::mat icePar, arma::mat landPar, arma::rowvec delta, IntegerVector aInd, bool zeroInflation, bool stationary);
RcppExport SEXP momentuHMM_nLogLike_rcpp(SEXP nbStatesSEXP, SEXP betaSEXP, SEXP covsSEXP, SEXP dataSEXP, SEXP stepDistSEXP, SEXP angleDistSEXP, SEXP omegaDistSEXP, SEXP dryDistSEXP, SEXP diveDistSEXP, SEXP iceDistSEXP, SEXP landDistSEXP, SEXP stepParSEXP, SEXP angleParSEXP, SEXP omegaParSEXP, SEXP dryParSEXP, SEXP diveParSEXP, SEXP iceParSEXP, SEXP landParSEXP, SEXP deltaSEXP, SEXP aIndSEXP, SEXP zeroInflationSEXP, SEXP stationarySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nbStates(nbStatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< std::string >::type stepDist(stepDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type angleDist(angleDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type omegaDist(omegaDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type dryDist(dryDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type diveDist(diveDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type iceDist(iceDistSEXP);
    Rcpp::traits::input_parameter< std::string >::type landDist(landDistSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type stepPar(stepParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type anglePar(angleParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type omegaPar(omegaParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type dryPar(dryParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type divePar(diveParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type icePar(iceParSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type landPar(landParSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type aInd(aIndSEXP);
    Rcpp::traits::input_parameter< bool >::type zeroInflation(zeroInflationSEXP);
    Rcpp::traits::input_parameter< bool >::type stationary(stationarySEXP);
    rcpp_result_gen = Rcpp::wrap(nLogLike_rcpp(nbStates, beta, covs, data, stepDist, angleDist, omegaDist, dryDist, diveDist, iceDist, landDist, stepPar, anglePar, omegaPar, dryPar, divePar, icePar, landPar, delta, aInd, zeroInflation, stationary));
    return rcpp_result_gen;
END_RCPP
}
// trMatrix_rcpp
arma::cube trMatrix_rcpp(int nbStates, arma::mat beta, arma::mat covs);
RcppExport SEXP momentuHMM_trMatrix_rcpp(SEXP nbStatesSEXP, SEXP betaSEXP, SEXP covsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nbStates(nbStatesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type covs(covsSEXP);
    rcpp_result_gen = Rcpp::wrap(trMatrix_rcpp(nbStates, beta, covs));
    return rcpp_result_gen;
END_RCPP
}
