#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
arma::cube getDM_rcpp(arma::cube DM, NumericVector covs, arma::cube tmpDM, int nr, int nc, std::string cov, int nbObs)
{
  for(unsigned int k=0; k<nbObs; k++){
    for(unsigned int j=0; j<nr; j++){
      for(unsigned int l=0; l<nc; l++){
        if(tmpDM(j,l,k)) DM(j,l,k) = covs(k);
      }
    }
  }
  return DM;
}
