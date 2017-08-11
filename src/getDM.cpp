#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
// [[Rcpp::export]]
arma::cube getDM_rcpp(arma::cube DM, arma::mat covs, CharacterVector tmpDM, int nr, int nc, CharacterVector cov, int nbObs)
{
  int covsize = cov.size();
  for(unsigned int v=0; v<covsize; v++){
    for(unsigned int j=0; j<nr; j++){
      for(unsigned int l=0; l<nc; l++){
        if(tmpDM(l*nr+j)==cov[v]){
          for(unsigned int k=0; k<nbObs; k++){
            DM(j,l,k) = covs(k,v);
          }
        }
      }
    }
  }
  return DM;
}
