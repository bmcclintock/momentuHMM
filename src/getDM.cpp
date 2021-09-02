#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
//' Get design matrix
//'
//' Loop for creating full design matrix (X) from pseudo-design matrix (DM). Written in C++. Used in \code{getDM}.
//'
//' @param X full design matrix
//' @param covs matrix of covariates
//' @param DM pseudo design matrix
//' @param nr number of rows in design matrix
//' @param nc number of column in design matrix
//' @param cov covariate names
//' @param nbObs number of observations
//'
//' @return full design matrix (X)
// [[Rcpp::export]]
arma::cube getDM_rcpp(arma::cube X, arma::mat covs, CharacterVector DM, unsigned int nr, unsigned int nc, CharacterVector cov, unsigned int nbObs)
{
  unsigned int covsize = (unsigned int) cov.size();
  for(unsigned int v=0; v<covsize; v++){
    for(unsigned int j=0; j<nr; j++){
      for(unsigned int l=0; l<nc; l++){
        if(DM(l*nr+j)==cov[v]){
          for(unsigned int k=0; k<nbObs; k++){
            X(j,l,k) = covs(k,v);
          }
        }
      }
    }
  }
  return X;
}
