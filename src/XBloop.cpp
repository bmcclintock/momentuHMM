#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

//' Get XB
//'
//' Loop for computation of design matrix (X) times the working scale parameters (B). Written in C++. Used in \code{\link{w2n}}.
//'
//' @param DM design matrix
//' @param Xvec working parameters
//' @param nbObs number of observations
//' @param nr number of rows in design matrix
//' @param nc number of column in design matrix
//'
//' @return XB matrix
// [[Rcpp::export]]
arma::mat XBloop_rcpp(List DM, NumericVector Xvec, int nbObs, int nr, int nc)
{
  unsigned int i;
  unsigned int j;
  unsigned int k;
  arma::mat XB(nr,nbObs);
  XB.zeros();
  for(i=0; i<nr; i++){
    for(j=0; j<nc; j++){
      NumericVector DMelem = DM[j*nr+i];
      int DMsize = DMelem.size();
      for(k=0; k<nbObs; k++){
        if(DMsize>1){
          XB(i,k) += DMelem[k] * Xvec[j];
        } else {
          XB(i,k) += DMelem[0] * Xvec[j];
        }
      }
    }
  }
  return XB;
}