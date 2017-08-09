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
//' @param circularAngleMean indicator for whether or not circular-circular regression model
//' @param rindex row index for design matrix
//' @param cindex column index for design matrix
//'
//' @return XB matrix
// [[Rcpp::export]]
arma::mat XBloop_rcpp(List DM, NumericVector Xvec, int nbObs, int nr, int nc, bool circularAngleMean, IntegerVector rindex, arma::mat cindex)
{
  unsigned int i;
  unsigned int j;
  unsigned int k;
  arma::mat XB(nr,nbObs);
  XB.zeros();
  arma::mat XB1(nr,nbObs);
  XB1.zeros();
  arma::mat XB2(nr,nbObs);
  XB2.zeros();
  for(i=0; i<rindex.size(); i++){
    for(j=0; j<nc; j++){
      if(cindex(rindex[i],j)){
        NumericVector DMelem = DM[j*nr+rindex[i]];
        int DMsize = DMelem.size();
        for(k=0; k<nbObs; k++){
          if(DMsize>1){
            if(!circularAngleMean){
              XB(rindex[i],k) += DMelem[k] * Xvec[j];
            } else {
              XB1(rindex[i],k) += sin(DMelem[k])*Xvec[j];
              XB2(rindex[i],k) += cos(DMelem[k])*Xvec[j];
            }
          } else {
            if(!circularAngleMean){
              XB(rindex[i],k) += DMelem[0] * Xvec[j];
            } else {
              XB1(rindex[i],k) += sin(DMelem[0])*Xvec[j];
              XB2(rindex[i],k) += cos(DMelem[0])*Xvec[j];
            }
          }
        }
      }
    }
  }
  if(circularAngleMean){
    XB = atan2(XB1,1.+XB2);
  }
  return XB;
}