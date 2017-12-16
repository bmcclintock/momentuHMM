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
arma::mat XBloop_rcpp(List DM, NumericVector Xvec, unsigned int nbObs, unsigned int nr, unsigned int nc, bool circularAngleMean, bool consensus, IntegerVector rindex, arma::mat cindex, int nbStates)
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
  
  unsigned int nrindex = rindex.size();
  unsigned int incr = (circularAngleMean) ? 2 : 1;
  
  NumericVector DMelem2;
  unsigned int Xcount;
  double Xsum;
  
  for(i=0; i<nrindex; i++){
    Xcount = 0;
    Xsum = 1.;
    for(j=0; j<nc; j+=incr){
      //Rprintf("i %d j %d Xcount %d \n",i,j,Xcount);
      if(cindex(rindex[i],j)){
        NumericVector DMelem = DM[j*nr+rindex[i]];
        if(circularAngleMean) {
          DMelem2 = DM[(j+1)*nr+rindex[i]];
          Xsum += abs(Xvec[Xcount]);
        }
        int DMsize = DMelem.size();
        for(k=0; k<nbObs; k++){
          if(DMsize>1){
            if(!circularAngleMean){
              XB(rindex[i],k) += DMelem[k] * Xvec[j];
            } else {
              XB1(rindex[i],k) += DMelem[k]*Xvec[Xcount];
              XB2(rindex[i],k) += DMelem2[k]*Xvec[Xcount];
            }
          } else {
            if(!circularAngleMean){
              XB(rindex[i],k) += DMelem[0] * Xvec[j];
            } else {
              XB1(rindex[i],k) += DMelem[0]*Xvec[Xcount];
              XB2(rindex[i],k) += DMelem2[0]*Xvec[Xcount];
              //Rprintf("i %d j %d rindex %d cindex %f Xcount %d DMelem %f DMelem2 %f Xvec %f XB1 %f XB2 %f \n",i,j,rindex[i],cindex(rindex[i],j),Xcount,DMelem[0],DMelem2[0],Xvec[Xcount],XB1(rindex[i],k),XB2(rindex[i],k));
            }
          }
        }
      }
      Xcount += 1;
    }
    if(circularAngleMean){
      XB.row(rindex[i]) = atan2(XB1.row(rindex[i]),1.+XB2.row(rindex[i]));
      if(consensus) XB.row(rindex[i]+nbStates) = sqrt(pow(XB1.row(rindex[i]),2.)+pow(1.+XB2.row(rindex[i]),2.)) / Xsum;
    }
  }
  return XB;
}