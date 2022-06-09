#include <RcppArmadillo.h>
#include <iostream>
#include "expmatrix.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;

//' Transition probability matrix
//'
//' Computation of the transition probability matrix, as a function of the covariates and the regression
//' parameters. Written in C++. Used in \code{\link{viterbi}}.
//'
//' @param nbStates Number of states
//' @param beta Matrix of regression parameters
//' @param covs Matrix of covariate values
//' @param betaRef Indices of reference elements for t.p.m. multinomial logit link.
//' @param CT logical indicating discrete-time approximation of a continuous-time model
//' @param dt numeric vector of length \code{nrow(covs)} indicating the time difference between observations. Ignored unless \code{CT=TRUE}.
//' @param aInd Vector of indices of the rows at which the data (i.e. \code{covs}) switches to another animal. Ignored unless \code{CT=TRUE}.
//' @param rateMatrix logical indicating whether or not to return the transition rate matrix. Ignored unless \code{CT=TRUE}.
//' @param maxRate maximum allowable value for the off-diagonal state transition rate parameters. Default: \code{Inf}. Setting less than \code{Inf} can help avoid numerical issues during optimization.
//' @param check logical indicating whether or not to check transition probability matrix for issues. Ignored unless \code{CT=TRUE}.
//'
//' @return Three dimensional array \code{trMat}, such that \code{trMat[,,t]} is the transition matrix at
//' time t.
// [[Rcpp::export]]
arma::cube trMatrix_rcpp(int nbStates, arma::mat beta, arma::mat covs, IntegerVector betaRef, bool CT = false, NumericVector dt = NumericVector::create(), IntegerVector aInd = IntegerVector::create(), bool rateMatrix = false, double maxRate = NA_REAL, bool check = true)
{
  int nbObs = covs.n_rows;
  bool maxInd;
  if(!R_FINITE(maxRate)){
    maxRate = R_PosInf;
    maxInd = false;
  } else maxInd = true;
  arma::cube trMat(nbStates,nbStates,nbObs);
  trMat.zeros();
  arma::mat rowSums(nbStates,nbObs);
  rowSums.zeros();
  
  arma::mat Gamma(nbStates,nbStates);
  
  arma::mat g(nbObs,nbStates*(nbStates-1));
  g = covs*beta;

  if(CT){
    
    unsigned int ani=0; // animal index
    //g = exp(g);
    for(int k=0;k<nbObs;k++) {
      Gamma.zeros();
      int cpt=0; // counter for diagonal elements
      for(int i=0;i<nbStates;i++) {
        for(int j=0;j<nbStates;j++) {
          if(j==(betaRef(i)-1)) {
            cpt++;
          } else {
            if(maxInd){
              if(i!=j) Gamma(i,j) = maxRate * R::plogis(g(k,i*nbStates+j-cpt),0.0,1.0,true,false);
              else Gamma(i,j) = -maxRate * R::plogis(g(k,i*nbStates+j-cpt),0.0,1.0,true,false);
            } else {
              if(i!=j) Gamma(i,j) = exp(g(k,i*nbStates+j-cpt));//1.e+5 * R::plogis(g(k,i*nbStates+j-cpt),0.0,1.0,true,false);
              else Gamma(i,j) = -exp(g(k,i*nbStates+j-cpt));//-1.e+5 * R::plogis(-g(k,i*nbStates+j-cpt),0.0,1.0,true,false);
            }
          }
        }
        for(int l=0;l<nbStates;l++){
          if((betaRef(i)-1)!=l) Gamma(i,(betaRef(i)-1)) -= Gamma(i,l);
        }
      }
      if(!rateMatrix){
        if(ani<aInd.size() && k==(unsigned)(aInd(ani)-1)){
          if(nbObs==1) {
            try {
              trMat.slice(k) = expmatrix_rcpp(Gamma * dt(k),maxRate,check);
            }
            catch(std::exception &ex) {	
              Rprintf("Observation %d: ",k+1);
              forward_exception_to_r(ex);
            }
          }
          else trMat.slice(k) = expmatrix_rcpp(Gamma * 0,maxRate,check);
        } else {
          try {
            trMat.slice(k) = expmatrix_rcpp(Gamma * dt(k-1),maxRate,check);
          }
          catch(std::exception &ex) {	
            Rprintf("Observation %d: ",k+1);
            forward_exception_to_r(ex);
          }

        }
      } else trMat.slice(k) = Gamma;
      
      if((ani+1<aInd.size() && k==(unsigned)(aInd(ani+1)-2))) ani++;
      //for(int i=0;i<nbStates;i++) {
      //  for(int j=0;j<nbStates;j++) {
      //    trMat(i,j,k) = Gamma(i,j);
      //  }
      //}
    }  
  } else {
    for(int k=0;k<nbObs;k++) {
      int cpt=0;
      for(int i=0;i<nbStates;i++) {
        for(int j=0;j<nbStates;j++) {
          if(j==(betaRef(i)-1)) {
            trMat(i,j,k)=1;
            cpt++;
          }
          else trMat(i,j,k) = exp(g(k,i*nbStates+j-cpt));
          rowSums(i,k)=rowSums(i,k)+trMat(i,j,k);
        }
      }
    }
  
    // normalization
    for(int k=0;k<nbObs;k++)
      for(int i=0;i<nbStates;i++)
        for(int j=0;j<nbStates;j++)
          trMat(i,j,k) = trMat(i,j,k)/rowSums(i,k);
  }
  return trMat;
}
