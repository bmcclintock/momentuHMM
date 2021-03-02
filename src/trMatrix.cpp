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
//'
//' @return Three dimensional array \code{trMat}, such that \code{trMat[,,t]} is the transition matrix at
//' time t.
// [[Rcpp::export]]
arma::cube trMatrix_rcpp(int nbStates, arma::mat beta, arma::mat covs, IntegerVector betaRef, bool CT = false, NumericVector dt = NumericVector::create())
{
  int nbObs = covs.n_rows;
  arma::cube trMat(nbStates,nbStates,nbObs);
  trMat.zeros();
  arma::mat rowSums(nbStates,nbObs);
  rowSums.zeros();
  
  arma::mat Gamma(nbStates,nbStates);
  
  arma::mat g(nbObs,nbStates*(nbStates-1));
  g = covs*beta;

  if(CT){
    g = exp(g);
    for(int k=0;k<nbObs;k++) {
      Gamma.zeros();
      int cpt=0; // counter for diagonal elements
      for(int i=0;i<nbStates;i++) {
        for(int j=0;j<nbStates;j++) {
          if(j==(betaRef(i)-1)) {
            if(i!=j){
              for(int l=0;l<(nbStates-1);l++){
                Gamma(i,j) +=  g(k,i*nbStates+l-cpt) * dt(k);
              }
            }
            cpt++;
          } else {
            if(i!=j) Gamma(i,j) = g(k,i*nbStates+j-cpt) * dt(k);
          }
        }
        for(int l=0;l<nbStates;l++){
          if(i!=l) Gamma(i,i) -= Gamma(i,l);
        }
      }
      trMat.slice(k) = expmatrix_rcpp(Gamma);
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
