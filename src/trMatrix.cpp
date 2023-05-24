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
//' @param kappa maximum allowed value for the row sums of the off-diagonal elements in the state transition rate matrix, such that the minimum value for the diagonal elements is \code{-kappa}. Default: \code{Inf}. Setting less than \code{Inf} can help avoid numerical issues during optimization, in which case the transition rate parameters \code{beta} are on the logit scale (instead of the log scale).
//' @param check logical indicating whether or not to check transition probability matrix for issues. Ignored unless \code{CT=TRUE}.
//'
//' @return Three dimensional array \code{trMat}, such that \code{trMat[,,t]} is the transition matrix at
//' time t.
// [[Rcpp::export]]
arma::cube trMatrix_rcpp(int nbStates, arma::mat beta, arma::mat covs, IntegerVector betaRef, bool CT = false, NumericVector dt = NumericVector::create(), IntegerVector aInd = IntegerVector::create(), bool rateMatrix = false, double kappa = NA_REAL, bool check = true)
{
  int nbObs = covs.n_rows;
  bool maxInd;
  if(!R_FINITE(kappa)){
    kappa = R_PosInf;
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
    for(unsigned int k=0;k<(unsigned)nbObs;k++) {
      Gamma.zeros();
      int cpt=0; // counter for diagonal elements
      for(int i=0;i<nbStates;i++) {
        for(int j=0;j<nbStates;j++) {
          if(j==(betaRef(i)-1)) {
            cpt++;
          } else {
            if(i!=j){
              Gamma(i,j) = exp(g(k,i*nbStates+j-cpt));//1.e+5 * R::plogis(g(k,i*nbStates+j-cpt),0.0,1.0,true,false);
              Gamma(i,i) -= Gamma(i,j);
            } else Gamma(i,j) -= exp(g(k,i*nbStates+j-cpt));//-1.e+5 * R::plogis(-g(k,i*nbStates+j-cpt),0.0,1.0,true,false);
          }
          //Rprintf("Gamma %d %d %f Gamma %d %d %f betaRef %d \n",i+1,j+1,Gamma(i,j),i+1,i+1,Gamma(i,i),betaRef(i));
        }
        if(i!=(betaRef(i)-1)){
          for(int l=0;l<nbStates;l++){
            if(l!=(betaRef(i)-1)) Gamma(i,betaRef(i)-1) -= Gamma(i,l);
            //if(i!=l) rowSums(i,k) += Gamma(i,l);
            //Rprintf("Gamma %d %d %f Gamma %d %d %f rowSums %f \n",i+1,l+1,Gamma(i,l),i+1,betaRef(i),Gamma(i,betaRef(i)-1),rowSums(i,k));
          }
        }
        for(int l=0;l<nbStates;l++){
          //if(l!=(betaRef(i)-1)) Gamma(i,betaRef(i)-1) -= Gamma(i,l);
          if(i!=l) rowSums(i,k) += Gamma(i,l);
          //Rprintf("Gamma %d %d %f Gamma %d %d %f rowSums %f \n",i+1,l+1,Gamma(i,l),i+1,betaRef(i),Gamma(i,betaRef(i)-1),rowSums(i,k));
        }
        if(maxInd){
          Gamma(i,i) = 0.;
          for(int l=0;l<nbStates;l++){
            if(i!=l){
              Gamma(i,l) = kappa * Gamma(i,l) / (1.+rowSums(i,k)); // normalization
              Gamma(i,i) -= Gamma(i,l);
              //Rprintf("k %d Gamma %d %d %f Gamma %d %d %f kappa %f rowSums %f \n",k+1,i+1,l+1,Gamma(i,l),i+1,i+1,Gamma(i,i),kappa,rowSums(i,k));
            }
          }
          //if((abs(Gamma(i,i))-kappa)>1.e-8){
          //  for(int l=0;l<nbStates;l++){
          //    Rprintf("k %d Gamma %d %d %f Gamma %d %d %f kappa %f rowSums %f diff %f \n",k+1,i+1,l+1,Gamma(i,l),i+1,i+1,Gamma(i,i),kappa,rowSums(i,k),Gamma(i,i)+kappa);
          //  }
          //}
        }
      }
      if(!rateMatrix){
        if(ani<aInd.size() && k==(unsigned)(aInd(ani)-1)){
          if(nbObs==1) {
            try {
              trMat.slice(k) = expmatrix_rcpp(Gamma * dt(k),kappa,check);
            }
            catch(std::exception &ex) {	
              Rprintf("Observation %d: ",k+1);
              forward_exception_to_r(ex);
            }
          }
          else trMat.slice(k) = expmatrix_rcpp(Gamma * 0,kappa,check);
        } else {
          try {
            trMat.slice(k) = expmatrix_rcpp(Gamma * dt(k-1),kappa,check);
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
