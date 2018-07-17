#include <RcppArmadillo.h>
#include <iostream>
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
//'
//' @return Three dimensional array \code{trMat}, such that \code{trMat[,,t]} is the transition matrix at
//' time t.
// [[Rcpp::export]]
arma::cube trMatrix_rcpp(int nbStates, arma::mat beta, arma::mat covs, IntegerVector betaRef)
{
  int nbObs = covs.n_rows;
  arma::cube trMat(nbStates,nbStates,nbObs);
  trMat.zeros();
  arma::mat rowSums(nbStates,nbObs);
  rowSums.zeros();

  arma::mat g(nbObs,nbStates*(nbStates-1));
  g = covs*beta;

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

  return trMat;
}
