#include <RcppArmadillo.h>
#include <iostream>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Stationary distribution for a continuous-time Markov chain
//'
//' Written in C++. 
//'
//' @param A transition rate matrix (of dimension nbStates x nbStates)
//'
//' @return row vector of stationary distribution probabilities
// [[Rcpp::export]]
arma::rowvec stationary_rcpp(arma::mat A) {
  
  int nbStates = A.n_rows;
  arma::mat Q, R, B;
  
  A = A.t();
  arma::mat v1(1,nbStates);
  v1.ones();
  A.insert_rows(nbStates,v1);
  arma::qr(Q,R,A);
  
  arma::rowvec pi;
  B = Q.t();
  arma::mat b(nbStates+1,1);
  b.zeros();
  b(nbStates,0) = 1;
  pi = arma::solve(R,B*b).t();
  return pi;
}