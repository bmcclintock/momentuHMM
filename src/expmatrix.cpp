#include <RcppArmadillo.h>
using namespace Rcpp;

typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;
extern void (*expm)(double *x, int n, double *z, precond_type precond_kind);

//' Matrix Exponential
//'
//' This function computes the exponential of a square matrix \code{A}, defined as the sum from \code{r=0} to infinity of \code{A^r/r!}, using the default method of \code{\link[expm]{expm}}.
//'
//' @param x a square matrix.
//' @param kappa maximum allowed value for the row sums of the off-diagonal elements in the state transition rate matrix, such that the minimum value for the diagonal elements is \code{-kappa}. Default: \code{Inf}. Setting less than \code{Inf} can help avoid numerical issues during optimization.
//' @param check logical indicating whether or not to check transition probability matrix for issues
//'
//' @return The matrix exponential of x.
// [[Rcpp::export]]
arma::mat expmatrix_rcpp(arma::mat x, double kappa = NA_REAL, bool check=false) {
  int nbStates = x.n_rows;
  if(!R_FINITE(kappa)) kappa = R_PosInf;
  for(int i=0; i<nbStates; i++){
    if (!R_FINITE(x(i,i))){
      /*		*err = -1; return; */
      stop("numerical overflow in calculating TPM\n");
    }
    x(i,i) = 0.;
    for (int j=0; j<nbStates; j++) {
      if(i!=j){
        x(i,i) -= x(i,j);
      }
    }
  }
  arma::mat z(nbStates,nbStates);
  expm(x.begin(), nbStates, z.begin(), Ward_2);
  double zrowSum = 0.0;
  for(int i=0; i<nbStates; i++){
    zrowSum = 0.0;
    for (int j=0; j<nbStates; j++) {
      if (!R_FINITE(z(i,j))){
        /*		*err = -1; return; */
        stop("numerical overflow in calculating TPM\n");
      }
      if (z(i,j) < DBL_EPSILON) z(i,j) = 0.;
      if (z(i,j) > 1. - DBL_EPSILON) z(i,j) = 1.;
      zrowSum += z(i,j);
    }
    if(check){
      if(abs(1. - zrowSum) > 1.e-5) {
        warning("Transition probability matrix rows for state %d do not sum to one, most likely due to numerical overflow or underflow. These should not be trusted (check 'kappa'; see output for 'beta' and 'Q' parameters)\n",i+1);
      }
    }
  }
  return z;
}
