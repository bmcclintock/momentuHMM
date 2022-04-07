#include <RcppArmadillo.h>
using namespace Rcpp;

typedef enum {Ward_2, Ward_1, Ward_buggy_octave} precond_type;
extern void (*expm)(double *x, int n, double *z, precond_type precond_kind);

//' Matrix Exponential
//'
//' This function computes the exponential of a square matrix \code{A}, defined as the sum from \code{r=0} to infinity of \code{A^r/r!}, using the default method of \code{\link[expm]{expm}}.
//'
//' @param x a square matrix.
//'
//' @return The matrix exponential of x.
// [[Rcpp::export]]
arma::mat expmatrix_rcpp(arma::mat x) {
  int nbStates = x.n_rows;
  for(int i=0; i<nbStates; i++){
    if (!R_FINITE(x(i,i))){
      /*		*err = -1; return; */
      stop("numerical overflow in calculating TPM\n");
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
      zrowSum += z(i,j);
      if (z(i,j) < DBL_EPSILON) z(i,j) = 0;
      if (z(i,j) > 1 - DBL_EPSILON) z(i,j) = 1;
    }
    if(abs(1 - zrowSum) > 1.e-8) stop("numerical overflow in calculating TPM\n");
  }
  return z;
}
