#include <Rcpp.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector combine(const List& list, double NAvalue)
{
  std::size_t n = list.size();
  
  // Figure out the length of the output vector
  std::size_t total_length = 0;
  for (std::size_t i = 0; i < n; ++i)
    total_length += Rf_length(list[i]);
  
  // Allocate the vector
  NumericVector output = no_init(total_length);
  
  // Loop and fill
  std::size_t index = 0;
  for (std::size_t i = 0; i < n; ++i)
  {
    NumericVector el = list[i];
    std::copy(el.begin(), el.end(), output.begin() + index);
    
    // Update the index
    index += el.size();
  }
  
  for(int i=0;i<output.size();i++) {
    if(!arma::is_finite(output(i))) {
      output(i) = NAvalue;
    }
  }
  
  return output;
  
}

// [[Rcpp::export]]
arma::mat cbindmean2(arma::mat x, arma::mat y){
  NumericMatrix out(2,x.size());
  out(0,_) = as<NumericVector>(wrap(x));
  out(1,_) = as<NumericVector>(wrap(y));
  return as<arma::mat>(out);
}

// [[Rcpp::export]]
arma::mat cbindmean3(arma::mat x, arma::mat y, arma::mat z){
  NumericMatrix out(3,x.size());
  out(0,_) = as<NumericVector>(wrap(x));
  out(1,_) = as<NumericVector>(wrap(y));
  out(2,_) = as<NumericVector>(wrap(z));  
  return as<arma::mat>(out);
}

// [[Rcpp::export]]
arma::mat cbindsigma2(arma::mat x, arma::mat y, arma::mat xy){
  NumericMatrix out(3,x.size());
  out(0,_) = as<NumericVector>(wrap(x));
  out(1,_) = as<NumericVector>(wrap(y));
  out(2,_) = as<NumericVector>(wrap(xy));
  return as<arma::mat>(out);
}

// [[Rcpp::export]]
arma::mat cbindsigma3(arma::mat x, arma::mat y, arma::mat z, arma::mat xy, arma::mat xz, arma::mat yz){
  NumericMatrix out(6,x.size());
  out(0,_) = as<NumericVector>(wrap(x));
  out(1,_) = as<NumericVector>(wrap(y));
  out(2,_) = as<NumericVector>(wrap(z));
  out(3,_) = as<NumericVector>(wrap(xy));
  out(4,_) = as<NumericVector>(wrap(xz));
  out(5,_) = as<NumericVector>(wrap(yz));
  return as<arma::mat>(out);
}

// // [[Rcpp::export]]
//arma::mat concat(arma::mat x, arma::mat y){
//  NumericMatrix out(1,x.size() + y.size());
//  std::copy(x.begin(),x.end(), out.begin());
//  std::copy(y.begin(),y.end(), out.begin() + x.size());
//  return as<arma::mat>(out);
//}