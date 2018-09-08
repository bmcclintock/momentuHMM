#ifndef _DENSITIES_
#define _DENSITIES_

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace Rcpp;
using namespace std;

//' Gamma density function
//'
//' Probability density function of the gamma distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean
//' @param sigma Standard deviation
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dgamma_rcpp(NumericVector x, arma::mat mu, arma::mat sigma)
{
  arma::colvec res(x.size());
  double shape;
  double scale;
  
  for(int i=0;i<x.size();i++) {
    //if(!arma::is_finite(x(i)))
    //  res(i) = 1; // if missing observation
    //else {
      // convert mean and sd to shape and scale
      shape = pow(mu(i),2)/pow(sigma(i),2);
      scale = pow(sigma(i),2)/mu(i);
      
      res(i) = R::dgamma(x(i),shape,scale,0);
    //}
  }
  
  return res;
}

//' Weibull density function
//'
//' Probability density function of the Weibull distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param shape Shape
//' @param scale Scale
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dweibull_rcpp(NumericVector x, arma::mat shape, arma::mat scale)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dweibull(x(i),shape(i),scale(i),0);
  }
  
  return res;
}

//' Normal density function
//'
//' Probability density function of the normal distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param mean Mean of the distribution 
//' @param sd Standard deviation of the distribution 
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dnorm_rcpp(NumericVector x, arma::mat mean, arma::mat sd)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dnorm(x(i),mean(i),sd(i),0);
  }
  
  return res;
}

//' Log-normal density function
//'
//' Probability density function of the log-normal distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param meanlog Mean of the distribution on the log-scale
//' @param sdlog Standard deviation of the distribution on the log-scale
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dlnorm_rcpp(NumericVector x, arma::mat meanlog, arma::mat sdlog)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dlnorm(x(i),meanlog(i),sdlog(i),0);
  }
  
  return res;
}

//' Exponential density function
//'
//' Probability density function of the exponential distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param rate Rate
//' @param foo Unused (for compatibility with template)
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dexp_rcpp(NumericVector x, arma::mat rate, arma::mat foo)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dexp(x(i),1./rate(i),0); // R::dexp expects scale=1/rate
  }
  
  return res;
}

//' Von Mises density function
//'
//' Probability density function of the Von Mises distribution, defined as a function
//' of the modified Bessel function of order 0 (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean
//' @param kappa Concentration
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dvm_rcpp(NumericVector x, arma::mat mu, arma::mat kappa)
{
  arma::colvec res(x.size());
  double b;
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // is missing observation
    else {
      b = R::bessel_i(kappa(i),0,2);
      res(i) = 1/(2*M_PI*b)*pow((exp(cos(x(i)-mu(i))-1)),kappa(i));
    }
  }
  
  return res;
}

//' Wrapped Cauchy density function
//'
//' Probability density function of the wrapped Cauchy distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean
//' @param rho Concentration
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dwrpcauchy_rcpp(NumericVector x, arma::mat mu, arma::mat rho)
{
  arma::colvec res(x.size());
  double num;
  double den;
  
  for(int i=0;i<x.size();i++) {
    //if(!arma::is_finite(x(i)))
    //  res(i) = 1; // if missing observation
    //else {
      num = 1-rho(i)*rho(i);
      den = (2*M_PI)*(1+rho(i)*rho(i)-2*rho(i)*cos(x(i)-mu(i)));
      res(i) = num/den;
    //}
  }
  
  return res;
}

//' Probability density function of the beta distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param shape1 Shape1
//' @param shape2 Shape2
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dbeta_rcpp(NumericVector x, arma::mat shape1, arma::mat shape2)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dbeta(x(i),shape1(i),shape2(i),0);
  }
  
  return res;
}

//' Poisson density function
//'
//' Probability density function of the Poisson distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param rate Rate
//' @param foo Unused (for compatibility with template)
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dpois_rcpp(NumericVector x, arma::mat rate, arma::mat foo)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dpois(x(i),rate(i),0);
  }
  
  return res;
}

//' Bernoulli density function
//'
//' Probability density function of the Bernoulli distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param prob success probability
//' @param foo Unused (for compatibility with template)
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dbern_rcpp(NumericVector x, arma::mat prob, arma::mat foo)
{
  arma::colvec res(x.size());
  
  for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dbinom(x(i),1,prob(i),0);
  }
  
  return res;
}

const double log2pi2 = log(2.0 * M_PI)/2.0;

//' C++ implementation of multivariate Normal probability density function for multiple inputs
//'
//'@param x data matrix of dimension \code{p x n}, \code{p} being the dimension of the
//'data and n the number of data points.
//'@param mean mean vectors matrix of dimension \code{p x n}
//'@param varcovM list of length \code{n} of variance-covariance matrices,
//'each of dimensions \code{p x p}.
//'@param Log logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}.
//'
//'@return matrix of densities of dimension \code{K x n}.
// [[Rcpp::export]]
arma::colvec dmvnorm_rcpp(NumericVector x,
                       const arma::mat mean,
                       const arma::mat varcovM){
  
  //const bool & Log=true;
  int p = mean.n_rows;
  int n = mean.n_cols;
  arma::vec xvec(p);
  arma::vec out(n);
  arma::mat sigma(p,p);
  double constant = - p*log2pi2;
  
  for (int i=0; i < n; i++) {
    for(int k=0; k < p; k++){
      sigma(k,k) = varcovM(k*p+k,i);
      for(int j=0; j < k; j++){
        sigma(k,j) = varcovM(k*p+j,i);
        sigma(j,k) = varcovM(k*p+j,i);
      }
      xvec(k) = x(k*n+i);
    }
    arma::mat Rinv = inv(trimatu(chol(sigma)));
    //mat R = chol(as<arma::mat>(varcovM[i]));
    double logSqrtDetvarcovM = sum(log(Rinv.diag()));
    arma::colvec mtemp = mean.col(i);
    arma::colvec x_i = xvec - mtemp;
    arma::rowvec xRinv = trans(x_i)*Rinv;
    //vec xRinv = solve(trimatl(R.t()), x_i);
    double quadform = sum(xRinv%xRinv);
    //if (!Log) {
      out(i) = exp(-0.5*quadform + logSqrtDetvarcovM + constant);
    //} else{
    //  out(i) = -0.5*quadform + logSqrtDetvarcovM + constant;
    //}
    //Rprintf("i %d x %f y %f mean %f %f sigma %f %f %f %f out %f \n",i,xvec(0),xvec(1),mtemp(0),mtemp(1),sigma(0,0),sigma(1,0),sigma(0,1),sigma(1,1),out(i));
  }
  
  return out;
  
}

// used in nLogLike_rcpp to map the functions' names to the functions
typedef arma::colvec (*FunPtr)(NumericVector, arma::mat, arma::mat);

#endif
