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
//' @param x quantile
//' @param mu Mean
//' @param sigma Standard deviation
//'
//' @return density
// [[Rcpp::export]]
double dgamma_rcpp(double x, double mu, double sigma)
{
  //arma::colvec res(x.size());
  double res;
  
  // convert mean and sd to shape and scale
  double shape = pow(mu,2)/pow(sigma,2);
  double scale = pow(sigma,2)/mu;
  
  //for(int i=0;i<x.size();i++) {
  if(!arma::is_finite(x))
    res = 1; // if missing observation
  else
    res = R::dgamma(x,shape,scale,0);
  //}
  
  return res;
}

//' Weibull density function
//'
//' Probability density function of the Weibull distribution (written in C++)
//'
//' @param x quantile
//' @param shape Shape
//' @param scale Scale
//'
//' @return density
// [[Rcpp::export]]
double dweibull_rcpp(double x, double scale, double shape)
{
  //arma::colvec res(x.size());
  double res;
  
  //for(int i=0;i<x.size();i++) {
  if(!arma::is_finite(x))
    res = 1; // if missing observation
  else
    res = R::dweibull(x,shape,scale,0);
  //}
  
  return res;
}

//' Log-normal density function
//'
//' Probability density function of the log-normal distribution (written in C++)
//'
//' @param x quantile
//' @param meanlog Mean of the distribution on the log-scale
//' @param sdlog Standard deviation of the distribution on the log-scale
//'
//' @return density
// [[Rcpp::export]]
double dlnorm_rcpp(double x, double meanlog, double sdlog)
{
  //arma::colvec res(x.size());
  double res;
  
  //for(int i=0;i<x.size();i++) {
  if(!arma::is_finite(x))
    res = 1; // if missing observation
  else
    res = R::dlnorm(x,meanlog,sdlog,0);
  //}
  
  return res;
}

//' Exponential density function
//'
//' Probability density function of the exponential distribution (written in C++)
//'
//' @param x quantile
//' @param rate Rate
//' @param foo Unused (for compatibility with template)
//'
//' @return density
// [[Rcpp::export]]
double dexp_rcpp(double x, double rate, double foo=0)
{
  //arma::colvec res(x.size());
  double res;
  
  //for(int i=0;i<x.size();i++) {
  if(!arma::is_finite(x))
    res = 1; // if missing observation
  else
    res = R::dexp(x,1/rate,0); // R::dexp expects scale=1/rate
  //}
  
  return res;
}

//' Von Mises density function
//'
//' Probability density function of the Von Mises distribution, defined as a function
//' of the modified Bessel function of order 0 (written in C++)
//'
//' @param x quantile
//' @param mu Mean
//' @param kappa Concentration
//'
//' @return density
// [[Rcpp::export]]
double dvm_rcpp(double x, double mu, double kappa)
{
  //arma::colvec res(x.size());
  double res;
  double b = R::bessel_i(kappa,0,2);
  
  //for(int i=0;i<x.size();i++) {
  if(!arma::is_finite(x))
    res = 1; // is missing observation
  else
    res = 1/(2*M_PI*b)*pow((exp(cos(x-mu)-1)),kappa);
  //}
  
  return res;
}

//' Wrapped Cauchy density function
//'
//' Probability density function of the wrapped Cauchy distribution (written in C++)
//'
//' @param x quantile
//' @param mu Mean
//' @param rho Concentration
//'
//' @return density
// [[Rcpp::export]]
double dwrpcauchy_rcpp(double x, double mu, double rho)
{
  //arma::colvec res(x.size());
  double res;
  double num = 1-rho*rho;
  double den;
  
  //for(int i=0;i<x.size();i++) {
  if(!arma::is_finite(x))
    res = 1; // if missing observation
  else {
    den = (2*M_PI)*(1+rho*rho-2*rho*cos(x-mu));
    res = num/den;
  }
  //}
  
  return res;
}

//' Probability density function of the beta distribution (written in C++)
//'
//' @param x quantile
//' @param shape1 Shape1
//' @param shape2 Shape2
//'
//' @return density
// [[Rcpp::export]]
double dbeta_rcpp(double x, double shape1, double shape2)
{
  //arma::colvec res(x.size());
  double res;

  //for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x))
      res = 1; // if missing observation
    else
      res = R::dbeta(x,shape1,shape2,0);
  //}

  return res;
}

//' Poisson density function
//'
//' Probability density function of the Poisson distribution (written in C++)
//'
//' @param x quantile
//' @param rate Rate
//' @param foo Unused (for compatibility with template)
//'
//' @return density
// [[Rcpp::export]]
double dpois_rcpp(double x, double rate, double foo=0)
{
  //arma::colvec res(x.size());
  double res;

  //for(int i=0;i<x.size();i++) {
    if(!arma::is_finite(x))
      res = 1; // if missing observation
    else
      res = R::dpois(x,rate,0);
  //}

  return res;
}

// used in nLogLike_rcpp to map the functions' names to the functions
typedef double (*FunPtr)(double,double,double);
//typedef arma::colvec (*FunPtr2)(NumericVector,double,double);

#endif
