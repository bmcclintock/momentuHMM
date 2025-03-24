#ifndef _DENSITIES_
#define _DENSITIES_
#define ARMA_WARN_LEVEL 1
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  double shape;
  double scale;
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  double b;
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  double num;
  double den;
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
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
//'@param varcovM (nbStates * (p + (p * p - p)/2 )) x n matrix of all sds (sd.x.state1, sd.x.state2, ..., sd.y.state1, ...) and then correlations (corr.xy.state1, corr.xy.state2, ..., corr.xz.state1, corr.xz.state2, ...) from lower triangle of variance-covariance matrix.
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
  int kk = 0;
  
  for (int i=0; i < n; i++) {
    kk = 0;
    for(int k=0; k < p; k++){
      sigma(k,k) = varcovM(k,i) * varcovM(k,i);
      for(int j=0; j < k; j++){
        sigma(k,j) = varcovM(k,i) * varcovM(j,i) * varcovM(p+kk,i);
        sigma(j,k) = varcovM(j,i) * varcovM(k,i) * varcovM(p+kk,i);
        ++kk;
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
      //if(out(i) < DBL_MIN) out(i) = DBL_MIN;
    //} else{
    //  out(i) = -0.5*quadform + logSqrtDetvarcovM + constant;
    //}
    //Rprintf("i %d x %f y %f mean %f %f sigma %f %f %f %f out %f \n",i,xvec(0),xvec(1),mtemp(0),mtemp(1),sigma(0,0),sigma(1,0),sigma(0,1),sigma(1,1),out(i));
  }
  
  return out;
  
}

bool isInteger(double x, bool warn = true) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

inline bool is_large_int(double x) {
  if (x > std::numeric_limits<int>::max())
    return true;
  return false;
}

inline double to_dbl(int x) {
  return static_cast<double>(x);
}

inline int to_pos_int(double x) {
  if (x < 0.0 || ISNAN(x))
    Rcpp::stop("value cannot be coerced to integer");
  if (is_large_int(x))
    Rcpp::stop("value out of integer range");
  return static_cast<int>(x);
}

//bool isInteger(double x, bool warn = true);
//inline bool is_large_int(double x); 
//inline double to_dbl(int x);
//inline int to_pos_int(double x);

#define GETV(x, i)      x[i % x.size()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.n_rows, j)   // wrapped indexing of matrix

//' Categorical density function
//'
//' Probability density function of the categorical distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param prob success probability
//' @param foo Unused (for compatibility with template)
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dcat_rcpp(const NumericVector x, const arma::mat prob, const arma::mat foo) 
{
  
  if (x.size() < 1 || prob.n_rows < 1) {
    return NumericVector(0);
  }
  
  int Nmax = std::max(
    static_cast<int>(x.size()),
    static_cast<int>(prob.n_cols)
  );
  int k = prob.n_rows;
  int l = prob.n_cols;
  arma::colvec p(Nmax);
  double p_tot;
  
  //bool throw_warning = false;
  
  //if (k < 2)
  //  Rcpp::stop("number of columns in prob is < 2");
  
  arma::mat prob_tab = prob.t();
  
  for (int i = 0; i < l; i++) {
    p_tot = 0.0;
    for (int j = 0; j < k; j++) {
      p_tot += prob_tab(i, j);
#ifdef IEEE_754
      if (ISNAN(p_tot))
        break;
#endif
      if (prob_tab(i, j) < 0.0) {
        p_tot = NAN;
        //throw_warning = true;
        break;
      }
    }
    for (int j = 0; j < k; j++)
      prob_tab(i, j) /= p_tot;
  }
  
  for (int i = 0; i < Nmax; i++) {
#ifdef IEEE_754
    if (ISNAN(GETV(x, i))) {
      p[i] = GETV(x, i);
      continue;
    }
#endif
    if (!isInteger(GETV(x, i)) || GETV(x, i) < 1.0 ||
        GETV(x, i) > to_dbl(k)) {
      p[i] = 0.0;
      continue;
    }
    if (is_large_int(GETV(x, i))) {
      //Rcpp::warning("NAs introduced by coercion to integer range in dcat_rcpp");
      p[i] = NA_REAL;
    }
    p[i] = GETM(prob_tab, i, to_pos_int(GETV(x, i)) - 1);
  }
  
  //if (log_prob)
  //  p = Rcpp::log(p);
  
  //if (throw_warning)
  //  Rcpp::warning("NaNs produced in dcat_rcpp");
  
  return p;
}

//' negative binomial density function
//'
//' Probability density function of the negative binomial distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param mu Mean of the distribution 
//' @param size Dispersion parameter
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dnbinom_rcpp(NumericVector x, arma::mat mu, arma::mat size)
{
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dnbinom_mu(x(i),size(i),mu(i),0);
  }
  
  return res;
}

//' logistic density function
//'
//' Probability density function of the logistic distribution (written in C++)
//'
//' @param x Vector of quantiles
//' @param location mean of the distribution 
//' @param scale Dispersion parameter
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dlogis_rcpp(NumericVector x, arma::mat location, arma::mat scale)
{
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dlogis(x(i),location(i),scale(i),0);
  }
  
  return res;
}

//' student t density function
//'
//' Probability density function of non-central student t (written in C++)
//'
//' @param x Vector of quantiles
//' @param df degrees of freedom 
//' @param ncp non-centrality parameter
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dt_rcpp(NumericVector x, arma::mat df, arma::mat ncp)
{
  unsigned int xs = (unsigned int) x.size();
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 1; // if missing observation
    else
      res(i) = R::dnt(x(i),df(i),ncp(i),0);
  }
  
  return res;
}

//' Correlated random walk Rice distribution
//'
//' Probability density function of Rice distribution under the CTCRW model of Johnson et al. (2008) (written in C++)
//'
//' @param x numeric data vector of length \code{n + n} where the first \code{n} entries correspond to the step lengths and the next \code{n} entries to the corresponding previous step lengths (the first of which is NA and ignored).
//' @param beta correlation parameter 
//' @param sigma speed parameter
//'
//' @return Vector of densities
// [[Rcpp::export]]
arma::colvec dcrwrice_rcpp(NumericVector x, arma::mat beta, arma::mat sigma)
{  
  double mu;
  double var;
  double sd;
  double xabs;
  double step_tm1;
  double besI;
  unsigned int xs = (unsigned int) beta.n_cols;
  arma::colvec res(xs);
  
  for(unsigned int i=0;i<xs;i++) {
    if(!arma::is_finite(x(i)))
      res(i) = 0.; // is missing observation
    else {
      step_tm1 = x(xs+i);
      mu = step_tm1*(1.-exp(-beta(i)))/beta(i);
      var = sigma(i)*sigma(i)/(beta(i)*beta(i)) * (1. - 2. / beta(i) * (1.-exp(-beta(i)))+1. / (2.*beta(i))*(1.-exp(-2.*beta(i))));
      sd = sqrt(var);
      xabs = abs(x(i)*mu/var);
      //if(xabs>709.){
      //  res(i) =  R_PosInf;
      //} else if(xabs<0.){
      //  res(i) =  R_NaN;
      //} else {
        besI = R::bessel_i(xabs,0,2);
        res(i) = log(x(i)) - 2.0 * log(sd) +
          (-(x(i)*x(i)+mu*mu)/(2.0*var)) +
          log(besI) + xabs;
        //if (!ISNAN(xabs) && !R_FINITE(besI)) {
        //  Rf_warning("Potential precision loss in dcrwrice R::bessel_i(%f, nu=0)",xabs);
        //}
      //}
    }
    //Rprintf("xs %d length(x) %d length(beta) %d length(sigma) %d i %d step_tm1 %f step %f beta %f sigma %f res %f \n",xs,x.size(),beta.n_cols,sigma.n_cols,i,step_tm1,x(i),beta(i),sigma(i),res(i));
  }
  return exp(res);
}

//' Correlated random walk von Mises density function
//'
//' Probability density function of the Von Mises distribution under the CTCRW model of Johnson et al. (2008) (written in C++)
//'
//' @param x numeric data vector of length \code{n + n + n} where the first \code{n} entries correspond to angles (von Mises distribution), the next \code{n} entries to the corresponding step lengths, and the last \code{n} entries to the corresponding previous step lengths.
//' @param beta correlation parameter 
//' @param sigma speed parameter
//'
//' @return Vector of densities
// [[Rcpp::export]]
 arma::colvec dcrwvm_rcpp(NumericVector x, arma::mat beta, arma::mat sigma)
{
 unsigned int xs = (unsigned int) beta.n_cols;
 arma::colvec res(xs);
 res(0) = 1.;
 
 double b;
 double kappa;
 double var;
 double step, step_tm1;
 for(unsigned int i=0;i<xs;i++) {
   if(!arma::is_finite(x(i)))
     res(i) = 1.; // is missing observation
   else {
     step = x(xs+i);
     step_tm1 = x(2*xs+i);
     var = sigma(i)*sigma(i)/(beta(i)*beta(i)) * (1. - 2. / beta(i) * (1.-exp(-beta(i)))+1. / (2.*beta(i))*(1-exp(-2*beta(i))));
     kappa = step*step_tm1*exp(-beta(i))/var;
     //if(kappa>709.){
     //  res(i) =  R_PosInf;
     //} else if(kappa<0.){
     //  res(i) =  R_NaN;
     //} else {
         b = R::bessel_i(kappa,0,2);
         res(i) = 1/(2*M_PI*b)*pow((exp(cos(x(i))-1)),kappa);
     //    if (!ISNAN(kappa) && !R_FINITE(b)) {
     //      Rf_warning("Potential precision loss in dcrwvm R::bessel_i(%f, nu=0)",kappa);
     //    }
     //  }
   }
 }
 
 return res;
}

// used in nLogLike_rcpp to map the functions' names to the functions
typedef arma::colvec (*FunPtr)(NumericVector, arma::mat, arma::mat);

#endif
