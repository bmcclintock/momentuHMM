#include "densities.h"

//' Negative log-likelihood
//'
//' Computation of the negative log-likelihood (forward algorithm - written in C++)
//'
//' @param nbStates Number of states
//' @param beta Matrix of regression coefficients for the transition probabilities
//' @param covs Covariates
//' @param data A \code{\link{moveData}} object of the observations
//' @param stepDist The name of the step length distribution
//' @param angleDist The name of the turning angle distribution
//' @param stepPar State-dependent parameters of the step length distribution
//' @param anglePar State-dependent parameters of the turning angle distribution
//' @param delta Stationary distribution
//' @param aInd Vector of indices of the rows at which the data switches to another animal
//' @param zeroInflation \code{true} if zero-inflation is included in the step length distribution,
//' \code{false} otherwise.
//' @param stationary \code{false} if there are covariates. If \code{true}, the initial distribution is considered
//' equal to the stationary distribution. Default: \code{false}.
//'
//' @return Negative log-likelihood
// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, arma::mat covs, DataFrame data, CharacterVector dataNames, List dist,
                     List Par,
                     IntegerVector aInd, List zeroInflation,
                     bool stationary=false)
{
  int nbObs = data.nrows();

  //=======================================================//
  // 1. Computation of transition probability matrix trMat //
  //=======================================================//

  arma::cube trMat(nbStates,nbStates,nbObs);
  trMat.zeros();
  arma::mat rowSums(nbStates,nbObs);
  rowSums.zeros();

  arma::mat g(nbObs,nbStates*(nbStates-1));
  
  arma::mat beta = Par["beta"];
  arma::rowvec delta = Par["delta"];
  NumericVector genData(nbObs);
  arma::mat genPar;
  std::string genDist;
  std::string genname;
  
  if(nbStates>1) {
    g = covs*beta;

    for(int k=0;k<nbObs;k++) {
      int cpt=0; // counter for diagonal elements
      for(int i=0;i<nbStates;i++) {
        for(int j=0;j<nbStates;j++) {
          if(i==j) {
            // if diagonal element, set to one and increment counter
            trMat(i,j,k)=1;
            cpt++;
          }
          else
            trMat(i,j,k) = exp(g(k,i*nbStates+j-cpt));

          // keep track of row sums, to normalize in the end
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

  //==========================================================//
  // 2. Computation of matrix of joint probabilities allProbs //
  //==========================================================//

  // map the functions names with the actual functions
  // (the type FunPtr and the density functions are defined in densities.h)
  map<std::string,FunPtr> funMap;
  funMap["gamma"] = dgamma_rcpp;
  funMap["weibull"] = dweibull_rcpp;
  funMap["lnorm"] = dlnorm_rcpp;
  funMap["exp"] = dexp_rcpp;
  funMap["vm"] = dvm_rcpp;
  funMap["wrpcauchy"] = dwrpcauchy_rcpp;
  funMap["beta"] = dbeta_rcpp;
  funMap["pois"] = dpois_rcpp;

  if(nbStates==1)
    delta = 1; // no distribution if only one state
  else if(stationary) {
    // compute stationary distribution delta

    arma::mat diag(nbStates,nbStates);
    diag.eye(); // diagonal of ones
    arma::mat Gamma = trMat.slice(0); // all slices are identical if stationary
    arma::colvec v(nbStates);
    v.ones(); // vector of ones
    try {
      delta = arma::solve(diag-Gamma+1,v).t();
    }
    catch(...) {
      throw std::runtime_error("A problem occurred in the calculation of "
                                 "the stationary distribution. You may want to "
                                 "try different initial values and/or the option "
                                 "stationary=FALSE");
    }
  }

  arma::mat allProbs(nbObs,nbStates);
  allProbs.ones();

  arma::colvec genProb(nbObs);
  NumericVector genArgs(2); 
  NumericVector gendata(nbObs);
  bool genzeroInflation;
  arma::rowvec zeromass(nbStates);

  // extract zero-mass parameters from step parameters if necessary


  for(int state=0;state<nbStates;state++)
  {
    for(unsigned int i=0;i<dist.size();i++){
      genProb.ones();
      genname = as<std::string>(dataNames[i]);
      genData = as<NumericVector>(data[genname]);
      genDist = as<std::string>(dist[genname]);
      genPar = as<arma::mat>(Par[genname]);
      genzeroInflation = as<bool>(zeroInflation[genname]);
      
      if(genzeroInflation) {
        zeromass = genPar.row(genPar.n_rows-1);
        genPar = genPar.submat(0,0,genPar.n_rows-2,genPar.n_cols-1);
      }
      
      for(unsigned int j=0;j<genPar.n_rows;j++)
        genArgs(j) = genPar(j,state);
      
      if(genzeroInflation) {
        // remove the NAs from step (impossible to subset a vector with NAs)
        for(int i=0;i<nbObs;i++) {
          if(!arma::is_finite(genData(i))) {
            genData(i) = -1;
            genProb(i) = 1;
          }
        }
        
        // compute probability of non-zero observations
        genProb.elem(arma::find(as<arma::vec>(genData)>0)) = (1-zeromass(state))*funMap[genDist](genData[genData>0],genArgs(0),genArgs(1));
        
        // compute probability of zero observations
        int nbZeros = as<NumericVector>(genData[genData==0]).size();
        arma::vec zm(nbZeros);
        for(int i=0;i<nbZeros;i++)
          zm(i) = zeromass(state);
        genProb.elem(arma::find(as<arma::vec>(genData)==0)) = zm;
        
        // put the NAs back
        for(int i=0;i<nbObs;i++) {
          if(genData(i)<0)
            genData(i) = NA_REAL;
        }
      } else {
        genProb = funMap[genDist](data[genname],genArgs(0),genArgs(1));
      }
      allProbs.col(state) = allProbs.col(state)%genProb;
    }
  }

  //======================//
  // 3. Forward algorithm //
  //======================//

  arma::mat Gamma(nbStates,nbStates); // transition probability matrix
  double lscale = 0; // scaled log-likelihood
  int k=1; // animal index
  arma::rowvec alpha = delta%allProbs.row(0);

  for(unsigned int i=1;i<allProbs.n_rows;i++) {
    if(k<aInd.size() && i==(unsigned)(aInd(k)-1)) {
      // if 'i' is the 'k'-th element of 'aInd', switch to the next animal
      k++;
      alpha = delta%allProbs.row(i);
    }

    if(nbStates>1)
      Gamma = trMat.slice(i);
    else
      Gamma = 1; // no transition if only one state

    alpha = alpha*Gamma%allProbs.row(i);

    lscale = lscale + log(sum(alpha));
    alpha = alpha/sum(alpha);
  }

  return -lscale;
}
