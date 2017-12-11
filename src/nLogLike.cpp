#include "densities.h"

//' Negative log-likelihood
//'
//' Computation of the negative log-likelihood (forward algorithm - written in C++)
//'
//' @param nbStates Number of states,
//' @param covs Covariates,
//' @param data A \code{\link{momentuHMMData}} object of the observations,
//' @param dataNames Character vector containing the names of the data streams,
//' @param dist Named list indicating the probability distributions of the data streams. 
//' @param Par Named list containing the state-dependent parameters of the data streams, matrix of regression coefficients 
//' for the transition probabilities ('beta'), and initial distribution ('delta').
//' @param aInd Vector of indices of the rows at which the data switches to another animal
//' @param zeroInflation Named list of logicals indicating whether the probability distributions of the data streams are zero-inflated.
//' @param oneInflation Named list of logicals indicating whether the probability distributions of the data streams are one-inflated.
//' @param stationary \code{false} if there are covariates. If \code{true}, the initial distribution is considered
//' equal to the stationary distribution. Default: \code{false}.
//' @param knownStates Vector of values of the state process which are known prior to fitting the
//' model (if any). Default: NULL (states are not known). This should be a vector with length the number
//' of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
//' the state is not known.
//' 
//' @return Negative log-likelihood
// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, arma::mat covs, DataFrame data, CharacterVector dataNames, List dist,
                     List Par,
                     IntegerVector aInd, List zeroInflation, List oneInflation,
                     bool stationary, IntegerVector knownStates)
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
  //arma::mat delta = Par["delta"];
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
  
  //=======================================================//
  // 1. Computation of initial distribution(s)             //
  //=======================================================//
  unsigned int nbAnimals = aInd.size();
  arma::mat delta(nbAnimals,nbStates);
  
  if(nbStates==1)
    delta.ones(); // no distribution if only one state
  else if(stationary) {
    // compute stationary distribution delta
    
    arma::mat diag(nbStates,nbStates);
    diag.eye(); // diagonal of ones
    arma::mat Gamma = trMat.slice(0).t(); // all slices are identical if stationary
    arma::colvec v(nbStates);
    v.ones(); // vector of ones
    arma::rowvec deltatmp(nbStates);
    try {
      deltatmp = arma::solve(diag-Gamma+1,v).t();
    }
    catch(...) {
      throw std::runtime_error("A problem occurred in the calculation of "
                                 "the stationary distribution. You may want to "
                                 "try different initial values and/or the option "
                                 "stationary=FALSE");
    }
    for(unsigned int k=0; k<nbAnimals; k++){
      delta.row(k) = deltatmp;
    }
  } else {
    arma::mat init = Par["delta"];
    delta = init;
    //arma::mat d(nbAnimals,nbStates);
    //arma::rowvec drowSums(nbAnimals);
    //drowSums.zeros();
    
    //d = covsDelta*init;
    
    //for(unsigned int k=0;k<nbAnimals;k++) {
    //  for(int i=0;i<nbStates;i++) {
    //    delta(k,i) = exp(d(k,i));
          
    //    // keep track of row sums, to normalize in the end
    //    drowSums(k) = drowSums(k) + delta(k,i);
    //  }
    //}
    
    //// normalization
    //for(unsigned int k=0;k<nbAnimals;k++) 
    //  for(int i=0;i<nbStates;i++) 
    //    delta(k,i) = delta(k,i)/drowSums(k);
  }

  //==========================================================//
  // 2. Computation of matrix of joint probabilities allProbs //
  //==========================================================//

  // map the functions names with the actual functions
  // (the type FunPtr and the density functions are defined in densities.h)
  map<std::string,FunPtr> funMap;
  funMap["bern"] = dbern_rcpp;
  funMap["beta"] = dbeta_rcpp;
  funMap["exp"] = dexp_rcpp;
  funMap["gamma"] = dgamma_rcpp;
  funMap["lnorm"] = dlnorm_rcpp;
  funMap["norm"] = dnorm_rcpp;
  funMap["pois"] = dpois_rcpp;
  funMap["vm"] = dvm_rcpp;
  funMap["weibull"] = dweibull_rcpp;
  funMap["wrpcauchy"] = dwrpcauchy_rcpp;

  arma::mat allProbs(nbObs,nbStates);
  allProbs.ones();

  arma::colvec genProb(nbObs);
  arma::rowvec zerom(nbObs);
  arma::rowvec onem(nbObs);
  bool genzeroInflation;
  bool genoneInflation;
  arma::mat genArgs1;
  arma::mat genArgs2;
  arma::mat zeromass(nbStates,nbObs);
  arma::uvec noZeros;
  arma::uvec nbZeros;
  arma::mat onemass(nbStates,nbObs);
  arma::uvec noOnes;
  arma::uvec nbOnes;
  arma::uvec noZerosOnes;
  
  double NAvalue = -99999999; // value designating NAs in data
  int zeroInd = 0;
  int oneInd = 0;
  
  unsigned int nDists = dist.size();

  for(unsigned int k=0;k<nDists;k++){
    genname = as<std::string>(dataNames[k]);
    genData = as<NumericVector>(data[genname]);
    genDist = as<std::string>(dist[genname]);
    genPar = as<arma::mat>(Par[genname]);
    genzeroInflation = as<bool>(zeroInflation[genname]);
    genoneInflation = as<bool>(oneInflation[genname]);
    
    if(genoneInflation) 
      oneInd = nbStates;
    else 
      oneInd = 0;
    
    if(genzeroInflation) 
      zeroInd = nbStates;
    else 
      zeroInd = 0;
    
    // remove the NAs from step (impossible to subset a vector with NAs)
    for(int i=0;i<nbObs;i++) {
      if(!arma::is_finite(genData(i))) {
        genData(i) = NAvalue;
      }
    }
    arma::uvec noNAs = arma::find(as<arma::vec>(genData)!=NAvalue);
    
    // extract zero-mass and one-mass parameters if necessary
    if(genzeroInflation || genoneInflation) {
      
      if(genzeroInflation){
        zeromass = genPar.rows(genPar.n_rows-oneInd-nbStates,genPar.n_rows-oneInd-1);   //genPar(arma::span(genPar.n_rows-1),arma::span(),arma::span());
        
        noZeros = arma::find(as<arma::vec>(genData)>0);
        nbZeros = arma::find(as<arma::vec>(genData)==0);
      }
      
      if(genoneInflation){
        onemass = genPar.rows(genPar.n_rows-nbStates,genPar.n_rows-1);   //genPar(arma::span(genPar.n_rows-1),arma::span(),arma::span());
        noOnes = arma::find(as<arma::vec>(genData)!=NAvalue && as<arma::vec>(genData)<1);
        nbOnes = arma::find(as<arma::vec>(genData)==1);
      }
      
      if(genzeroInflation && genoneInflation)
        noZerosOnes = arma::find(as<arma::vec>(genData)>0 && as<arma::vec>(genData)<1);
      
      genPar = genPar.rows(0,genPar.n_rows-oneInd-zeroInd-1); //genPar.tube(0, 0, genPar.n_rows-2, genPar.n_cols-1);
    }

    for(int state=0;state<nbStates;state++){
      
      genProb.ones();
      
      genArgs1 = genPar.row(state); //genPar(arma::span(0),arma::span(state),arma::span());
      genArgs2 = genPar.row(genPar.n_rows - nbStates + state); //genPar(arma::span(genPar.n_rows-1),arma::span(state),arma::span());
    
      if(genzeroInflation && !genoneInflation) {
        
        zerom = zeromass.row(state);

        // compute probability of non-zero observations
        genProb.elem(noZeros) = (1. - zerom.elem(noZeros)) % funMap[genDist](genData[genData>0],genArgs1.elem(noZeros),genArgs2.elem(noZeros));

        // compute probability of zero observations
        genProb.elem(nbZeros) = zerom.elem(nbZeros);
        
      } else if(genoneInflation && !genzeroInflation){
        
        onem = onemass.row(state);
        
        // compute probability of non-one observations
        genProb.elem(noOnes) = (1. - onem.elem(noOnes)) % funMap[genDist](genData[(genData!=NAvalue) & (genData<1)],genArgs1.elem(noOnes),genArgs2.elem(noOnes));
        
        // compute probability of one observations
        genProb.elem(nbOnes) = onem.elem(nbOnes);
        
      } else if(genzeroInflation && genoneInflation){
        
        zerom = zeromass.row(state);
        onem = onemass.row(state);
        
        // compute probability of non-zero and non-one observations
        genProb.elem(noZerosOnes) = (1. - zerom.elem(noZerosOnes) - onem.elem(noZerosOnes)) % funMap[genDist](genData[(genData>0) & (genData<1)],genArgs1.elem(noZerosOnes),genArgs2.elem(noZerosOnes));
        
        // compute probability of zero observations
        genProb.elem(nbZeros) = zerom.elem(nbZeros);
        
        // compute probability of one observations
        genProb.elem(nbOnes) = onem.elem(nbOnes);
        
      } else {
        
        genProb.elem(noNAs) = funMap[genDist](genData[genData!=NAvalue],genArgs1.elem(noNAs),genArgs2.elem(noNAs));
        
      }
      
      allProbs.col(state) = allProbs.col(state) % genProb;
    }
    
    // put the NAs back
    genData[genData==NAvalue] = NA_REAL;
  }
  
  // deal with states known a priori
  double prob = 0;
  if(knownStates(0) != -1) {
    // loop over rows
    for(unsigned int i=0 ; i<allProbs.n_rows ; i++) {
      if(knownStates(i)>0) {
        prob = allProbs(i,knownStates(i)-1); // save non-zero probability
        allProbs.row(i).zeros(); // set other probabilities to zero
        allProbs(i,knownStates(i)-1) = prob;
      }
    }
  }

  //======================//
  // 3. Forward algorithm //
  //======================//

  arma::mat Gamma(nbStates,nbStates); // transition probability matrix
  double lscale = 0; // scaled log-likelihood
  unsigned int k=0; // animal index
  arma::rowvec alpha(nbStates);
  arma::rowvec delt(nbStates);
  for(unsigned int i=0;i<allProbs.n_rows;i++) {
    
    if(nbStates>1)
      Gamma = trMat.slice(i);
    else
      Gamma = 1; // no transition if only one state
    
    if(k<aInd.size() && i==(unsigned)(aInd(k)-1)) {
      // if 'i' is the 'k'-th element of 'aInd', switch to the next animal
      delt = delta.row(k);
      alpha = (delt * Gamma) % allProbs.row(i);
      k++;
    } else {
      alpha = (alpha * Gamma) % allProbs.row(i);
    }
    
    lscale = lscale + log(sum(alpha));
    alpha = alpha/sum(alpha);
  }

  return -lscale;
}
