#include "densities.h"
#include "combine.h"
#include "expmatrix.h"
#include "stationary.h"

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
//' @param stationary \code{false} if there are time-varying covariates in \code{formula} or any covariates in \code{formulaDelta}. If \code{true}, the initial distribution is considered
//' equal to the stationary distribution. Default: \code{false}.
//' @param knownStates Vector of values of the state process which are known prior to fitting the
//' model (if any). Default: NULL (states are not known). This should be a vector with length the number
//' of rows of 'data'; each element should either be an integer (the value of the known states) or NA if
//' the state is not known.
//' @param betaRef Indices of reference elements for t.p.m. multinomial logit link.
//' @param mixtures Number of mixtures for the state transition probabilities
//' @param dtIndex time difference index for calculating transition probabilities of hierarchical continuous-time models
//' @param CT logical indicating whether to fit discrete-time approximation of a continuous-time model
//' 
//' @return Negative log-likelihood
// [[Rcpp::export]]
double nLogLike_rcpp(int nbStates, arma::mat covs, DataFrame data, CharacterVector dataNames, List dist,
                     List Par,
                     IntegerVector aInd, List zeroInflation, List oneInflation,
                     bool stationary, IntegerVector knownStates, IntegerVector betaRef, int mixtures, IntegerVector dtIndex, bool CT = false)
{
  int nbObs = data.nrows();

  //=======================================================//
  // 1. Computation of transition probability matrix trMat //
  //=======================================================//

  std::vector<arma::cube> trMat(mixtures);
  for(int mix=0; mix<mixtures; mix++){
    trMat[mix] = arma::cube(nbStates,nbStates,nbObs);
    trMat[mix].zeros();
  }

  arma::mat rowSums(nbStates,nbObs);
  rowSums.zeros();

  arma::mat g(nbObs,nbStates*(nbStates-1));
  
  arma::mat betaMix = Par["beta"];
  //arma::mat delta = Par["delta"];
  NumericVector genData;
  List L;
  arma::mat genPar;
  std::string genDist;
  std::string genname;
  int nbCovs = 0;
  NumericVector dt;
  
  if(nbStates>1) {
    
    nbCovs = betaMix.n_rows * betaMix.n_cols / (mixtures * nbStates * (nbStates-1)) - 1;
    arma::mat beta(nbCovs+1,nbStates*(nbStates-1));
    
    for(int mix=0; mix<mixtures; mix++){    
      
      beta = betaMix.rows(mix*(nbCovs+1)+0,mix*(nbCovs+1)+nbCovs);
      g = covs*beta;
      
      if(CT){ // continuous-time HMM approximation
        dt = data["dt"];
        //g = exp(g);
        for(int k=0;k<nbObs;k++) {
          int cpt=0; // counter for diagonal elements
          for(int i=0;i<nbStates;i++) {
            for(int j=0;j<nbStates;j++) {
              if(j==(betaRef(i)-1)) {
                cpt++;
              } else {
                if(i!=j) trMat[mix](i,j,k) = exp(g(k,i*nbStates+j-cpt));
                else trMat[mix](i,j,k) = -exp(-g(k,i*nbStates+j-cpt));
                //Rprintf("k %d i %d j %d i*nbStates+j-cpt %d g %f trMat[mix](i,j,k) %f dt %f \n",k,i,j,i*nbStates+j-cpt,g(k,i*nbStates+j-cpt),trMat[mix](i,j,k),dt(k-kStart));
              }
            }
            for(int l=0;l<nbStates;l++){
              if((betaRef(i)-1)!=l) trMat[mix](i,(betaRef(i)-1),k) -= trMat[mix](i,l,k);
            }
            //Rprintf("k %d i %d j %d trMat[mix](i,j,k) %f dt %f \n",k,i,i,trMat[mix](i,i,k),dt(k-kStart));
          }
        }        
      } else { // standard discrete-time HMM
        for(int k=0;k<nbObs;k++) {
          int cpt=0; // counter for diagonal elements
          for(int i=0;i<nbStates;i++) {
            for(int j=0;j<nbStates;j++) {
              if(j==(betaRef(i)-1)) {
                // if reference element, set to one and increment counter
                trMat[mix](i,j,k)=1;
                cpt++;
              }
              else
                trMat[mix](i,j,k) = exp(g(k,i*nbStates+j-cpt));
    
              // keep track of row sums, to normalize in the end
              rowSums(i,k)=rowSums(i,k)+trMat[mix](i,j,k);
            }
          }
        }
      
        // normalization
        for(int k=0;k<nbObs;k++)
          for(int i=0;i<nbStates;i++)
            for(int j=0;j<nbStates;j++)
              trMat[mix](i,j,k) = trMat[mix](i,j,k)/rowSums(i,k);
        
        rowSums.zeros();
      }
    }
  }
  
  //=======================================================//
  // 1. Computation of initial distribution(s)             //
  //=======================================================//
  unsigned int nbAnimals = (unsigned int) aInd.size();
  std::vector<arma::mat> delta(mixtures);
  for(int mix=0; mix<mixtures; mix++)
    delta[mix] = arma::mat(nbAnimals,nbStates);
  
  if(nbStates==1){
    for(int mix=0; mix<mixtures; mix++)
      delta[mix].ones(); // no distribution if only one state
  } else if(stationary) {
    // compute stationary distribution delta
    arma::mat diag(nbStates,nbStates);
    diag.eye(); // diagonal of ones
    arma::mat Gamma(nbStates,nbStates); // all slices are identical if stationary
    arma::colvec v(nbStates);
    v.ones(); // vector of ones
    arma::rowvec deltatmp(nbStates);
    for(int mix=0; mix<mixtures; mix++){
      if(!nbCovs){
        Gamma = trMat[mix].slice(0); // all slices are identical if stationary
        try {
          if(!CT) deltatmp = arma::solve(diag-Gamma.t()+1,v).t();
          else deltatmp = stationary_rcpp(Gamma);
        }
        catch(...) {
          throw std::runtime_error("A problem occurred in the calculation of "
                                     "the stationary distribution. You may want to "
                                     "try different initial values and/or the option "
                                     "stationary=FALSE");
        }
        for(unsigned int k=0; k<nbAnimals; k++){
          delta[mix].row(k) = deltatmp;
        }
      } else{
        for(unsigned int k=0; k<nbAnimals; k++){
          Gamma = trMat[mix].slice(aInd[k]); // all slices are identical for each individual if stationary
          try {
            if(!CT) deltatmp = arma::solve(diag-Gamma.t()+1,v).t();
            else deltatmp = stationary_rcpp(Gamma);
          }
          catch(...) {
            throw std::runtime_error("A problem occurred in the calculation of "
                                       "the stationary distribution. You may want to "
                                       "try different initial values and/or the option "
                                       "stationary=FALSE");
          }
          delta[mix].row(k) = deltatmp;
        }
      }
    }
    
  } else {
    arma::mat init = Par["delta"];
    for(int mix=0; mix<mixtures; mix++)
      delta[mix] = init.rows(mix*nbAnimals,mix*nbAnimals+nbAnimals-1);
  }

  //==========================================================//
  // 2. Computation of matrix of joint probabilities allProbs //
  //==========================================================//

  // map the functions names with the actual functions
  // (the type FunPtr and the density functions are defined in densities.h)
  map<std::string,FunPtr> funMap;
  funMap["bern"] = dbern_rcpp;
  funMap["beta"] = dbeta_rcpp;
  funMap["cat"] = dcat_rcpp;
  funMap["ctds"] = dcat_rcpp;
  funMap["exp"] = dexp_rcpp;
  funMap["gamma"] = dgamma_rcpp;
  funMap["logis"] = dlogis_rcpp;
  funMap["lnorm"] = dlnorm_rcpp;
  funMap["negbinom"] = dnbinom_rcpp;
  funMap["norm"] = dnorm_rcpp;
  funMap["rw_norm"] = dnorm_rcpp;
  funMap["mvnorm2"] = dmvnorm_rcpp;
  funMap["rw_mvnorm2"] = dmvnorm_rcpp;
  funMap["mvnorm3"] = dmvnorm_rcpp;
  funMap["rw_mvnorm3"] = dmvnorm_rcpp;
  funMap["pois"] = dpois_rcpp;
  funMap["t"] = dt_rcpp;
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
  
  double NAvalue = -99999999999; // value designating NAs in data
  int zeroInd = 0;
  int oneInd = 0;
  bool mvn = false;
  
  unsigned int nDists = (unsigned int) dist.size();

  for(unsigned int k=0;k<nDists;k++){
    genname = as<std::string>(dataNames[k]);
    genDist = as<std::string>(dist[genname]);
    if(genDist=="mvnorm2" || genDist=="rw_mvnorm2"){
      L = List::create(as<NumericVector>(data[genname+".x"]) ,as<NumericVector>(data[genname+".y"]));
      //NumericVector genData = combine(L);
      mvn = true;
    } else if(genDist=="mvnorm3" || genDist=="rw_mvnorm3"){
      L = List::create(as<NumericVector>(data[genname+".x"]) ,as<NumericVector>(data[genname+".y"]),as<NumericVector>(data[genname+".z"]));
      //NumericVector genData = combine(L);
      mvn = true;
    } else {
      L = List::create(as<NumericVector>(data[genname]));
      mvn = false;
    }
    //genData = combine(L);
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
    NumericVector tmpL;
    for(int l=0; l<L.size(); l++){
      tmpL = L[l];
      for(int i=0;i<nbObs;i++) {
        if(!arma::is_finite(tmpL(i))) {
          //genData(i) = NAvalue;
          tmpL(i) = NAvalue;
        }
      }
      L[l] = tmpL;
    }
    arma::uvec noNAs = arma::find(as<arma::vec>(tmpL)!=NAvalue);
    genData = combine(L);
    
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
      
      if(mvn){
        if(genDist=="mvnorm2" || genDist=="rw_mvnorm2"){
          genArgs1 = cbindmean2(genPar.row(state),genPar.row(nbStates+state));
          genArgs2 = cbindsigma2(genPar.row(nbStates*2+state),genPar.row(nbStates*3+state),genPar.row(nbStates*4+state));
          //for(int i=0; i<genArgs1.n_cols; i++){
            //for(int j=0; j<genArgs1.n_rows; j++){
              //Rprintf("i %d (%f,%f) mean(%f,%f) sigma(%f,%f,%f,%f) noNAs %d \n",i,genData(i),genData(nbObs+i),genArgs1(0,i),genArgs1(1,i),genArgs2(0,i),genArgs2(1,i),genArgs2(2,i),genArgs2(3,i));
            //}
          //}
        } else if(genDist=="mvnorm3" || genDist=="rw_mvnorm3"){
          genArgs1 = cbindmean3(genPar.row(state),genPar.row(nbStates+state),genPar.row(nbStates*2+state));
          genArgs2 = cbindsigma3(genPar.row(nbStates*3+state),genPar.row(nbStates*4+state),genPar.row(nbStates*5+state),genPar.row(nbStates*6+state),genPar.row(nbStates*7+state),genPar.row(nbStates*8+state));          
        }
      } else if((genDist=="cat") | (genDist=="ctds")){
        int catDim = genPar.n_rows / nbStates;
        arma::mat tmpPar1(catDim,nbObs);
        for(int l=0; l<catDim; l++){
          tmpPar1.row(l) = genPar.row(l*nbStates+state);
        }
        genArgs1 = tmpPar1;
        genArgs2 = tmpPar1;
      } else {
        genArgs1 = genPar.row(state); //genPar(arma::span(0),arma::span(state),arma::span());
        genArgs2 = genPar.row(genPar.n_rows - nbStates + state); //genPar(arma::span(genPar.n_rows-1),arma::span(state),arma::span());
      }
      
      if(genzeroInflation && !genoneInflation) {
        
        zerom = zeromass.row(state);

        // compute probability of non-zero observations
        genProb.elem(noZeros) = (1. - zerom.elem(noZeros)) % funMap[genDist](genData[genData>0.],genArgs1.elem(noZeros),genArgs2.elem(noZeros));

        // compute probability of zero observations
        genProb.elem(nbZeros) = zerom.elem(nbZeros);
        
      } else if(genoneInflation && !genzeroInflation){
        
        onem = onemass.row(state);
        
        // compute probability of non-one observations
        genProb.elem(noOnes) = (1. - onem.elem(noOnes)) % funMap[genDist](genData[(genData!=NAvalue) & (genData<1.)],genArgs1.elem(noOnes),genArgs2.elem(noOnes));
        
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
        genProb.elem(noNAs) = funMap[genDist](genData[genData!=NAvalue],genArgs1.cols(noNAs),genArgs2.cols(noNAs));
      }
      
      allProbs.col(state) = allProbs.col(state) % genProb;
      //for(int i=0; i<nbObs; i++){
      //  Rprintf("allProbs state %d i %d %f \n",state,i,allProbs(i,state));
      //}
    }
    
    // put the NAs back and check for underflow to zero
    int badInd;
    for(int l=0; l<L.size(); l++){
      tmpL = L[l];
      for(int i=0;i<nbObs;i++) {
        if(tmpL(i)==NAvalue) {
          tmpL(i) = NA_REAL;
        } else {
          badInd = 0;
          for(int state=0; state<nbStates; state++){
            if(allProbs(i,state) < DBL_MIN){
              badInd+=1;
            }
          }
          if(badInd==nbStates){
            for(int state=0; state<nbStates; state++){
              allProbs(i,state) = DBL_MIN;
            }
            //Rprintf("Warning: '%s' probability density for observation %d is zero and results in numerical underflow; this probability has been set to %f\n",genDist.c_str(),i+1,DBL_MIN);
          }
        }
      }
      L[l] = tmpL;
    }
    genData = combine(L);
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
  arma::rowvec mixlscale(mixtures);
  mixlscale.zeros();
  unsigned int k=0; // animal index
  arma::mat alpha(mixtures,nbStates);
  arma::rowvec delt(nbStates);
  arma::mat pie = Par["pi"];
  
  double A; 
  double lscale = 0.0;
  arma::rowvec maxscale(mixtures);
  arma::mat Gammat(nbStates,nbStates); 
  arma::mat Gamma0(nbStates,nbStates); 
  Gamma0.eye(); // diagonal of ones
  
  for(unsigned int i=0;i<allProbs.n_rows;i++) {
      
    for(int mix=0; mix<mixtures; mix++){
      
      if(nbStates>1){
        if(CT) {
          if(k<nbAnimals && i==(unsigned)(aInd(k)-1)) {
            Gamma = Gamma0;
          } else {
            try {
              Gammat = trMat[mix].slice(i) * dt(dtIndex(i-1));
              Gamma = expmatrix_rcpp(Gammat);
              //Rprintf("i %d dtIndex %d dt %f Gamma %f %f %f %f \n",i,index,dt(dtIndex(i-1)),Gamma(0,0),Gamma(0,1),Gamma(1,0),Gamma(1,1));
            }
            catch(std::exception &ex) {	
              forward_exception_to_r(ex);
            }
          }
        } else Gamma = trMat[mix].slice(i);
      } else
        Gamma = 1; // no transition if only one state
      
      if(k<nbAnimals && i==(unsigned)(aInd(k)-1)) {
        // if 'i' is the 'k'-th element of 'aInd', switch to the next animal
        delt = delta[mix].row(k);
        alpha.row(mix) = (delt * Gamma) % allProbs.row(i);
        mixlscale(mix) = 0;
        //if(mix==(mixtures-1)) k++;
      } else {
        alpha.row(mix) = (alpha.row(mix) * Gamma) % allProbs.row(i);
      }
      
      mixlscale(mix) += log(sum(alpha.row(mix)));
      alpha.row(mix) = alpha.row(mix)/sum(alpha.row(mix));
    }
    if((k+1<nbAnimals && i==(unsigned)(aInd(k+1)-2)) || (i==(allProbs.n_rows-1))){
      maxscale = mixlscale + log(pie.row(k));
      A = maxscale.max(); // cancels out below; helps prevent numerical issues
      lscale += A + log(sum(exp(maxscale - A)));
      k++;
    }
  }
  //mixlscale += log(pie) * nbAnimals;


  return -lscale;
}
