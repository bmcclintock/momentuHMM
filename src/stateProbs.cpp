#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat stateProbs_cpp(int nbObs, int nbStates, int mixtures,
                         arma::vec aInd, arma::mat eta,
                         Rcpp::List laList, Rcpp::List lbList) {
  
  arma::mat stateProbs = arma::zeros<arma::mat>(nbObs, nbStates);
  
  std::vector<arma::mat> la, lb;
  for(int mix = 0; mix < mixtures; mix++) {
    la.push_back(Rcpp::as<arma::mat>(laList[mix]));
    lb.push_back(Rcpp::as<arma::mat>(lbList[mix]));
  }
  
  arma::vec aInd0 = aInd - 1;
  
  arma::uvec animal_idx(nbObs);
  int current_animal = 0;
  for (int i = 0; i < nbObs; i++) {
    animal_idx(i) = current_animal;
    if (i == aInd0(current_animal) && current_animal < aInd0.n_elem - 1) {
      current_animal++;
    }
  }

  for (int i = 0; i < nbObs; i++) {
    int a = animal_idx(i); 
    
    for (int mix = 0; mix < mixtures; mix++) {
      double ieta = eta(a, mix);
      
      arma::rowvec log_alpha = la[mix].row(i);
      arma::rowvec log_beta = lb[mix].row(i);
      
      arma::rowvec log_prob = log_alpha + log_beta;
      
      double max_val = log_prob.max();
      
      if (!std::isfinite(max_val)) continue;
      
      arma::rowvec prob = arma::exp(log_prob - max_val);
      
      arma::rowvec norm_prob = prob / arma::sum(prob);
      stateProbs.row(i) += norm_prob * ieta;
    }
  }
  
  return stateProbs;
}