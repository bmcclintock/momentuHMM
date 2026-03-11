#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::uvec viterbi_cpp(int nbStates, int nbAnimals, int mixtures,
                            arma::vec aInd, arma::mat mixProbs, arma::mat delta,
                            arma::mat probs, Rcpp::List trMatList) {
  
  int nTotal = probs.n_rows;
  arma::uvec allStates(nTotal);

  std::vector<arma::cube> trMat;
  for(int mix = 0; mix < mixtures; mix++) {
    trMat.push_back(Rcpp::as<arma::cube>(trMatList[mix]));
  }
  
  for (int zoo = 0; zoo < nbAnimals; zoo++) {

    int start_idx = aInd(zoo) - 1; 
    int end_idx;
    
    if (zoo != nbAnimals - 1) {
      end_idx = aInd(zoo + 1) - 2; 
    } else {
      end_idx = nTotal - 1;
    }
    
    int nbObs = end_idx - start_idx + 1;
    
    arma::mat tmxi = arma::zeros<arma::mat>(nbObs, nbStates);
    arma::mat xi_mix = arma::zeros<arma::mat>(nbObs, nbStates);
    
    arma::rowvec foo = arma::zeros<arma::rowvec>(nbStates);
    for(int mix = 0; mix < mixtures; mix++) {
      arma::rowvec d = delta.row(mix * nbAnimals + zoo);
      arma::mat tm_1 = trMat[mix].slice(start_idx);
      arma::rowvec p_1 = probs.row(start_idx);
      foo += mixProbs(zoo, mix) * (d * tm_1) % p_1;
    }
    xi_mix.row(0) = foo / arma::sum(foo);
    
    for (int i = 1; i < nbObs; i++) {
      arma::mat foo_mat = arma::zeros<arma::mat>(nbStates, nbStates);
      for(int mix = 0; mix < mixtures; mix++) {
        arma::mat tm_i = trMat[mix].slice(start_idx + i);
        arma::mat step_mat = tm_i;

        step_mat.each_col() %= xi_mix.row(i - 1).t();
        foo_mat += mixProbs(zoo, mix) * step_mat;
      }
      arma::rowvec max_vals = arma::max(foo_mat, 0); 
      arma::rowvec p_i = probs.row(start_idx + i);
      foo = max_vals % p_i;
      xi_mix.row(i) = foo / arma::sum(foo);
    }
    
    arma::uvec stSeq(nbObs);
    tmxi.row(nbObs - 1) = xi_mix.row(nbObs - 1);
    stSeq(nbObs - 1) = tmxi.row(nbObs - 1).index_max();
    
    for (int i = nbObs - 2; i >= 0; i--) {
      arma::rowvec current_tmxi = arma::zeros<arma::rowvec>(nbStates);
      for (int mix = 0; mix < mixtures; mix++) {
        arma::vec tm_col = trMat[mix].slice(start_idx + i + 1).col(stSeq(i + 1));
        current_tmxi += xi_mix.row(i) % tm_col.t() * mixProbs(zoo, mix);
      }
      tmxi.row(i) = current_tmxi;
      stSeq(i) = current_tmxi.index_max();
    }
    
    for (int i = 0; i < nbObs; i++) {
      allStates(start_idx + i) = stSeq(i) + 1;
    }
  }
  
  return allStates;
}