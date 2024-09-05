// code modified from the hmmTMB package: Michelot T, Glennie R (2023). hmmTMB: Fit Hidden Markov Models using Template Model Builder. R package version 1.0.2, <https://CRAN.R-project.org/package=hmmTMB>.

/// @file momentuhmmTMB.hpp

#ifndef momentuhmmTMB_hpp
#define momentuhmmTMB_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

#include "include/dist.hpp"

// negative log-likelihood of the HMM
template<class Type>
Type momentuhmmTMB(objective_function<Type>* obj) {
  // DATA
  DATA_VECTOR(ID); // vector of time series IDs
  DATA_MATRIX(data); // data stream
  DATA_IVECTOR(datadim); // dimension of observations for each variable 
  DATA_MATRIX(known_states); // known states n observations x states 1 = possible, 0 = impossible
  DATA_INTEGER(n_states); // number of states
  DATA_INTEGER(statdist); // use stationary distribution with respect to first tpm 
  DATA_IVECTOR(distcode); // codes of observation distributions
  DATA_IVECTOR(distpar); // number of parameters for observation distributions
  DATA_INTEGER(CT); // indicator for continuous time
  DATA_VECTOR(delta_t); // time interval between observations
  DATA_IVECTOR(dtIndex); // time index for dt
  DATA_SCALAR(kappa); // max transition rate
  // model matrices for observation process
  DATA_SPARSE_MATRIX(X_fe_obs); // design matrix for fixed effects
  DATA_SPARSE_MATRIX(X_re_obs); // design matrix for random effects
  DATA_SPARSE_MATRIX(S_obs); // penalty matrix
  DATA_IMATRIX(ncol_re_obs); // number of columns of S and X_re for each random effect
  // model matrices for hidden state process
  DATA_SPARSE_MATRIX(X_fe_hid); // design matrix for fixed effects
  DATA_SPARSE_MATRIX(X_re_hid); // design matrix for random effects
  DATA_SPARSE_MATRIX(S_hid); // penalty matrix
  DATA_IMATRIX(ncol_re_hid); // number of columns of S and X_re for each random effect
  DATA_INTEGER(include_smooths); // > 0 = include penalty in likelihood evaluation
  DATA_IVECTOR(ref_tpm); // indices of reference transition probabilities
  // prior information 
  DATA_MATRIX(coeff_fe_obs_prior); // means, sds for prior on fixed effects for obs 
  DATA_MATRIX(coeff_fe_hid_prior); // means, sds for prior on fixed effects for hidden 
  DATA_MATRIX(log_lambda_obs_prior); // means, sds for prior on smoothing parameters for obs 
  DATA_MATRIX(log_lambda_hid_prior); // means, sds for prior on smoothing parameters for hidden 
  // working parameter bounds 
  DATA_VECTOR(lower_fe_obs); // observation parameters (fixed effects)
  DATA_VECTOR(lower_lambda_obs); // smoothness parameters
  DATA_VECTOR(lower_fe_hid); // state process parameters (fixed effects)
  DATA_VECTOR(lower_lambda_hid); // smoothness parameters
  DATA_VECTOR(lower_delta0); // initial distribution
  DATA_VECTOR(lower_re_obs); // observation parameters (random effects)
  DATA_VECTOR(lower_re_hid); // state process parameters (random effects)
  DATA_VECTOR(upper_fe_obs); // observation parameters (fixed effects)
  DATA_VECTOR(upper_lambda_obs); // smoothness parameters
  DATA_VECTOR(upper_fe_hid); // state process parameters (fixed effects)
  DATA_VECTOR(upper_lambda_hid); // smoothness parameters
  DATA_VECTOR(upper_delta0); // initial distribution
  DATA_VECTOR(upper_re_obs); // observation parameters (random effects)
  DATA_VECTOR(upper_re_hid); // state process parameters (random effects)
  
  // PARAMETERS (fixed effects first, then random effects)
  PARAMETER_VECTOR(coeff_fe_obs); // observation parameters (fixed effects)
  PARAMETER_VECTOR(log_lambda_obs); // smoothness parameters
  PARAMETER_VECTOR(coeff_fe_hid); // state process parameters (fixed effects)
  PARAMETER_VECTOR(log_lambda_hid); // smoothness parameters
  PARAMETER_VECTOR(log_delta0); // initial distribution
  PARAMETER_VECTOR(coeff_re_obs); // observation parameters (random effects)
  PARAMETER_VECTOR(coeff_re_hid); // state process parameters (random effects)
  
  // Number of observed variables
  int n_var = distcode.size();
  // Number of data rows
  int n = data.rows();
  // Number of time series
  int n_ID;
  if(n_states>1){
    n_ID = log_delta0.size()/(n_states - 1);
  } else n_ID = log_delta0.size();
  
  bool maxInd;
  if(!R_finite(asDouble(kappa))){
    kappa = R_PosInf;
    maxInd = false;
  } else maxInd = true;
  
  //==============================//
  // Transform working parameters //
  //==============================//
  // Observation parameters
  vector<Type> wcoeff_fe_obs(coeff_fe_obs.size());
  for(int p = 0; p < coeff_fe_obs.size(); p++){
    if(R_finite(asDouble(lower_fe_obs(p))) && !R_finite(asDouble(upper_fe_obs(p)))){
      wcoeff_fe_obs(p) = exp(coeff_fe_obs(p))+lower_fe_obs(p);
    } else if(R_finite(asDouble(lower_fe_obs(p))) && R_finite(asDouble(upper_fe_obs(p)))){
      wcoeff_fe_obs(p) = (upper_fe_obs(p)-lower_fe_obs(p)) * invlogit(coeff_fe_obs(p))+lower_fe_obs(p);
    } else if(!R_finite(asDouble(lower_fe_obs(p))) && R_finite(asDouble(upper_fe_obs(p)))){
      wcoeff_fe_obs(p) = -(exp(-coeff_fe_obs(p)) - upper_fe_obs(p));
    } else {
      wcoeff_fe_obs(p) = coeff_fe_obs(p);
    }
  }
  // Transition probabilities
  vector<Type> wcoeff_fe_hid(coeff_fe_hid.size());
  for(int p = 0; p < coeff_fe_hid.size(); p++){
    if(R_finite(asDouble(lower_fe_hid(p))) && !R_finite(asDouble(upper_fe_hid(p)))){
      wcoeff_fe_hid(p) = exp(coeff_fe_hid(p))+lower_fe_hid(p);
    } else if(R_finite(asDouble(lower_fe_hid(p))) && R_finite(asDouble(upper_fe_hid(p)))){
      wcoeff_fe_hid(p) = (upper_fe_hid(p)-lower_fe_hid(p)) * invlogit(coeff_fe_hid(p))+lower_fe_hid(p);
    } else if(!R_finite(asDouble(lower_fe_hid(p))) && R_finite(asDouble(upper_fe_hid(p)))){
      wcoeff_fe_hid(p) = -(exp(-coeff_fe_hid(p)) - upper_fe_hid(p));
    } else {
      wcoeff_fe_hid(p) = coeff_fe_hid(p);
    }
  }  
  
  //======================//
  // Transform parameters //
  //======================//
  // Observation parameters
  vector<Type> par_vec = X_fe_obs * wcoeff_fe_obs + X_re_obs * coeff_re_obs;
  matrix<Type> par_mat(n, par_vec.size()/n);
  for(int i = 0; i < par_mat.cols(); i++) {
    // Matrix with one row for each time step and
    // one column for each parameter
    par_mat.col(i) = par_vec.segment(i*n, n);
  }
  
  // Transition probabilities
  vector<Type> ltpm_vec = X_fe_hid * wcoeff_fe_hid + X_re_hid * coeff_re_hid;
  matrix<Type> ltpm_mat(n, ltpm_vec.size()/n);
  matrix<Type> I = matrix<Type>::Identity(n_states, n_states);
  for(int i = 0; i < ltpm_mat.cols(); i++) {
    // Matrix with one row for each time step and
    // one column for each transition probability
    ltpm_mat.col(i) = ltpm_vec.segment(i*n, n);
  }
  
  // Initial distribution
  matrix<Type> delta0(n_ID, n_states); 
  int nbID = -1;
  // Create array of transition probability matrices
  vector<matrix<Type> > tpm_array(n);
  if(CT==0){
    for(int i = 0; i < n; i++) {
      matrix<Type> tpm(n_states, n_states);
      int cur = 0;
      for (int j = 0; j < n_states; j++) {
        tpm(j, ref_tpm(j) - 1) = 1;
        for (int k = 0; k < n_states; k++) {
          if (k != ref_tpm(j) - 1) {
            tpm(j, k) = exp(ltpm_mat(i, cur));
            cur++;
          }
        }
        tpm.row(j) = tpm.row(j)/tpm.row(j).sum();
      }
      if (statdist == 1) {
        // Case 1: the initial distribution is the stationary distribution
        if(i == 0 || ID(i-1) != ID(i)) {
          nbID++;
          matrix<Type> tpminv = I - tpm; 
          tpminv = (tpminv.array() + 1).matrix(); 
          matrix<Type> ivec(1, n_states); for (int j = 0; j < n_states; ++j) ivec(0, j) = 1;
          // if tpm is ill-conditioned then just use uniform initial distribution 
          try {
            tpminv = tpminv.inverse();
            delta0.row(nbID) = ivec * tpminv;
          } catch(...) {
            for(int j = 0; j < n_states; j++) {
              delta0(nbID, j) = 1.0 / n_states;           
            } 
          }
          //nbID++;
        }
      }
      tpm_array(i) = tpm;
    }
  } else {
    for(int i = 0; i < n; i++) {
      matrix<Type> tpm(n_states, n_states);
      if(i == 0 || ID(i-1) != ID(i)) {
        nbID++;
        tpm_array(i) = I;
        if (statdist == 1){
          for(int j = 0; j < n_states; j++) {
            delta0(nbID, j) = 1.0 / n_states; // just in case an individual has <2 observations          
          } 
        }
      } else {
        tpm.setZero();
        int cur = 0;
        for (int j = 0; j < n_states; j++) {
          for (int k = 0; k < n_states; k++) {
            if (k != ref_tpm(j) - 1) {
              if(j!=k) {
                tpm(j, k) = exp(ltpm_mat(i, cur));
                tpm(j, j) -= tpm(j, k);
              } else {
                tpm(j, k) -= exp(ltpm_mat(i, cur));
              }
              cur++;
            }
          }
          if(j!=(ref_tpm(j)-1)){
            for(int l=0;l<n_states;l++){
              if(l!=(ref_tpm(j)-1)) tpm(j,(ref_tpm(j)-1)) -= tpm(j,l);
            }
          }
          if(maxInd){
            tpm(j,j) = 0.;
            Type rowSums = tpm.row(j).sum();
            for(int l=0;l<n_states;l++){
              if(j!=l){
                tpm(j,l) = kappa * tpm(j,l) / (1.+rowSums); // normalization
                tpm(j,j) -= tpm(j,l);
              }
            }
          }
        }
        if (statdist == 1) {
          if(((nbID==0) & (i==1)) || ((i > 1) && (ID(i-2) != ID(i)))) {
            // Case 1: the initial distribution is the stationary distribution
            matrix<Type> tpminv = -tpm; 
            tpminv = (tpminv.array() + 1).matrix(); 
            matrix<Type> ivec(1, n_states); for (int j = 0; j < n_states; ++j) ivec(0, j) = 1;
            // if tpm is ill-conditioned then just use uniform initial distribution 
            try {
              tpminv = tpminv.inverse();
              delta0.row(nbID) = ivec * tpminv;
            } catch(...) {
              for(int j = 0; j < n_states; j++) {
                delta0(nbID, j) = 1.0 / n_states;           
              } 
            }
            //nbID++;
          }
        }
        tpm *= delta_t(dtIndex(i-1) - 1);
        tpm_array(i) = atomic::expm(tpm);
      }
    }
  }
  
  if (statdist == 0) {
    // Case 2: the initial distribution is estimated
    // transform working parameters
    vector<Type> wlog_delta0(log_delta0.size());
    for(int p = 0; p < log_delta0.size(); p++){
      if(R_finite(asDouble(lower_delta0(p))) && !R_finite(asDouble(upper_delta0(p)))){
        wlog_delta0(p) = exp(log_delta0(p))+lower_delta0(p);
      } else if(R_finite(asDouble(lower_delta0(p))) && R_finite(asDouble(upper_delta0(p)))){
        wlog_delta0(p) = (upper_delta0(p)-lower_delta0(p)) * invlogit(log_delta0(p))+lower_delta0(p);
      } else if(!R_finite(asDouble(lower_delta0(p))) && R_finite(asDouble(upper_delta0(p)))){
        wlog_delta0(p) = -(exp(-log_delta0(p)) - upper_delta0(p));
      } else {
        wlog_delta0(p) = log_delta0(p);
      }
    } 
    delta0.setOnes();
    for(int i = 0; i < n_ID; i++) {
      for(int j = 0; j < n_states - 1; j++) {
        delta0(i, j + 1) = exp(wlog_delta0(j * n_ID + i));
      }
      delta0.row(i) = delta0.row(i)/delta0.row(i).sum();
    }
  }
  
  //===================================//  
  // Compute observation probabilities //
  //===================================//
  // Initialise matrix of probabilities to 1
  matrix<Type> prob(n, n_states);
  for(int i = 0; i < n; i++) {
    if (!R_IsNA(asDouble(known_states(i)))) {
      for (int s = 0; s < n_states; ++s) {
        prob(i, s) = (known_states(i, s) == 1) ? 1 : 0; 
      }
    } else {
      for(int s = 0; s < n_states; s++) {
        prob(i, s) = 1;
      }
    }
  }
  
  // Counter to subset parameter vector
  int par_count = 0;
  int var_count = 0; 
  
  // Loop over observed variables
  for(int var = 0; var < n_var; var++) {
    // Define observation distribution
    std::unique_ptr<Dist<Type>> obsdist = dist_generator<Type>(distcode(var));
    
    // Loop over observations (rows)
    for (int i = 0; i < n; ++i) {
      // Don't update likelihood if the observation is missing
      if(!R_IsNA(asDouble(data(i, var_count)))) {
        // Subset and transform observation parameters
        vector<Type> sub_wpar = 
          par_mat.row(i).segment(par_count, distpar(var) * n_states);
        matrix<Type> par = obsdist->invlink(sub_wpar, n_states);
        // Loop over states (columns)
        for (int s = 0; s < n_states; ++s) {
          // Vector of parameters for state s
          vector<Type> subpar = par.row(s);
          if (datadim(var) > 1) {
            vector<Type> subdat = data.row(i).segment(var_count, datadim(var)); 
            prob(i, s) = prob(i, s) * obsdist->pdf(subdat, subpar, delta_t(i), false);
          } else {
            prob(i, s) = prob(i, s) * obsdist->pdf(data(i, var_count), subpar, delta_t(i), false);
          }
        }        
      }
    }
    var_count = var_count + datadim(var); 
    par_count = par_count + distpar(var) * n_states;
  }
  
  //======================//
  // Priors               //
  //======================//
  Type llk = 0; 
  // fixed effects for observation 
  for (int i = 0; i < coeff_fe_obs.size(); ++i) {
    if (!R_IsNA(asDouble(coeff_fe_obs_prior(i, 0)))) {
      llk += dnorm(coeff_fe_obs(i), coeff_fe_obs_prior(i, 0), coeff_fe_obs_prior(i, 1), 1.0); 
    }
  }
  // fixed effects for hidden  
  for (int i = 0; i < coeff_fe_hid.size(); ++i) {
    if (!R_IsNA(asDouble(coeff_fe_hid_prior(i, 0)))) {
      llk += dnorm(coeff_fe_hid(i), coeff_fe_hid_prior(i, 0), coeff_fe_hid_prior(i, 1), 1.0); 
    }
  }
  // smoothing parameters for observation
  if (ncol_re_obs(0, 0) > -1) {
    for (int i = 0; i < log_lambda_obs.size(); ++i) {
      if (!R_IsNA(asDouble(log_lambda_obs_prior(i, 0)))) {
        llk += dnorm(log_lambda_obs(i), log_lambda_obs_prior(i, 0), log_lambda_obs_prior(i, 1), 1.0); 
      }
    }
  }
  // smoothing parameters for hidden 
  if (ncol_re_hid(0, 0) > -1) {
    for (int i = 0; i < log_lambda_hid.size(); ++i) {
      if (!R_IsNA(asDouble(log_lambda_hid_prior(i, 0)))) {
        llk += dnorm(log_lambda_hid(i), log_lambda_hid_prior(i, 0), log_lambda_hid_prior(i, 1), 1.0); 
      }
    }
  }
  
  //========================//
  // Compute log-likelihood //
  //========================//
  // Initialise log-likelihood
  matrix<Type> phi(delta0.row(0));
  Type sumphi = 0;
  
  // Forward algorithm
  int id = 0;
  for (int i = 0; i < n; ++i) {
    // Re-initialise phi at first observation of each time series
    if(i == 0 || ID(i-1) != ID(i)) {
      phi = delta0.row(id) * tpm_array(i);
      phi = (phi.array() * prob.row(i).array()).matrix();
      id = id + 1;
    } else {
      phi = phi * tpm_array(i);
      phi = (phi.array() * prob.row(i).array()).matrix();
    }
    sumphi = phi.sum();
    llk = llk + log(sumphi);
    phi = phi / sumphi;
  }
  
  // Negative log-likelihood
  Type nllk = -llk;
  
  //===================//
  // Smoothing penalty //
  //===================//
  // Are there smooths in the observation model?
  if((include_smooths > 0) & (ncol_re_obs(0, 0) > -1)) {
    // Index in matrix S
    int S_start = 0;
    
    // Loop over smooths
    for(int i = 0; i < ncol_re_obs.cols(); i++) {
      // Size of penalty matrix for this smooth
      int Sn = ncol_re_obs(1, i) - ncol_re_obs(0, i) + 1;
      
      // Penalty matrix for this smooth
      // (dense matrix for matinvpd and sparse matrix for Quadform)
      matrix<Type> this_S_dense = S_obs.block(S_start, S_start, Sn, Sn);
      Eigen::SparseMatrix<Type> this_S = asSparseMatrix(this_S_dense);
      
      // Coefficients for this smooth
      vector<Type> this_coeff_re = coeff_re_obs.segment(ncol_re_obs(0, i) - 1, Sn);
      
      // Get log-determinant of S^(-1) for additive constant
      Type log_det = 0;
      matrix<Type> inv_this_S = atomic::matinvpd(this_S_dense, log_det);
      log_det = - log_det; // det(S^(-1)) = 1/det(S)
      
      // Add penalty
      nllk = nllk +
        Type(0.5) * Sn * log(2*M_PI) +
        Type(0.5) * log_det -
        Type(0.5) * Sn * log_lambda_obs(i) +
        Type(0.5) * exp(log_lambda_obs(i)) * density::GMRF(this_S).Quadform(this_coeff_re);
      
      // Increase index
      S_start = S_start + Sn;
    }
  }
  
  // Are there smooths in the hidden state model?
  if((include_smooths > 0) & (ncol_re_hid(0, 0) > -1)) {
    // Index in matrix S
    int S_start = 0;
    
    // Loop over smooths
    for(int i = 0; i < ncol_re_hid.cols(); i++) {
      // Size of penalty matrix for this smooth
      int Sn = ncol_re_hid(1, i) - ncol_re_hid(0, i) + 1;
      
      // Penalty matrix for this smooth
      // (dense matrix for matinvpd and sparse matrix for Quadform)
      matrix<Type> this_S_dense = S_hid.block(S_start, S_start, Sn, Sn);
      Eigen::SparseMatrix<Type> this_S = asSparseMatrix(this_S_dense);
      
      // Coefficients for this smooth
      vector<Type> this_coeff_re = coeff_re_hid.segment(ncol_re_hid(0, i) - 1, Sn);
      
      // Get log-determinant of S^(-1) for additive constant
      Type log_det = 0;
      matrix<Type> inv_this_S = atomic::matinvpd(this_S_dense, log_det);
      log_det = - log_det; // det(S^(-1)) = 1/det(S)
      
      // Add penalty
      nllk = nllk +
        Type(0.5) * Sn * log(2*M_PI) +
        Type(0.5) * log_det -
        Type(0.5) * Sn * log_lambda_hid(i) +
        Type(0.5) * exp(log_lambda_hid(i)) * density::GMRF(this_S).Quadform(this_coeff_re);
      
      // Increase index
      S_start = S_start + Sn;
    }
  }
  
  return nllk;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
