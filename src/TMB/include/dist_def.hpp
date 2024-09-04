// code modified from the hmmTMB package: Michelot T, Glennie R (2023). hmmTMB: Fit Hidden Markov Models using Template Model Builder. R package version 1.0.2, <https://CRAN.R-project.org/package=hmmTMB>.

#ifndef _DIST_
#define _DIST_

// Defines the abstract basic distribution Type and then a number of derived
// distribution types, e.g., Poisson, Normal, and Gamma. 
// Each distribution has a link function, an inverse link function, and a pdf.  

#ifndef _HMMTMB_
#define _HMMTMB_

#endif 

// Abstract Distribution Class 
template <class Type>
class Dist {
public:
  // Constructor
  Dist() {};
  // Link function
  virtual vector<Type> link(const vector<Type>& par, const int& n_states) = 0; 
  // Inverse link function
  virtual matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) = 0;
  
  // Both a univariate and vector input can be given to a pdf() function. 
  // Derived distributions can either have both options or only one. 
  
  // Probability density/mass function
  virtual Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    return(0.0); 
  };
  // Vector input Probability density/mass function
  virtual Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    return(0.0); 
  }
};

// DISCRETE DISTRIBUTIONS ----------------------

template<class Type> 
class Poisson : public Dist<Type> {
public:
  // Constructor
  Poisson() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // rate
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // rate
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type lambda = par(0) * delta_t;
    Type val = dpois(x, lambda, logpdf);
    return(val); 
  }
};

template<class Type> 
class ZeroInflatedPoisson : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedPoisson() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // rate
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i));
    // z
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i))); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // rate
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    // z
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i+ n_states))); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type lambda = par(0) * delta_t;
    Type val = dzipois(x, lambda, par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class ZeroTruncatedPoisson : public Dist<Type> {
public:
  // Constructor
  ZeroTruncatedPoisson() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // rate
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // rate
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type lambda = par(0) * delta_t;
    Type val = dpois(x, lambda) / (1 - dpois(Type(0), lambda));
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class Bernoulli : public Dist<Type> {
public:
  // Constructor
  Bernoulli() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    // prob
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // prob
    for(int i = 0; i < n_states; ++i) par(i, 0) = 1.0 / (1.0 + exp(-wpar(i)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dbinom(x, Type(1), par(0), logpdf);
    return(val); 
  }
};


template<class Type> 
class ZeroInflatedBinomial : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // size
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // prob
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i))); 
    // z
    for (int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // size
    for(int i = 0; i < n_states; ++i) par(i, 0) = wpar(i);
    // prob
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states))); 
    // z
    for(int i = 0; i < n_states; ++i) par(i, 2) = 1.0 / (1.0 + exp(-wpar(i + 2 * n_states))); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val;
    if (x == Type(0)) val = par(2) + (1 - par(2)) * dbinom(x, par(0), par(1)); 
    else val = (1 - par(2)) * dbinom(x, par(0), par(1)); 
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class ZeroInflatedNegativeBinomial : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedNegativeBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // size
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i));
    // prob
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i))); 
    // z
    for (int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // size
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    // prob
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states))); 
    // z
    for(int i = 0; i < n_states; ++i) par(i, 2) = 1.0 / (1.0 + exp(-wpar(i + 2 * n_states))); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val;
    if (x == Type(0)) val = par(2) + (1 - par(2)) * dnbinom(x, par(0), par(1)); 
    else val = (1 - par(2)) * dnbinom(x, par(0), par(1)); 
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class ZeroTruncatedNegativeBinomial : public Dist<Type> {
public:
  // Constructor
  ZeroTruncatedNegativeBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    // size
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i)); 
    // prob
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // size
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    // prob
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dnbinom(x, par(0), par(1)) / (1.0 - dnbinom(Type(0.0), par(0), par(1)));
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class NegativeBinomial : public Dist<Type> {
public:
  // Constructor
  NegativeBinomial() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    // size
    for (int i = 0; i < n_states; ++i) wpar(i) = log(par(i)); 
    // prob
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // size
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    // prob
    for(int i = 0; i < n_states; ++i) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dnbinom(x, par(0), par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class Categorical : public Dist<Type> {
public:
  // Constructor
  Categorical() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    int n_par = par.size() / n_states; 
    matrix<Type> wparmat(n_states, n_par);
    for (int i = 0; i < n_par; ++i) {
      wparmat.col(i) = par.segment(n_states * i, n_states * i + n_par); 
    }
    vector<Type> rowsums = wparmat.rowwise().sum(); 
    vector<Type> wpar(n_states * n_par); 
    for (int j = 1;  j < n_par; ++j) {
      for (int i = 0; i < n_states; ++i) {
        wpar(i + j * n_states) = log(wparmat(i, j) / (1 - rowsums(i))); 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size() / n_states;
    matrix<Type> par(n_states, n_par);
    vector<Type> ewpar = exp(wpar);
    matrix<Type> wparmat(n_states, n_par);
    for (int i = 0; i < n_par; ++i) {
      wparmat.col(i) = ewpar.segment(n_states * i, n_states * i + n_par); 
    }
    vector<Type> etmp = wparmat.rowwise().sum(); 
    for (int i = 0; i < n_states; ++i) {
      Type s = 1.0 / (1.0 + etmp(i)); 
      for (int j = 0; j < n_par; ++j) {
        par(i, j) = exp(wpar(i + n_states * j)) * s;  
      }
    }
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    int obs = int(asDouble(x)); 
    Type val; 
    if (obs == 1) {
      val = 1.0 - par.sum(); 
    } else {
      val = par(obs - 2);
    }
    if (logpdf) val = log(val); 
    return(val); 
  }
};

// CONTINUOUS DISTRIBUTIONS --------------------

template<class Type> 
class Normal : public Dist<Type> {
public:
  // Constructor
  Normal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // sd
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i);
    // sd
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dnorm(x, par(0), par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class TruncatedNormal : public Dist<Type> {
public:
  // Constructor
  TruncatedNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // sd
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    // min
    for (int i = n_states * 2; i < 3 * n_states; ++i) wpar(i) = par(i);
    // max
    for (int i = n_states * 3; i < 4 * n_states; ++i) wpar(i) = par(i);
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i);
    // sd
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    // min
    for (int i = 0; i < n_states; ++i) par(i, 2) = wpar(i + 2 * n_states);
    // max
    for (int i = 0; i < n_states; ++i) par(i, 3) = wpar(i + 3 * n_states);
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type left = pnorm(par(2), par(0), par(1)); 
    Type right = pnorm(par(3), par(0), par(1)); 
    Type val = dnorm(x, par(0), par(1)) / (right - left);
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class FoldedNormal : public Dist<Type> {
public:
  // Constructor
  FoldedNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // sd
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); 
    // sd
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dnorm(x, par(0), par(1)) + dnorm(-x, par(0), par(1));
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class Studentst : public Dist<Type> {
public:
  // Constructor
  Studentst() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // scale
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); ; 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i);
    // scale
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type df = 2 * par(1) * par(1) / (par(1) * par(1) - 1); 
    Type val = dt(x - par(0), df, 0); 
    if (logpdf) val = log(val); 
    return(val); 
  }
};

template<class Type> 
class LogNormal : public Dist<Type> {
public:
  // Constructor
  LogNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // lmean
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // lsd
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // lmean
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i);  
    // lsd
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dnorm(log(x), par(0), par(1), logpdf) / x;
    return(val); 
  }
};

template<class Type> 
class Gamma : public Dist<Type> {
public:
  // Constructor
  Gamma() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // shape and scale
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // shape 
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    // scale
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));  
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dgamma(x, par(0), par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class Gamma2 : public Dist<Type> {
public:
  // Constructor
  Gamma2() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean and sd
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean 
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    // sd 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type scale = par(1) * par(1) / par(0);
    Type shape = par(0) / scale; 
    Type val = dgamma(x, shape, scale, logpdf);
    return(val); 
  }
};

template<class Type> 
class ZeroInflatedGamma : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedGamma() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // shape and scale
    for(int i = 0; i < 2 * n_states; ++i) wpar(i) = log(par(i));
    // z
    for(int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i))); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // shape
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    // scale
    for(int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    // z
    for(int i = 0; i < n_states; ++i) par(i, 2) = 1.0 / (1.0 + exp(-wpar(i + 2*n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = 0.0;
    if(x == Type(0)) {
      val = par(2);
    } else {
      val = (1 - par(2)) * dgamma(x, par(0), par(1), 0);
    }
    if(logpdf) {
      val = log(val);
    }
    return(val); 
  }
};

template<class Type> 
class ZeroInflatedGamma2 : public Dist<Type> {
public:
  // Constructor
  ZeroInflatedGamma2() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean and sd
    for(int i = 0; i < 2 * n_states; ++i) wpar(i) = log(par(i));
    // z
    for(int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i) / (1.0 - par(i))); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for(int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    // sd
    for(int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    // z
    for(int i = 0; i < n_states; ++i) par(i, 2) = 1.0 / (1.0 + exp(-wpar(i + 2*n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = 0.0;
    Type shape = par(0) * par(0) / (par(1) * par(1));
    Type scale = par(1) * par(1) / par(0);
    if(x == Type(0)) {
      val = par(2);
    } else {
      val = (1 - par(2)) * dgamma(x, shape, scale, 0);
    }
    if(logpdf) {
      val = log(val);
    }
    return(val); 
  }
};

template<class Type> 
class Weibull : public Dist<Type> {
public:
  // Constructor
  Weibull() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // shape and scale
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // shape
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    // scale
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dweibull(x, par(0), par(1), logpdf);
    return(val); 
  }
};

template<class Type> 
class Exponential : public Dist<Type> {
public:
  // Constructor
  Exponential() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // rate
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // rate 
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dexp(x, par(0), logpdf);
    return(val); 
  }
};

template<class Type> 
class Beta : public Dist<Type> {
public:
  // Constructor
  Beta() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // shape and scale
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // shape 
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i));
    // scale
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dbeta(x, par(0), par(1), logpdf);
    return(val); 
  }
};

// MIXED DISTRIBUTIONS -------------------------

template<class Type> 
class Tweedie : public Dist<Type> {
public:
  // Constructor
  Tweedie() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean
    for (int i = 0; i < n_states; ++i) wpar(i) = par(i);
    // power
    for (int i = n_states; i < 2 * n_states; ++i) wpar(i) = log(par(i) / (1 - par(i))); 
    // dispersion
    for (int i = 2 * n_states; i < 3 * n_states; ++i) wpar(i) = log(par(i)); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean 
    for (int i = 0; i < n_states; ++i) par(i, 0) = wpar(i); 
    // power
    for (int i = 0; i < n_states; ++i) par(i, 1) = 1 / (1 + exp(-wpar(i + n_states)));
    // dispersion 
    for (int i = 0; i < n_states; ++i) par(i, 2) = exp(wpar(i + 2 * n_states)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dtweedie(x, par(0), par(2), par(1) + 1.0, logpdf);
    return(val); 
  }
};

// ANGULAR DISTRIBUTIONS -----------------------
template<class Type> 
class VonMises : public Dist<Type> {
public:
  // Constructor
  VonMises() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean in (-pi, pi]
    for(int i = 0; i < n_states; i++) wpar(i) = logit((par(i) + M_PI) / (2 * M_PI));
    // concentration > 0
    for(int i = n_states; i < 2*n_states; i++) wpar(i) = log(par(i));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for (int i = 0; i < n_states; i++) par(i, 0) = 2 * M_PI * invlogit(wpar(i)) - M_PI; 
    // concentration
    for (int i = 0; i < n_states; i++) par(i, 1) = exp(wpar(i + n_states));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type b = besselI(Type(par(1)), Type(0));
    Type val = 0;
    if(!logpdf)
      val = 1/(2 * M_PI * b) * exp(par(1) * cos(x - par(0)));
    else
      val = - log(2 * M_PI * b) + par(1) * cos(x - par(0));
    return(val); 
  }
};

template<class Type> 
class WrpCauchy : public Dist<Type> {
public:
  // Constructor
  WrpCauchy() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean in (-pi, pi]
    for(int i = 0; i < n_states; i++) wpar(i) = logit((par(i) + M_PI) / (2 * M_PI));
    // 0 < concentration < 1
    for(int i = n_states; i < 2*n_states; i++) wpar(i) = log(par(i) / (1.0 - par(i)));
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean
    for (int i = 0; i < n_states; i++) par(i, 0) = 2 * M_PI * invlogit(wpar(i)) - M_PI; 
    // concentration
    for (int i = 0; i < n_states; i++) par(i, 1) = 1.0 / (1.0 + exp(-wpar(i + n_states)));
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = (1 - par(1)*par(1)) / (2 * M_PI * (1 + par(1)*par(1) - 2 * par(1) * cos(x - par(0)))); 
    if(logpdf) val = log(val); 
    return(val); 
  }
};

// MULTIVARIATE DISTRIBUTIONS ------------------

template<class Type> 
class MultivariateNormal : public Dist<Type> {
public:
  // Constructor
  MultivariateNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    int n_par = wpar.size()/n_states;
    int dim = this->dim(n_par); 
    int k = 0; 
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = log(par(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        Type tmp = 0.5 * (par(k) + 1);
        wpar(k) = log(tmp / (1 - tmp));
        ++k; 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    int dim = this->dim(n_par);
    int k = 0; 
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d) = wpar(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + dim) = exp(wpar(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 2 * dim) = 1 / (1 + exp(-wpar(k)));
        par(i, d + 2 * dim) = 2 * par(i, d + 2 * dim) - 1;
        ++k; 
      }
    }
    return(par); 
  }
  // Univariate Probability density function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = dnorm(x, par(0), par(1), logpdf); 
    return(val); 
  }
  
  // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    int dim = this->dim(par.size()); 
    vector<Type> y(dim);
    for (int i = 0; i < dim; ++i) y(i) = x(i) - par(i); 
    vector<Type> sds(dim); 
    for (int i = 0; i < dim; ++i) sds(i) = par(i + dim); 
    vector<Type> corr((dim * dim - dim) / 2); 
    for (int i = 0; i < (dim * dim - dim) / 2; ++i) corr(i) = par(i + 2 * dim);
    matrix<Type> Sigma(dim, dim); 
    int k = 0; 
    for (int i = 0; i < dim; ++i) {
      for (int j = i; j < dim; ++j) {
        Sigma(j, i) = sds(j) * sds(i);
        if (i != j) {
          Sigma(j, i) *= corr(k);
          ++k; 
        }
      }
    }
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < i; ++j) {
        Sigma(j, i) = Sigma(i, j); 
      }
    }
    Type val = density::MVNORM(Sigma)(y); 
    val = -val; 
    if (!logpdf) val = exp(val); 
    return(val); 
  }
  
  // Solve for dimension 
  int dim(const double& l) {
    double a = 1; 
    double b = 3; 
    double c = - 2 * l; 
    double d = b * b - 4 * a * c; 
    double root = (-b + sqrt(d)) / (2 * a); 
    return(int(root)); 
  }
};

template<class Type> 
class RWMultivariateNormal : public Dist<Type> {
public:
  // Constructor
  RWMultivariateNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    int n_par = wpar.size()/n_states;
    int dim = this->dim(n_par); 
    int k = 0; 
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = log(par(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        Type tmp = 0.5 * (par(k) + 1);
        wpar(k) = log(tmp / (1 - tmp));
        ++k; 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    int dim = this->dim(n_par);
    int k = 0; 
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d) = wpar(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + dim) = exp(wpar(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 2 * dim) = 1 / (1 + exp(-wpar(k)));
        par(i, d + 2 * dim) = 2 * par(i, d + 2 * dim) - 1;
        ++k; 
      }
    }
    return(par); 
  }
  // Univariate Probability density function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type sdt = par(1) * sqrt(delta_t);
    Type val = dnorm(x, par(0), sdt, logpdf); 
    return(val); 
  }
  
  // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    int dim = this->dim(par.size()); 
    vector<Type> y(dim);
    for (int i = 0; i < dim; ++i) y(i) = x(i) - par(i); 
    vector<Type> sds(dim); 
    for (int i = 0; i < dim; ++i) sds(i) = par(i + dim); 
    vector<Type> corr((dim * dim - dim) / 2); 
    for (int i = 0; i < (dim * dim - dim) / 2; ++i) corr(i) = par(i + 2 * dim);
    matrix<Type> Sigma(dim, dim); 
    int k = 0; 
    for (int i = 0; i < dim; ++i) {
      for (int j = i; j < dim; ++j) {
        Sigma(j, i) = sds(j) * sds(i) * delta_t;
        //Rprintf("sds(j) %f sds(i) %f delta_t %f \n",asDouble(sds(j)),asDouble(sds(i)),asDouble(delta_t));
        if (i != j) {
          Sigma(j, i) *= corr(k);
          ++k; 
        }
      }
    }
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < i; ++j) {
        Sigma(j, i) = Sigma(i, j); 
      }
    }
    Type val = density::MVNORM(Sigma)(y); 
    val = -val; 
    if (!logpdf) val = exp(val); 
    return(val); 
  }
  
  // Solve for dimension 
  int dim(const double& l) {
    double a = 1; 
    double b = 3; 
    double c = - 2 * l; 
    double d = b * b - 4 * a * c; 
    double root = (-b + sqrt(d)) / (2 * a); 
    return(int(root)); 
  }
};

template<class Type> 
class RWcovMultivariateNormal : public Dist<Type> {
public:
  // Constructor
  RWcovMultivariateNormal() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    int n_par = wpar.size()/n_states;
    int dim;
    if(n_par<7) dim = 1;
    else if(n_par>7) dim = 3;
    else dim = 2;
    int k = 0; 
    // tm1
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = log(par(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        Type tmp = 0.5 * (par(k) + 1);
        wpar(k) = log(tmp / (1 - tmp));
        ++k; 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    int dim;
    if(n_par<7) dim = 1;
    else if(n_par>7) dim = 3;
    else dim = 2;
    int k = 0; 
    // tm1
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d) = wpar(k); 
        //Rprintf("tm1 dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d)));
        ++k; 
      }
    }
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + dim) = wpar(k); 
        //Rprintf("mean dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d + dim)));
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 2 * dim) = exp(wpar(k)); 
        //Rprintf("sd dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d + 2 * dim)));
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 3 * dim) = 1 / (1 + exp(-wpar(k)));
        par(i, d + 3 * dim) = 2 * par(i, d + 3 * dim) - 1;
        //Rprintf("corr dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d + 3 * dim)));
        ++k; 
      }
    }
    return(par); 
  }
  // Univariate Probability density function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type mean = par(0)+par(1)*delta_t;
    Type sdt = par(2) * sqrt(delta_t);
    Type val = dnorm(x, mean, sdt, logpdf); 
    return(val); 
  }
  
  // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    int dim;
    int parSize = par.size();
    if(parSize<7) dim = 1;
    else if(parSize>7) dim = 3;
    else dim = 2;
    //int dim = this->dim(par.size()); 
    vector<Type> tm1(dim);
    for (int i = 0; i < dim; ++i) tm1(i) = par(i); 
    vector<Type> sds(dim); 
    for (int i = 0; i < dim; ++i) sds(i) = par(i + 2*dim); 
    vector<Type> corr((dim * dim - dim) / 2); 
    for (int i = 0; i < (dim * dim - dim) / 2; ++i) corr(i) = par(i + 3 * dim);
    matrix<Type> Sigma(dim, dim); 
    int k = 0; 
    for (int i = 0; i < dim; ++i) {
      for (int j = i; j < dim; ++j) {
        Sigma(j, i) = sds(j) * sds(i) * delta_t;
        //Rprintf("sds(j) %f sds(i) %f delta_t %f \n",asDouble(sds(j)),asDouble(sds(i)),asDouble(delta_t));
        if (i != j) {
          Sigma(j, i) *= corr(k);
          ++k; 
        }
      }
    }
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < i; ++j) {
        Sigma(j, i) = Sigma(i, j); 
      }
    }
    
    vector<Type> y(dim);
    for (int i = 0; i < dim; ++i){
      y(i) = x(i) - (tm1(i) + par(i + dim) * delta_t);
      //Rprintf("i %d y %f x %f tm1 %f Sigma(i,i) %f mean(i) %f delta_t %f \n",i,asDouble(y(i)),asDouble(x(i)),asDouble(tm1(i)),asDouble(Sigma(i,i)),asDouble(par(i+dim)),asDouble(delta_t));
    }
    Type val = density::MVNORM(Sigma)(y); 
    //Rprintf("dim %d dens %f \n",dim,asDouble(val));
    val = -val; 
    if (!logpdf) val = exp(val); 
    return(val); 
  }
  
  // Solve for dimension 
  //int dim(const double& l) {
  //  double a = 1; 
  //  double b = 3; 
  //  double c = - 2 * l; 
  //  double d = b * b - 4 * a * c; 
  //  double root = (-b + sqrt(d)) / (2 * a); 
  //  return(int(root)); 
  //}
};

template<class Type> 
class Langevin : public Dist<Type> {
public:
  // Constructor
  Langevin() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    int n_par = wpar.size()/n_states;
    int dim;
    if(n_par<7) dim = 1;
    else if(n_par>7) dim = 3;
    else dim = 2;
    int k = 0; 
    // tm1
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = par(k); 
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        wpar(k) = log(par(k)); 
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        Type tmp = 0.5 * (par(k) + 1);
        wpar(k) = log(tmp / (1 - tmp));
        ++k; 
      }
    }
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    int dim;
    if(n_par<7) dim = 1;
    else if(n_par>7) dim = 3;
    else dim = 2;
    int k = 0; 
    // tm1
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d) = wpar(k); 
        //Rprintf("tm1 dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d)));
        ++k; 
      }
    }
    // means
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + dim) = wpar(k); 
        //Rprintf("mean dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d + dim)));
        ++k; 
      }
    }
    // sds 
    for (int d = 0; d < dim; ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 2 * dim) = exp(wpar(k)); 
        //Rprintf("sd dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d + 2 * dim)));
        ++k; 
      }
    }
    // corr 
    for (int d = 0; d < 0.5*(dim * dim - dim); ++d) {
      for (int i = 0; i < n_states; ++i) {
        par(i, d + 3 * dim) = 1 / (1 + exp(-wpar(k)));
        par(i, d + 3 * dim) = 2 * par(i, d + 3 * dim) - 1;
        //Rprintf("corr dim %d d %d state %d par(i,d) %f \n",dim,d,i,asDouble(par(i,d + 3 * dim)));
        ++k; 
      }
    }
    return(par); 
  }
  // Univariate Probability density function
  Type pdf(const Type& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type mean = par(0)+par(1)*par(2)*par(2)*Type(0.5)*delta_t;
    Type sdt = par(2) * sqrt(delta_t);
    Type val = dnorm(x, mean, sdt, logpdf); 
    return(val); 
  }
  
  // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    int dim;
    int parSize = par.size();
    if(parSize<7) dim = 1;
    else if(parSize>7) dim = 3;
    else dim = 2;
    //int dim = this->dim(par.size()); 
    vector<Type> tm1(dim);
    for (int i = 0; i < dim; ++i) tm1(i) = par(i); 
    vector<Type> sds(dim); 
    for (int i = 0; i < dim; ++i) sds(i) = par(i + 2*dim); 
    vector<Type> corr((dim * dim - dim) / 2); 
    for (int i = 0; i < (dim * dim - dim) / 2; ++i) corr(i) = par(i + 3 * dim);
    matrix<Type> Sigma(dim, dim); 
    int k = 0; 
    for (int i = 0; i < dim; ++i) {
      for (int j = i; j < dim; ++j) {
        Sigma(j, i) = sds(j) * sds(i) * delta_t;
        if (i != j) {
          Sigma(j, i) *= corr(k);
          ++k; 
        }
      }
    }
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < i; ++j) {
        Sigma(j, i) = Sigma(i, j); 
      }
    }
    
    vector<Type> y(dim);
    for (int i = 0; i < dim; ++i){
      y(i) = x(i) - (tm1(i) + Sigma(i,i) * Type(0.5) * par(i + dim));
      //Rprintf("i %d y %f x %f tm1 %f Sigma(i,i) %f mean(i) %f delta_t %f \n",i,asDouble(y(i)),asDouble(x(i)),asDouble(tm1(i)),asDouble(Sigma(i,i)),asDouble(par(i+dim)),asDouble(delta_t));
      for(int j=0; j < dim; ++j){
        if(i!=j) y(i) -= Sigma(i,j) * Type(0.5) * par(j + dim);
        //Rprintf("j %d y %f x %f tm1 %f Sigma(i,j) %f mean(j) %f \n",j,asDouble(y(i)),asDouble(x(i)),asDouble(tm1(i)),asDouble(Sigma(i,j)),asDouble(par(j+dim)));
      }
    }
    Type val = density::MVNORM(Sigma)(y); 
    //Rprintf("dim %d dens %f \n",dim,asDouble(val));
    val = -val; 
    if (!logpdf) val = exp(val); 
    return(val); 
  }
  
  // Solve for dimension 
  //int dim(const double& l) {
  //  double a = 1; 
  //  double b = 3; 
  //  double c = - 2 * l; 
  //  double d = b * b - 4 * a * c; 
  //  double root = (-b + sqrt(d)) / (2 * a); 
  //  return(int(root)); 
  //}
};

template<class Type> 
class Dirichlet : public Dist<Type> {
public:
  // Constructor
  Dirichlet() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size());
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    for (int i = 0; i < n_par; ++i) {
      for (int j = 0; j < n_states; ++j) {
        par(j, i) = exp(wpar(i * n_states + j)); 
      }
    }
    return(par); 
  }
  
  // Multivariate Probability density function 
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type val = 0; 
    for (int i = 0; i < x.size(); ++i) {
      val += (par(i) - 1) * log(x(i)); 
      val -= lgamma(par(i)); 
    }
    val += lgamma(par.sum()); 
    if (!logpdf) val = exp(val); 
    return(val); 
  }
  
};

template<class Type> 
class CRWRice : public Dist<Type> {
public:
  // Constructor
  CRWRice() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean and sd
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean 
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    // sd 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type beta = par(0);
    Type sigma = par(1);
    Type step_tm1 = x(0);
    Type step = x(1);
    Type mu = step_tm1*(Type(1.)-exp(-beta))/beta;
    Type var = sigma*sigma/(beta*beta) * (Type(1.) - Type(2.) / beta * (Type(1.)-exp(-beta))+Type(1.) / (Type(2.)*beta)*(Type(1.)-exp(-Type(2.)*beta)));
    Type sd = sqrt(var);
    Type xabs = fabs(step*mu/var);
    Type b = xabs;
    if(xabs<=Type(709)){ // maximum allowed value for besselI
      b += log(exp(-xabs)*besselI(xabs,Type(0)));
    }
    Type val = log(step) - Type(2.0) * log(sd) +
      (-(step*step+mu*mu)/(Type(2.0)*var)) + b;
    if(!logpdf) val = exp(val);
    return(val); 
  }
};

template<class Type> 
class CRWVonMises : public Dist<Type> {
public:
  // Constructor
  CRWVonMises() {}; 
  // Link function 
  vector<Type> link(const vector<Type>& par, const int& n_states) {
    vector<Type> wpar(par.size()); 
    // mean and sd
    wpar = log(par); 
    return(wpar); 
  } 
  // Inverse link function 
  matrix<Type> invlink(const vector<Type>& wpar, const int& n_states) {
    int n_par = wpar.size()/n_states;
    matrix<Type> par(n_states, n_par);
    // mean 
    for (int i = 0; i < n_states; ++i) par(i, 0) = exp(wpar(i)); 
    // sd 
    for (int i = 0; i < n_states; ++i) par(i, 1) = exp(wpar(i + n_states)); 
    return(par); 
  }
  // Probability density/mass function
  Type pdf(const vector<Type>& x, const vector<Type>& par, const Type& delta_t, const bool& logpdf) {
    Type beta = par(0);
    Type sigma = par(1);
    Type step_tm1 = x(0);
    Type step = x(1);
    Type angle = x(2);
    Type var = sigma*sigma/(beta*beta) * (Type(1.) - Type(2.) / beta * (Type(1.)-exp(-beta))+Type(1.) / (Type(2.)*beta)*(Type(1.)-exp(-Type(2.)*beta)));
    Type kappa = step_tm1*step*exp(-beta)/var;
    Type b = exp(-kappa) * besselI(kappa, Type(0));
    Type val = -log(Type(2.) * M_PI * b) + kappa * (cos(angle)-Type(1.));
    if(!logpdf) val = exp(val);
    return(val); 
  }
};
#endif
