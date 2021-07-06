# Fit HMM to VDBA striated caracara (Phalcoboenus australis), Fahlbusch & Harrington (2019; https://doi.org/10.1242/jeb.211136) using 
# momentuHMM. This script replicates the Supplementary Tutorial example in McClintock et al. (2020; https://doi.org/10.1111/ele.13610)
library(momentuHMM)                                                  

# Data --------------------------------------------------------------------
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/caracaraExample.RData"))
# prepare data for momentuHMM
vdba <- momentuHMM::prepData(data=as.data.frame(VDBA), coordNames=NULL) 
# summarize data 
summary(vdba, dataNames="VDBA")
# plot data
plot(vdba, dataNames="VDBA", ask=FALSE) 

# Fit Model ---------------------------------------------------------------
# set number of states 
N <- 4 
# name states 
stateNames <- c("resting", "activity1", "activity2", "flying") 
# gamma distribution for VDBA
dist <- list(VDBA = "gamma")                                            
# stationary = TRUE means use initial distribution to be equilibrium distribution
stationary <- TRUE                                                      

## starting values for optimization
# gamma means for states 1,...,N
mu <- c(0.005,0.01,0.05,0.1)                                     
# gamma standard deviations for states 1,...,N
sig <- c(0.01,0.02,0.03,0.06)               
Par0 <- list(VDBA = c(mu,sig))
# logit-scale off-diagonal elements of the state transition probability matrix
beta0 <- matrix(-1.5, nrow=1, ncol=N*(N-1)) # default values                              

## fit the 4-state HMM
fit <- momentuHMM::fitHMM(data = vdba,
                          nbStates = N,
                          dist = dist,
                          Par0 = Par0,
                          beta0 = beta0,
                          stationary = stationary,
                          stateNames = stateNames,
                          nlmPar = list(print.level=2))

# look at fit 
fit

# plot state-dependent distributions 
plot(fit, ask=FALSE)

# Model Checking ----------------------------------------------------------
# pseudo-residuals
momentuHMM::plotPR(fit, lag.max = 30)

## simulate from fitted model
set.seed(5282)
nsims <- 100
simobs <- matrix(0, nr = length(VDBA), nc = nsims)
for (i in 1:nsims) {
  cat(i, " / ", nsims, "\r")
  isim <- suppressMessages(momentuHMM::simData(model=fit, obsPerAnimal=length(VDBA)+1))
  simobs[,i] <- isim$VDBA[1:length(VDBA)]
}

quans <- seq(0, 0.6, 0.01)
qmu <- qlcl <- qucl <- rep(0, length(quans))
for (i in 1:length(quans)) {
  qsamp <- apply(simobs, 2, FUN = function(x) {quantile(x, prob = quans[i])})
  qmu[i] <- mean(qsamp)
  qlcl[i] <- quantile(qsamp, prob = 0.025)
  qucl[i] <- quantile(qsamp, prob = 0.975)
}
qobs <- quantile(VDBA, prob = quans)
plot(quans, qobs, type = "l", lwd = 1.5, xlab = "Probabiliy", ylab = "Quantile")
lines(quans, qmu, col = "steelblue", lwd = 1.5)
lines(quans, qlcl, col = "steelblue", lwd = 1.5, lty = "dashed")
lines(quans, qucl, col = "steelblue", lwd = 1.5, lty = "dashed")

maxlag <- 30
acfs <- apply(simobs, 2, FUN = function(x) {acf(x, lag.max = maxlag, plot = FALSE)$acf})
acfmu <- rowMeans(acfs)
acfci <- apply(acfs, 1, quantile, prob = c(0.025, 0.975))
acf(VDBA, lag.max = maxlag)
points(0:maxlag, acfmu, col = "blue", pch = 19)
points(0:maxlag, acfci[1,], col = "blue", pch = "-", cex = 2)
points(0:maxlag, acfci[2,], col = "blue", pch = "-", cex = 2)


# State Decoding ----------------------------------------------------------
# global state decoding
gstates <- momentuHMM::viterbi(fit)  
# local state decoding
stProbs <- momentuHMM::stateProbs(fit) 
lstates <- apply(stProbs, 1, which.max)

momentuHMM::plotStates(fit, ask=FALSE)

# compare global and local state decoding
cb <- c("#F0E442","#E69F00", "#56B4E9", "#009E73")
plot(VDBA, col = cb[gstates], pch = 19)
plot(VDBA, col = cb[lstates], pch = 19)
table(gstates, lstates)


# Appendix ----------------------------------------------------------------

### fitting alternative models 
## three states
stateNames <- c("resting", "activity", "flying")                        # arbitrary names for the states

### starting values for optimization
mu <- c(0.005, 0.02, 0.04)                                              # gamma means for states 1,...,N
sig <- c(0.01, 0.03, 0.06)                                              # gamma standard deviations for states 1,...,N
Par0 <- list(VDBA = c(mu,sig))

### fit the 3-state HMM
fit3 <- momentuHMM::fitHMM(data = vdba,
                          nbStates = 3,
                          dist = dist,
                          Par0 = Par0,
                          stationary = stationary,
                          stateNames = stateNames,
                          nlmPar = list(print.level=2))

momentuHMM::plotPR(fit3) # pseudo-residuals

## two states
### starting values for optimization
mu <- c(0.01, 0.04)                                                     # gamma means for states 1,...,N
sig <- c(0.02, 0.06)                                                    # gamma standard deviations for states 1,...,N
Par0 <- list(VDBA = c(mu,sig))

### fit the 2-state HMM
fit2 <- momentuHMM::fitHMM(data = vdba,
                           nbStates = 2,
                           dist = dist,
                           Par0 = Par0,
                           stationary = stationary,
                           nlmPar = list(print.level=2))

momentuHMM::plotPR(fit2) # pseudo-residuals

# AIC model selection
AIC(fit2,fit3,fit)
