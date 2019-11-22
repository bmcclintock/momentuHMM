library(momentuHMM)
library(setRNG)

# set seed
oldRNG<-setRNG::setRNG()
setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=30)

# load nfs data from github
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/nfsData.RData"))

nSims <- 100 # number of imputatons
retryFits <- 30 # number attempt to re-fit based on random perturbation
ncores <- 7 # number of CPU cores

# set time steps based on overlapping time period of location and dive data
predTimes <- seq(max(foragedives$time[1],nfsData$time[1]),min(foragedives$time[nrow(foragedives)],nfsData$time[nrow(nfsData)])+60*60,"hour")

# use default starting values (theta) and explore likelihood surface using retryFits
crwOut<-crawlWrap(nfsData,ncores=ncores,retryFits=20,err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                  fixPar=c(1,1,NA,NA),
                  attempts=100,predTime=predTimes)
plot(crwOut)

# merge crwData object with forage dive data
crwOut <- crawlMerge(crwOut,foragedives,"time")

nbStates <- 3
stateNames <- c("resting", "foraging", "transit")
dist <- list(step = "gamma", angle = "wrpcauchy", dive = "pois")

### construct pseudo-design matrix constraining parameters (to avoid label switching across imputations)
# constrain step length mean parameters: transit > resting
stepDM<-matrix(c(1,1,0,0,0,0,
                 0,0,1,0,0,0,
                 1,0,0,0,0,0,
                 0,0,0,1,0,0,
                 0,0,0,0,1,0,
                 0,0,0,0,0,1),2*nbStates,6,byrow=TRUE,
               dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates)),
                             c("mean_13:(Intercept)","mean_1","mean_2:(Intercept)",
                               paste0("sd_",1:nbStates,":(Intercept)"))))
stepworkBounds <- matrix(c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,Inf,0,rep(Inf,4)),6,2,
                         dimnames=list(colnames(stepDM),c("lower","upper")))

# constrain turning angle concentration parameters: transit > resting
angleDM<-matrix(c(1,1,0,
                  0,0,1,
                  1,0,0),nbStates,3,byrow=TRUE,
                dimnames=list(paste0("concentration_",1:nbStates),
                              c("concentration_13:(Intercept)","concentration_1","concentration_2:(Intercept)")))
angleworkBounds <- matrix(c(-Inf,-Inf,-Inf,Inf,0,Inf),3,2,dimnames=list(colnames(angleDM),c("lower","upper")))

# constrain dive lambda parameters: foraging > transit
diveDM<-matrix(c(1,0,0,
                 0,1,0,
                 0,1,1),nbStates,3,byrow=TRUE,
               dimnames=list(paste0("lambda_",1:nbStates),
                             c("lambda_1:(Intercept)","lambda_23:(Intercept)","lambda_3")))
diveworkBounds <- matrix(c(-Inf,-Inf,-Inf,rep(Inf,2),0),3,2,
                         dimnames=list(colnames(diveDM),c("lower","upper")))

DM<-list(step=stepDM,angle=angleDM,dive=diveDM)
workBounds<-list(step=stepworkBounds,angle=angleworkBounds,dive=diveworkBounds)

Par0 <- getParDM(nbStates = nbStates, dist = dist,
                 Par = Par, DM = DM, workBounds = workBounds,
                 estAngleMean = list(angle = FALSE))

fixPar <- list(dive = c(-100, NA, NA))

# set prior to help prevent working parameters from straying along boundary
prior <- function(par){sum(dnorm(par,0,10,log=TRUE))}

nfsFits <- MIfitHMM(crwOut, nSims = nSims, ncores = ncores, nbStates = nbStates, dist = dist,
                    Par0 = Par0, DM = DM, workBounds = workBounds,
                    estAngleMean = list(angle = FALSE), 
                    fixPar = fixPar, retryFits = retryFits,
                    prior = prior,
                    stateNames=stateNames)
plot(nfsFits,legend.pos="topright")

save.image("nfsExample.RData")

setRNG::setRNG(oldRNG)
