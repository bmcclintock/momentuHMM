library(momentuHMM)

# load nfs data from github
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/nfsData.RData"))

nSims <- 100 # number of imputatons
retryFits <- 30 # number attempt to re-fit based on random perturbation
ncores <- 7 # number of CPU cores

predTimes <- seq(as.POSIXct("2007-10-07 17:49:25",tz="GMT"),as.POSIXlt("2007-10-17 04:49:25",tz="GMT"),"hour")

inits<-list(a=c(nfsData$x[1],0,nfsData$y[1],0),P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 5000 ^ 2, 10 * 3600 ^ 2)))

crwOut<-crawlWrap(nfsData,ncores=ncores,retryFits=retryFits,err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                  initial.state=inits,
                  theta=c(4.422616, -7.938689),fixPar=c(1,1,NA,NA),
                  attempts=100,
                  predTime=predTimes)
plot(crwOut)

crwOut <- crawlMerge(crwOut,foragedives,"time")

nbStates <- 3
stateNames <- c("resting", "foraging", "transit")
dist <- list(step = "gamma", angle = "wrpcauchy", dive = "pois")

### construct pseudo-design matrix constraining parameters (to avoid label switching across imputations)
# constrain step length mean parameters: transit > resting
stepDM<-matrix(c(1,0,1,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1),nrow=2*nbStates)
stepcons<-c(1,2,1,1,1,1)
# constrain turning angle concentration parameters: transit > resting
angleDM<-matrix(c(1,0,1,-1,0,0,0,1,0),nrow=nbStates)
anglecons<-c(1,2,1)
# constrain dive lambda parameters: foraging > transit
diveDM<-matrix(c(1,0,0,0,1,1,0,0,-1),nrow=nbStates)
divecons<-c(1,1,2)

DM<-list(step=stepDM,angle=angleDM,dive=diveDM)
cons<-list(step=stepcons,angle=anglecons,dive=divecons)

Par0 <- getParDM(nbStates = nbStates, dist = dist,
                 Par = Par, DM = DM, cons = cons,
                 estAngleMean = list(angle = FALSE))

fixPar <- list(dive = c(-100, NA, NA))

nfsFits <- MIfitHMM(crwOut, nSims = nSims, ncores = ncores, nbStates = nbStates, dist = dist,
                    Par0 = Par0, DM = DM, cons = cons,
                    estAngleMean = list(angle = FALSE), 
                    fixPar = fixPar, retryFits = retryFits,
                    stateNames=stateNames)
plot(nfsFits)

save.image("nfsExample.RData")
