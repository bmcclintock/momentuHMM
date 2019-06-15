library(momentuHMM)
library(dtwclust)
library(setRNG)

# set seed
oldRNG<-setRNG::setRNG()
setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=10)

# load grey seal data from github
load(url("https://raw.github.com/bmcclintock/momentuHMM/develop/vignettes/greySealData.RData"))

nSims <- 400 # number of imputatons
ncores <- 7 # number of CPU cores

# fit crawl functions
crwOut<-crawlWrap(greySealData,err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                  fixPar=c(1,1,NA,NA),attempts=100,retryFits=20,
                  timeStep="2 hours")
plot(crwOut,ask=FALSE)

# Find cluster centroids corresponding to activity centers
clusters<-tsclust(greySealData[,c("x","y")], k = 3L, type = "t",
                  seed = 8, trace = TRUE,
                  control = tadpole_control(dc = 2500,
                                            window.size = 1))
centers<-matrix(unlist(clusters@centroids),ncol=2,byrow=TRUE)
centers<-centers[order(centers[,1]),]
rownames(centers)<-c("Abertay","Farne","Dogger")
nCenters<-nrow(centers)

# preprocess data for best predicted track
bestDat<-prepData(crwOut,center=centers)

# Define HMM
nbStates<-5 # number of states
stateNames<-c(rownames(centers),"low","high")

dist<-list(step="weibull",angle="wrpcauchy")
estAngleMean<-list(angle=TRUE)
circularAngleMean<-list(angle=0)

delta<-c(2500,2500,15000) # distance thresholds for activity center states

# specify design matrix formulas
distFormula<- formula(paste("~",paste(paste0("state1(I(Abertay.dist>",delta[1],"))"),
                                      paste0("state2(I(Farne.dist>",delta[2],"))"),
                                      paste0("state3(I(Dogger.dist>",delta[3],"))"),sep = " + ")))
angleFormula<- ~state1(Abertay.angle)+state2(Farne.angle)+state3(Dogger.angle)
stepDM<-list(shape=distFormula,scale=distFormula)
angleDM<-list(mean=angleFormula,concentration=distFormula)

# fix parameters such that exploratory states are simple random walks (by fixing working concentration parameters to large negative number)
# activity center states are biased random walks (because circularAngleMean=list(angle=0)) with a single angular covariate per angleMean, so these are simply fixed to +1 (positive=attraction, negative=repulsion)
fixPar<-list(angle=c(rep(1,nCenters),rep(NA,nCenters*2),rep(-100,nbStates-nCenters)))

# determine bestDat locations within activity center thresholds and specify knownStates accordingly
knownStates<-rep(NA,nrow(bestDat))
for(j in 1:nCenters){
  centerInd<-as.numeric(bestDat[[paste0(rownames(centers),".dist")[j]]]>delta[j])
  knownStates[which(centerInd==0)]<-j
}

## use prior on unconstrained t.p.m. parameters
# note that objects in the fitHMM environment (e.g. nbStates, optInd) can be used in prior function
ln.prior <- function(par) sum(dnorm(par[-c(27+80+1:(nbStates-1))][-optInd[-(1:5)]][-(1:27)],0,10,log=TRUE))

bestmod<-fitHMM(bestDat,nbStates=nbStates,dist=dist,Par0=Par0$Par,beta0=Par0$beta,delta0=Par0$delta,fixPar=fixPar,
                estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                DM=list(step=stepDM,angle=angleDM),formula=distFormula,
                stateNames=stateNames,knownStates=knownStates,prior=ln.prior)
plot(bestmod,plotCI=TRUE,ask=FALSE)

# draw nSims imputed data sets
crwSim<-MIfitHMM(crwOut,ncores=ncores,nSims=nSims,center=centers,fit=FALSE)
miData<-crwSim$miData

# determine miData locations within activity center thresholds and specify knownStates accordingly
knownStates<-list()
for(i in 1:nSims){
  knownStates[[i]]<-rep(NA,nrow(miData[[i]]))
  for(j in 1:nCenters){
    centerInd<-as.numeric(miData[[i]][[paste0(rownames(centers),".dist")[j]]]>delta[j])
    knownStates[[i]][which(centerInd==0)]<-j
  }
}

# specify pseudo-design matrices with step length constraints (to help avoid label switching between exploratory states across multiple imputations)
MIstepDM<-matrix(c(1,paste0("I(",paste0("Abertay.dist>",delta[1],")")),0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,1,paste0("I(",paste0("Farne.dist>",delta[2],")")),  0,0,0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,1,paste0("I(",paste0("Dogger.dist>",delta[3],")")), 0,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,                                                1,0,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,                                                1,1,0,0,0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,1,paste0("I(",paste0("Abertay.dist>",delta[1],")")),0,0,0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,1,paste0("I(",paste0("Farne.dist>",delta[2],")")),  0,0,0,0,
                   0,0,0,0,0,0,0,0,0,0,0,0,1,paste0("I(",paste0("Dogger.dist>",delta[3],")")), 0,0,
                   0,0,0,0,0,0,                                                0,0,0,0,0,0,0,0,1,0,
                   0,0,0,0,0,0,                                                0,0,0,0,0,0,0,0,1,1),nbStates*2,16,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates)),colnames(bestmod$CIbeta$step$est)))
colnames(MIstepDM)[c(7,8,15,16)]<-c("shape_45:(Intercept)","shape_5","scale_45:(Intercept)","scale_5")

# step length parameter constraints
stepworkBounds <- matrix(c(rep(-Inf,7),0,rep(-Inf,7),0,rep(Inf,ncol(MIstepDM))),nrow=ncol(MIstepDM),ncol=2,dimnames=list(colnames(MIstepDM),c("lower","upper")))

#extract parameter estimates from bestmod and adjust for MIstepDM and stepworkBounds
bestPar<-getPar(bestmod)
bestPar$Par$step[c(7,8,15,16)]<-c(bestPar$Par$step[7],log(bestPar$Par$step[8]-bestPar$Par$step[7]),bestPar$Par$step[15],log(bestPar$Par$step[16]-bestPar$Par$step[15]))

startTime<-proc.time()
greySealFits<-MIfitHMM(miData,nSims=nSims,ncores=ncores,poolEstimates=FALSE,nbStates=nbStates,dist=dist,
                       Par0=bestPar$Par,beta0=bestPar$beta,delta0=bestPar$delta,fixPar=fixPar,
                       estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                       DM=list(step=MIstepDM,angle=angleDM),workBounds=list(step=stepworkBounds),formula=distFormula,
                       stateNames=stateNames,knownStates=knownStates,prior=ln.prior)
endTime<-proc.time()-startTime

greySealPool<-MIpool(greySealFits,ncores=ncores)
plot(greySealPool,plotCI=TRUE,ask=FALSE)

setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=2484)
greySealSim<-simData(model=greySealPool,centers=centers,initialPosition = centers[1,],obsPerAnimal = 1515,states=TRUE)
setRNG::setRNG(oldRNG)

save.image("greySealExample.RData")
save(greySealPool,greySealSim,file="greySealResults.RData")
