library(momentuHMM)
library(dtwclust)

# load grey seal data from github
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/greySealData_TPM.RData"))

nSims <- 400 # number of imputatons
ncores <- 7 # number of CPU cores

# fit crawl functions
crwOut<-crawlWrap(greySealData,ncores=1,err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                  initial.state=list(a=c(greySealData$x[1],0,greySealData$y[1],0),P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 5000 ^ 2, 10 * 3600 ^ 2))),
                  theta=c(4.269033, -9.230745),fixPar=c(1,1,NA,NA),attempts=100,
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
circularAngleMean<-list(angle=TRUE)

delta<-c(2500,2500,15000) # distance thresholds for activity center states

# specify design matrix formulas
distFormula<- formula(paste("~",paste(paste0("state1(I(Abertay.dist>",delta[1],"))"),
                                      paste0("state2(I(Farne.dist>",delta[2],"))"),
                                      paste0("state3(I(Dogger.dist>",delta[3],"))"),sep = " + ")))
angleFormula<- ~state1(Abertay.angle)+state2(Farne.angle)+state3(Dogger.angle)
stepDM<-list(shape=distFormula,scale=distFormula)
angleDM<-list(mean=angleFormula,concentration=distFormula)

# fix parameters such that activity center states are biased random walks and exploratory states are simple random walks
fixPar<-list(angle=c(rep(100,nCenters),rep(NA,nCenters*2),rep(-100,2)))

# determine bestDat locations within activity center thresholds and specify knownStates accordingly
knownStates<-rep(NA,nrow(bestDat))
for(j in 1:nCenters){
  centerInd<-as.numeric(bestDat[[paste0(rownames(centers),".dist")[j]]]>delta[j])
  knownStates[which(centerInd==0)]<-j
}

bestmod<-fitHMM(bestDat,nbStates=nbStates,dist=dist,Par0=Par0$Par,beta0=Par0$beta,delta0=Par0$delta,fixPar=fixPar,
                estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                DM=list(step=stepDM,angle=angleDM),formula=distFormula,
                stateNames=stateNames,knownStates=knownStates)
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

# specify pseudo-design matrices with step length constraints (to avoid label switching between exploratory states across multiple imputations)
MIstepDM<-matrix(c(1,0,0,0,0,0,0,0,0,0,
                 paste0("I(",paste0("Abertay.dist>",delta[1],")")),0,0,0,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,0,0,0,
                 0,paste0("I(",paste0("Farne.dist>",delta[2],")")),0,0,0,0,0,0,0,0,
                 0,0,1,0,0,0,0,0,0,0,
                 0,0,paste0("I(",paste0("Dogger.dist>",delta[3],")")),0,0,0,0,0,0,0,
                 0,0,0,1,1,0,0,0,0,0,
                 0,0,0,0,1,0,0,0,0,0,
                 0,0,0,0,0,1,0,0,0,0,
                 0,0,0,0,0,paste0("I(",paste0("Abertay.dist>",delta[1],")")),0,0,0,0,
                 0,0,0,0,0,0,1,0,0,0,
                 0,0,0,0,0,0,paste0("I(",paste0("Farne.dist>",delta[2],")")),0,0,0,
                 0,0,0,0,0,0,0,1,0,0,
                 0,0,0,0,0,0,0,paste0("I(",paste0("Dogger.dist>",delta[3],")")),0,0,
                 0,0,0,0,0,0,0,0,1,1,
                 0,0,0,0,0,0,0,0,0,1),nrow=nbStates*2)
colnames(MIstepDM)<-colnames(bestmod$CIbeta$step$est)
colnames(MIstepDM)[c(7,8,15,16)]<-c("shape:(Intercept)","shape_5","scale:(Intercept)","scale_5")

# step length parameter constraints
stepcons<-c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2)

#extract parameter estimates from bestmod and adjust for MIstepDM and stepcons
bestPar<-getPar(bestmod)
bestPar$Par$step[c(7,8,15,16)]<-c(bestPar$Par$step[7],(bestPar$Par$step[8]-bestPar$Par$step[7])^(1/stepcons[8]),bestPar$Par$step[15],(bestPar$Par$step[16]-bestPar$Par$step[15])^(1/stepcons[16]))

startTime<-proc.time()
greySealFits<-MIfitHMM(miData,nSims=nSims,ncores=ncores,poolEstimates=FALSE,nbStates=nbStates,dist=dist,
                       Par0=bestPar$Par,beta0=bestPar$beta0,delta0=bestPar$delta,fixPar=fixPar,
                       estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,
                       DM=list(step=MIstepDM,angle=angleDM),cons=list(step=stepcons),formula=distFormula,
                       stateNames=stateNames,knownStates=knownStates)
endTime<-proc.time()-startTime

greySealPool<-MIpool(greySealFits,ncores=ncores)
plot(greySealPool,plotCI=TRUE,ask=FALSE)

oldRNG<-setRNG::setRNG()
setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=374)#seed=7)
greySealSim<-simData(model=greySealPool,centers=centers,initialPosition = centers[1,],obsPerAnimal = 1515,states=TRUE)
setRNG::setRNG(oldRNG)

save.image("greySealExample_TPM.RData")
save(greySealPool,greySealSim,file="greySealResults_TPM.RData")