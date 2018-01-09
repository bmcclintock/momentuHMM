oldRNG<-setRNG::setRNG()

nbStates<-2
stepDist<-"gamma"
angleDist<-"wrpcauchy"
stepPar<-c(100,1000,50,100)
anglePar<-c(0,0,0.25,0.75)
beta <- matrix(-1.5,1,nbStates * (nbStates - 1))

###### The basics
# Instead of individual arguments for step length (stepDist, stepPar, zeroInflation) and turning angle (angleDist, anglePar, angleMean), momentuHMM
# uses list arguments (dist, Par, zeroInflation, estAngleMean, circularAngleMean, etc.) where each element refers to a data stream. At a minimum, 'dist' and 'Par'
# must be specified, where 'dist' determines the names of the data streams (e.g., 'step', 'angle') and the probability distributions.

### 1. moveHMM-type model
setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1)
data1<-simData(nbAnimals=4,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar,angle=anglePar),beta=beta)
mod1<-fitHMM(data1,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par0=list(step=stepPar,angle=anglePar[3:4]))
print(mod1)
plot(mod1,ask=FALSE)

### 2. individual effects on transition probabilities
data2<-simData(nbAnimals=4,nbStates=nbStates,formula=~ID-1,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar,angle=anglePar))
mod2<-fitHMM(data2,nbStates=nbStates,formula=~ID-1,dist=list(step=stepDist,angle=angleDist),Par0=list(step=stepPar,angle=anglePar[3:4]))
print(mod2)
# gamma estimates for individual 1
mod2$CIreal$gamma
# gamma estimates for individual 2
CIreal(mod2,covs=data.frame(ID=2))$gamma
plot(mod2,plotCI=TRUE,ask=FALSE)

### 3. covariate effect on step parameters
wstepPar <- c(log(stepPar[1]),0.1,log(stepPar[2]),-0.1,log(stepPar[3]),0.1,log(stepPar[4]),-0.1) # must be on working (log) scale when specifying DM
stepDM <- list(mean=~cov1,sd=~cov1) # design matrix for step parameters
data3<-simData(nbAnimals=4,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par=list(step=wstepPar,angle=anglePar),beta=beta,DM=list(step=stepDM),nbCovs=1)

#get Par0 from natural scale parameters (based on DM for covariate means)
Par0 <- getParDM(data3,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar,angle=anglePar[3:4]),DM=list(step=stepDM))

mod3<-fitHMM(data3,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par0=Par0,DM=list(step=stepDM))
mod3
plot(mod3,plotCI=TRUE,ask=FALSE)
# specifiy covariate value for plots (step histogram in this case)
plot(mod3,covs=data.frame(cov1=1),plotCI=TRUE,ask=FALSE)

### 4. include additional data stream: omega~beta(shape1,shape2) with zero inflation
omegaDist="beta"
omegaPar=c(1,10,10,1,0.1,0.02)
data4<-simData(nbAnimals=4,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par=list(step=stepPar,angle=anglePar,omega=omegaPar),beta=beta,zeroInflation=list(omega=TRUE))
mod4<-fitHMM(data4,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=list(step=stepPar,angle=anglePar[3:4],omega=omegaPar))
mod4
plot(mod4,ask=FALSE)
plotPR(mod4)

### 5. constrain step mean_1 <= mean_2
stepDM<-matrix(c(1,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4,dimnames=list(NULL,c("mean:(Intercept)","mean_2","sd_1:(Intercept)","sd_2:(Intercept)")))
stepcons <- c(1,2,1,1) # raises "mean_2" working parameter to second power (i.e., constrain it to be positive)
Par0 <- getParDM(data1,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar,angle=anglePar),DM=list(step=stepDM),cons=list(step=stepcons),estAngleMean=list(angle=TRUE))
data5<-simData(nbAnimals=4,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par=Par0,beta=beta,DM=list(step=stepDM),cons=list(step=stepcons))
mod5<-fitHMM(data5,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par0=Par0,DM=list(step=stepDM),cons=list(step=stepcons),estAngleMean=list(angle=TRUE))

# same constraint, different DM
stepDM<-matrix(c(1,1,0,0,-1,0,0,0,0,0,1,0,0,0,0,1),4,4,dimnames=list(NULL,c("mean:(Intercept)","mean_1","sd_1:(Intercept)","sd_2:(Intercept)")))
Par0 <- getParDM(data1,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar,angle=anglePar),DM=list(step=stepDM),cons=list(step=stepcons),estAngleMean=list(angle=TRUE))
mod5b<-fitHMM(data5,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par0=Par0,DM=list(step=stepDM),cons=list(step=stepcons),estAngleMean=list(angle=TRUE))


###### More advanced stuff

### 6. use spatial covariate raster 'forest' loaded with package
setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=1) # movement happens to stay within raster extent with this seed
nbStates<-2
betaForest<-matrix(c(5,-10,-25,50),nrow=2,ncol=2,byrow=TRUE)
spatialCov<-list(forest=forest)
data6<-simData(nbAnimals=4,nbStates=nbStates,formula=~forest,dist=list(step=stepDist,angle=angleDist),Par=list(step=stepPar,angle=anglePar),beta=betaForest,spatialCovs=spatialCov,states=TRUE,retrySims=10)
plotSpatialCov(data6,forest,states=data6$states,ask=FALSE)

mod6<-fitHMM(data6,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),Par0=list(step=stepPar,angle=anglePar[3:4]),beta0=betaForest,formula=~forest)
plot(mod6,plotCI=TRUE,ask=FALSE)

### 7. activity centers
setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=10)
nbStates<-3
stepPar<-c(100,100,1000,50,50,100)
anglePar<-c(0.75,0.75,0.5)

# specify coordinates for 2 activity centers ('haulout' and 'fish') from which 'dist' and 'angle' covariates will be calculated
activityCenters=matrix(c(500,-500,500,-500),2,2,dimnames=list(c("haulout","fish")))

betaInt<-c(-2.5,-2.5,-2.5,-2.5,-0.5,-0.5)
dist1<-c(-1,-1,0,0,-0.05,0)
dist2<-c(0,0,-1,-1,0,-0.05)
dist12<-c(0.1,0,0.1,0,-0.05,-0.05)
betaCenters<-rbind(betaInt,dist1,dist2,dist12)
formula=~log1p(haulout.dist)*log1p(fish.dist)

# step DM
stepDM <- matrix(c(1,0,0,0,0,0,"log1p(haulout.dist)",0,0,0,0,0,0,1,0,0,0,0,0,"log1p(fish.dist)",0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1),6,8)
# step working parameters
wstepPar0 <- log(c(stepPar[1],1.1,stepPar[2],1.1,stepPar[3:6]))

# angle DM (using circular-circular regression on mean angle)
angleDM <- matrix(c("haulout.angle",0,0,0,0,0,0,"fish.angle",0,0,0,0,0,0,0,1,0,0,0,0,0,"log1p(haulout.dist)",0,0,0,0,0,0,1,0,0,0,0,0,"log1p(fish.dist)",0,0,0,0,0,0,1),6,7)
# angle working parameters
wanglePar0 <- c(10,10,boot::logit(c(0.01,anglePar[1],0.01,anglePar[2],anglePar[3])))

centersPar <- list(step=wstepPar0,angle=wanglePar0)
centersDM <- list(step=stepDM,angle=angleDM)

# takes a while to simulate...
data7<-simData(nbAnimals=4,nbStates=nbStates,beta=betaCenters,formula=formula,Par=centersPar,DM=centersDM,circularAngleMean=list(angle=TRUE),centers=activityCenters,dist=list(step="gamma",angle="wrpcauchy"),states=TRUE)
plot(data7,ask=FALSE)

# get Par0
covs<-data.frame(haulout.dist=0,fish.dist=0,haulout.angle=CircStats::circ.mean(data7$haulout.angle),fish.angle=CircStats::circ.mean(data7$fish.angle))
Par0<-getParDM(data=covs,nbStates=nbStates,dist=list(step="gamma",angle="wrpcauchy"),Par=list(step=stepPar,angle=c(atan2(sin(CircStats::circ.mean(data7$haulout.angle)),1+cos(CircStats::circ.mean(data7$haulout.angle))),atan2(sin(CircStats::circ.mean(data7$fish.angle)),1+cos(CircStats::circ.mean(data7$fish.angle))),0,0.01,0.01,0.5)),DM=centersDM,estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE))

# takes a while to fit...
mod7 <- fitHMM(data7,nbStates=nbStates,Par=Par0,DM=list(step=stepDM,angle=angleDM),estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),dist=list(step="gamma",angle="wrpcauchy"))
plot(mod7,plotCI=TRUE,ask=FALSE)

### 8. Same activity centers model, but with temporal irregularity and location measurement error
lambda <- 2 # expect 2 observations per time step based on exponential wait times
errorEllipse <- list(M=150,m=50,r=180) # semi-major, semi-minor, and orientation upper bounds for bivariate normal error ellipse

data8<-simObsData(data7,lambda,errorEllipse) # alternatively, lambda and errorEllipse could simply have been specified in simData()

# fit crawl::crwMLE and predict locations at regular intervals using crawl::crwPredict
err.model <- list(x=~ln.sd.x-1, y=~ln.sd.y-1, rho=~error.corr) # error ellipse model
inits <- list(a=c(0,0,0,0),P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 5000 ^ 2, 10 * 3600 ^ 2))) # initial state
crwOut<-crawlWrap(data8,ncores=1,retryFits=2,err.model=err.model,initial.state=inits,theta=c(4,0),fixPar=c(1,1,NA,NA),attempts=100)
plot(crwOut,ask=FALSE)

# fit best predicted locations
bestData <- prepData(crwOut,centers=activityCenters)
bestFit <- fitHMM(bestData,nbStates=nbStates,Par=Par0,DM=list(step=stepDM,angle=angleDM),estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),dist=list(step="gamma",angle="wrpcauchy"))
bestFit
plot(bestFit,plotCI=TRUE,ask=FALSE)
plotStates(bestFit,ask=FALSE)

# fit 4 realizations of the position process using multiple imputation
miSims <- MIfitHMM(crwOut,nSims=4,ncores=1,fit=FALSE)

miFits <- list()
miFits$HMMfits <- MIfitHMM(miSims$miData,nSims=4,ncores=4,poolEstimates=FALSE,
                           nbStates=nbStates,Par=Par0,DM=list(step=stepDM,angle=angleDM),estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),dist=list(step="gamma",angle="wrpcauchy"),
                           centers=activityCenters,
                           parIS=0,fullPost=FALSE)
withCallingHandlers(miFits$miSum <-  MIpool(miFits$HMMfits,ncores=1),warning=function(w){if(any(grepl("NaNs produced",w)) | any(grepl("Hessian is singular",w))) invokeRestart("muffleWarning")})
miFits <- miHMM(miFits)

# print pooled estimates
miFits
#plot pooled estimates
plot(miFits,plotCI=TRUE,ask=FALSE)
plotPR(miFits)

### 9. Exact same model as 8) but using state-dependent formulas

stepDM2 <- list(mean=~state1(log1p(haulout.dist))+state2(log1p(fish.dist)),sd=~1)

angleDM2 <- list(mean=~state1(haulout.angle)+state2(fish.angle),concentration=~state1(log1p(haulout.dist))+state2(log1p(fish.dist)))

centersDM2 <- list(step=stepDM2,angle=angleDM2)

# fit best predicted locations
bestFit2<-fitHMM(bestData,nbStates=nbStates,Par=Par0,DM=centersDM2,estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),dist=list(step="gamma",angle="wrpcauchy"))

if(!all.equal(bestFit$mod$minimum,bestFit2$mod$minimum)) stop("bestFit and bestFit2 minimum do not match")
if(!all.equal(bestFit$mod$estimate,bestFit2$mod$estimate)) stop("bestFit and bestFit2 estimate do not match")

tmp <- mget(c(paste0("mod",1:7),"miFits","bestFit","bestFit2"))

newExample <- lapply(tmp[c(paste0("mod",1:7),"bestFit","bestFit2")],function(x) x$mod$minimum)
newExample$miFits <- lapply(tmp$miFits$HMMfits,function(x) x$mod$minimum)

#oldExample <- newExample
#save(oldExample,file="extraExamples.RData")

setRNG::setRNG(oldRNG)
