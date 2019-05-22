library(momentuHMM)
library(sp)
library(doRNG)
library(doParallel)
library(setRNG)

# set seed
oldRNG <- setRNG::setRNG()
setRNG::setRNG(kind="L'Ecuyer-CMRG",normal.kind="Inversion",seed=4)

nSims <- 30 # number of imputations
ncores <- 7 # number of cores to use in parallel

timeStep <- 2 # time step (hours)

#################################################################################################################################
## Prepare data
#################################################################################################################################

# load McClintock et al 2013 data
load(url("http://www.esapubs.org/archive/ecol/E094/072/data.RData"))

# back-calculate observation times from j_ti and ssidx
time <- numeric(length(j_ti))
ID <- factor(length(j_ti),levels=1:17)
predTimes <- trunPredTimes <- list()

cumT<-c(0,cumsum(T_n)+1:N)
for(i in 1:N){
  for(t in (cumT[i]+1):(cumT[i+1]-1)){
    if((ssidx[t+1])>=(ssidx[t]+1)){
      for(s in (ssidx[t]+1):(ssidx[t+1])){
        time[s]<- (t-cumT[i]+j_ti[s])*timeStep
        ID[s]<-i
      }
    }
  }
  predTimes[[i]]<-((cumT[i]+1):(cumT[i+1]-1)-cumT[i])*timeStep
  trunPredTimes[[i]]<-predTimes[[i]][predTimes[[i]]>=time[ssidx[(cumT[i]+1)]+1]] # to deal with bug in crawl::crwPredict when any predTimes are before first observed location
}
names(predTimes)<-names(trunPredTimes)<-unique(ID)

# project location data
data<-data.frame(ID=ID,time=time,x=x_ti,y=y_ti)
coordinates(data)<-c("x","y")
proj4string(data)<-CRS("+proj=longlat +ellps=WGS84")
data<-spTransform(data,CRS="+init=epsg:27700 +units=m")

# assume isotropic error ellipse with 50m semi-major axis (see crawl::argosDiag2Cov)
data$ln.sd.x<-rep(3.565449,nrow(data))
data$ln.sd.y<-rep(3.565449,nrow(data))
data$error.corr<-rep(0,nrow(data))

# crawl initial states and parameters
theta<-fixPar<-ln.prior<-list()
for(i in unique(data$ID)){
  theta[[i]]<-c(7, 1)
  fixPar[[i]]<-c(1,1,NA,NA)
  ln.prior[[i]] <- function(theta){-abs(theta[2]+3)/0.5} # laplace prior with location = -3 and scale = 0.5 
}

# fit crawl model to all tracks
crwOut<-crawlWrap(data,ncores=ncores,retryParallel=TRUE,
                  err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                  theta=theta,fixPar=fixPar,prior=ln.prior,attempts=100,
                  predTime=trunPredTimes,retryFits = 10)
plot(crwOut,ask=FALSE)

# replace missing dive activity data with NA
omega_t[which(NAobsvect==1)]<-NA

# create sex covariate
sex<-rep(sexind,times=T_n)
sex[sex==0]<-"M"
sex[sex==1]<-"F"
sex<-factor(sex)

# merge location and dive activity data
omega_sex<-data.frame(ID=factor(rep(1:N,times=T_n),levels=1:N),time=unlist(predTimes),omega=omega_t,sex=sex)
crwOut<-crawlMerge(crwOut,omega_sex,Time.name="time")

# prepare HMM data
hsData<-prepData(crwOut,covNames=c("sex"))

nbStates <- 3
stateNames <- c("resting","foraging","transit")
stepDist<-"weibull"
angleDist<-"wrpcauchy"
omegaDist<-"beta"

#################################################################################################################################
## Fit three models (no effects, sex effects, and individual-level effects) to best predicted tracks using fitHMM
#################################################################################################################################

### Model 1: no individual- or sex-level effects

stepDM<-matrix(c(1,0,0,0,0,0,0,0,0,
                 0,1,0,0,0,0,0,0,0,
                 0,0,1,0,0,0,0,0,0,
                 0,0,0,1,0,0,0,0,0,
                 0,0,0,1,1,0,0,0,0,
                 0,0,0,1,1,1,0,0,0,
                 0,0,0,0,0,0,1,0,1,
                 0,0,0,0,0,0,0,1,1,
                 0,0,0,0,0,0,0,0,1),nrow=3*nbStates,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates),paste0("zeromass_",1:nbStates)),c(paste0("shape_",1:nbStates,":(Intercept)"),"scale:(Intercept)","scale_2","scale_3",paste0("zeromass_",1:nbStates,":(Intercept)"))))
stepworkBounds<-matrix(c(rep(-Inf,4),0,0,rep(-Inf,2),-Inf,rep(Inf,ncol(stepDM))),ncol(stepDM),2,dimnames=list(colnames(stepDM),c("lower","upper")))
stepBounds<-matrix(c(0,5,
                     0,5,
                     0,5,
                     0,14400,
                     0,14400,
                     0,14400,
                     0,1,
                     0,1,
                     0,1),nrow=3*nbStates,byrow=TRUE,dimnames=list(rownames(stepDM),c("lower","upper")))

angleDM<-matrix(c(1,0,0,
                  0,1,1,
                  0,1,0),nrow=nbStates,byrow=TRUE,dimnames=list(paste0("concentration_",1:nbStates),c("concentration_1:(Intercept)","concentration_23:(Intercept)","concentration_2")))
angleBounds<-matrix(c(0,0.95,
                      0,0.95,
                      0,0.95),nrow=nbStates,byrow=TRUE,dimnames=list(rownames(angleDM),c("lower","upper")))
angleworkBounds <- matrix(c(-Inf,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),-Inf,rep(Inf,2),0),ncol(angleDM),2,dimnames=list(colnames(angleDM),c("lower","upper")))

omegaDM<-matrix(c(1,0,0,0,0,0,
                  0,0,1,1,0,0,
                  0,0,1,1,0,0,
                  1,1,0,0,0,0,
                  0,0,1,0,0,0,
                  0,0,1,0,0,0,
                  0,0,0,0,1,1,
                  0,0,0,0,0,1,
                  0,0,0,0,0,1),nrow=nbStates*3,byrow=TRUE,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates)),c("shape_1:(Intercept)","shape2_1","shape_2:(Intercept)","shape1_2","zeromass_1:(Intercept)","zeromass_23:(Intercept)")))
omegaworkBounds <- matrix(c(-Inf,0,-Inf,0,-Inf,-Inf,rep(Inf,ncol(omegaDM))),ncol(omegaDM),2,dimnames=list(colnames(omegaDM),c("lower","upper")))
omegaBounds<-matrix(c(1,10,
                      1,10,
                      1,10,
                      1,10,
                      1,10,
                      1,10,
                      0,1,
                      0,1,
                      0,1),nrow=nbStates*3,byrow=TRUE,dimnames=list(rownames(omegaDM),c("lower","upper")))


stepPar0<-c(0.7,1.1,2.3,700,1200,5500,1.e-3,1.e-6,1.e-100)
anglePar0<-c(1.e-6,1.e-6,0.76)
omegaPar0<-c(1.01,9.99,9.99,3,3.3,3.3,c(0.4,1.e-100,1.e-100))
beta0 <- matrix(c(-2,-2.5,-2,-2.1,-1.1,-1),1)
delta0 <- c(0.99,0.005,0.005)

Par0<-getParDM(data=hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par=list(step=stepPar0,angle=anglePar0,omega=omegaPar0),DM=list(step=stepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=stepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds = list(step=stepBounds,angle=angleBounds,omega=omegaBounds))

fixPar<-list(step=c(rep(NA,nbStates*2),NA,NA,boot::logit(1.e-100)),
             omega=c(rep(NA,4),NA,boot::logit(1.e-100)))

# define prior function to help prevent working parameters from initially straying too far along boundaries
prior <- function(par){
  sum(dnorm(par,0,100,log=TRUE))
}

bestFit<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par0,beta0=beta0,delta0=delta0,DM=list(step=stepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=stepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar,stateNames=stateNames,nlmPar=list(steptol=1.e-8,hessian=FALSE),prior=prior)
# fit without prior constraints
Par<-getPar0(bestFit)
bestFit<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par$Par,beta0=Par$beta,delta0=Par$delta,DM=list(step=stepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=stepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar,stateNames=stateNames,nlmPar=list(steptol=1.e-8,hessian=FALSE))
# double check optimization
Par<-getPar0(bestFit)
bestFit<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par$Par,beta0=Par$beta,delta0=Par$delta,DM=list(step=stepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=stepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar,stateNames=stateNames,optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-9,hessian=FALSE))


### Model 2: sex-level effects
# Note that factor-level covariates must be individually specified (e.g., 'sexF', 'sexM') when using pseudo-design matrix

stepDM.sex<-matrix(c("sexF",0,0,"sexM",0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,"sexF",0,0,"sexM",0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,"sexF",0,0,"sexM",0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,"sexF",0,0,"sexM",0,0,0,0,0,0,0,
                     0,0,0,0,0,0,"sexF","sexF",0,"sexM","sexM",0,0,0,0,0,0,
                     0,0,0,0,0,0,"sexF","sexF","sexF","sexM","sexM","sexM",0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,"sexF",0,1,"sexM",0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,"sexF",1,0,"sexM",
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),nrow=3*nbStates,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates),paste0("zeromass_",1:nbStates)),c("shape_1:sexF","shape_2:sexF","shape_3:sexF","shape_1:sexM","shape_2:sexM","shape_3:sexM","scale_1:sexF","scale_2:sexF","scale_3:sexF","scale_1:sexM","scale_2:sexM","scale_3:sexM","zeromass_1:sexF","zeromass_2:sexF","zeromass_3:(Intercept)","zeromass_1:sexM","zeromass_2:sexM")))
stepworkBounds.sex<-matrix(c(rep(-Inf,7),0,0,-Inf,0,0,rep(-Inf,2),-Inf,rep(-Inf,2),rep(Inf,ncol(stepDM.sex))),ncol(stepDM.sex),2,dimnames=list(colnames(stepDM.sex),c("lower","upper")))
angleDM.sex<-matrix(c("sexF",0,0,"sexM",0,0,
                      0,"sexF","sexF",0,"sexM","sexM",
                      0,"sexF",0,0,"sexM",0),nrow=nbStates,byrow=TRUE,dimnames=list(paste0("concentration_",1:nbStates),c("concentration_1:sexF","concentration_23:sexF","concentration_2:sexF","concentration_1:sexM","concentration_23:sexM","concentration_2:sexM")))
angleworkBounds.sex<-matrix(c(-Inf,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),-Inf,-Inf,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),-Inf,Inf,Inf,0,Inf,Inf,0),ncol(angleDM.sex),2,dimnames=list(colnames(angleDM.sex),c("lower","upper")))

omegaDM.sex<-matrix(c("sexF",0,"sexM",0,0,0,0,0,0,0,0,
                      0,0,0,0,"sexF","sexF","sexM","sexM",0,0,0,
                      0,0,0,0,"sexF","sexF","sexM","sexM",0,0,0,
                      "sexF","sexF","sexM","sexM",0,0,0,0,0,0,0,
                      0,0,0,0,"sexF",0,"sexM",0,0,0,0,
                      0,0,0,0,"sexF",0,"sexM",0,0,0,0,
                      0,0,0,0,0,0,0,0,"sexF",1,"sexM",
                      0,0,0,0,0,0,0,0,0,1,0,
                      0,0,0,0,0,0,0,0,0,1,0),nrow=nbStates*3,byrow=TRUE,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates)),c("shape_1:sexF","shape2_1:sexF","shape_1:sexM","shape2_1:sexM","shape_2:sexF","shape1_2:sexF","shape_2:sexM","shape1_2:sexM","zeromass_1:sexF","zeromass_23:(Intercept)","zeromass_1:sexM")))
omegaworkBounds.sex<-matrix(c(rep(c(-Inf,0),4),-Inf,-Inf,-Inf,rep(Inf,ncol(omegaDM.sex))),ncol(omegaDM.sex),2,dimnames=list(colnames(omegaDM.sex),c("lower","upper")))

fixPar.sex<-list(step=c(rep(NA,nbStates*2*2),NA,NA,boot::logit(1.e-100),NA,NA),
                 omega=c(rep(NA,4*2),NA,boot::logit(1.e-100),NA))

bfPar<-getPar0(bestFit)
Par0.sex<-getPar0(bestFit,formula=~sex,formulaDelta=~sex,DM=list(step=stepDM.sex,angle=angleDM.sex,omega=omegaDM.sex))
Par0.sex$Par$step[c("shape_1:sexF","shape_2:sexF","shape_3:sexF","shape_1:sexM","shape_2:sexM","shape_3:sexM")]<-bfPar$Par$step[c("shape_1:(Intercept)","shape_2:(Intercept)","shape_3:(Intercept)")]
Par0.sex$Par$step[c("scale_1:sexF","scale_2:sexF","scale_3:sexF","scale_1:sexM","scale_2:sexM","scale_3:sexM")]<-bfPar$Par$step[c("scale:(Intercept)","scale_2","scale_3")]
Par0.sex$Par$step[c("zeromass_1:sexF","zeromass_2:sexF","zeromass_1:sexM","zeromass_2:sexM")] <- bfPar$Par$step[c("zeromass_1:(Intercept)","zeromass_2:(Intercept)")]

Par0.sex$Par$angle[c("concentration_1:sexF", "concentration_23:sexF",  "concentration_2:sexF",  "concentration_1:sexM", "concentration_23:sexM",  "concentration_2:sexM" )]<-bfPar$Par$angle
Par0.sex$Par$omega[c("shape_1:sexF", "shape2_1:sexF","shape_2:sexF", "shape1_2:sexF","zeromass_1:sexF")]<-bfPar$Par$omega[c("shape_1:(Intercept)","shape2_1","shape_2:(Intercept)","shape1_2", "zeromass_1:(Intercept)" )]
Par0.sex$Par$omega[c("shape_1:sexM", "shape2_1:sexM","shape_2:sexM", "shape1_2:sexM","zeromass_1:sexM")]<-bfPar$Par$omega[c("shape_1:(Intercept)","shape2_1","shape_2:(Intercept)","shape1_2", "zeromass_1:(Intercept)" )]

bestFit.sex<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par0.sex$Par,beta0=Par0.sex$beta,delta0=Par0.sex$delta,formula=~sex,formulaDelta=~sex,DM=list(step=stepDM.sex,angle=angleDM.sex,omega=omegaDM.sex),workBounds=list(step=stepworkBounds.sex,angle=angleworkBounds.sex,omega=omegaworkBounds.sex),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar.sex,stateNames=stateNames,nlmPar=list(steptol=1.e-8,hessian=FALSE),prior=prior)
# fit without prior constraints
Par<-getPar0(bestFit.sex)
bestFit.sex<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par$Par,beta0=Par$beta,delta0=Par$delta,formula=~sex,formulaDelta=~sex,DM=list(step=stepDM.sex,angle=angleDM.sex,omega=omegaDM.sex),workBounds=list(step=stepworkBounds.sex,angle=angleworkBounds.sex,omega=omegaworkBounds.sex),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar.sex,stateNames=stateNames,nlmPar=list(steptol=1.e-8,hessian=FALSE))
# double check optimization
Par<-getPar0(bestFit.sex)
bestFit.sex<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par$Par,beta0=Par$beta,delta0=Par$delta,formula=~sex,formulaDelta=~sex,DM=list(step=stepDM.sex,angle=angleDM.sex,omega=omegaDM.sex),workBounds=list(step=stepworkBounds.sex,angle=angleworkBounds.sex,omega=omegaworkBounds.sex),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar.sex,stateNames=stateNames,optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-8,hessian=FALSE))

### Model 3: individual-level effects

stepzm<-unique(hsData$ID[which(hsData$step==0)])
stepDM.ind<-matrix(0,nrow=3*nbStates,ncol=N*2*nbStates+length(stepzm)*(nbStates-1)+1,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates),paste0("zeromass_",1:nbStates)),c(paste0("shape_1:(Intercept)ID",1:N),paste0("shape_2:(Intercept)ID",1:N),paste0("shape_3:(Intercept)ID",1:N),paste0("scale:(Intercept)ID",1:N),paste0("scale_2ID",1:N),paste0("scale_3ID",1:N),paste0("zeromass_1:(Intercept)ID",stepzm),paste0("zeromass_2:(Intercept)ID",stepzm),"zeromass_3:(Intercept)")))
stepDM.ind[1,1:N]<-paste0("ID",1:N)
stepDM.ind[2,N+1:N]<-paste0("ID",1:N)
stepDM.ind[3,2*N+1:N]<-paste0("ID",1:N)
stepDM.ind[4,3*N+1:N]<-stepDM.ind[5,3*N+1:N]<-stepDM.ind[6,3*N+1:N]<-paste0("ID",1:N)
stepDM.ind[5,4*N+1:N]<-stepDM.ind[6,4*N+1:N]<-paste0("ID",1:N)
stepDM.ind[6,5*N+1:N]<-paste0("ID",1:N)
stepDM.ind[7,6*N+1:length(stepzm)]<-paste0("ID",stepzm)
stepDM.ind[8,6*N+length(stepzm)+1:length(stepzm)]<-paste0("ID",stepzm)
stepDM.ind[7:9,ncol(stepDM.ind)]<-1

stepworkBounds.ind<-matrix(c(rep(-Inf,3*N),rep(-Inf,N),rep(0,2*N),rep(-Inf,length(stepzm)*(nbStates-1)),-Inf,rep(Inf,3*N),rep(Inf,N),rep(Inf,2*N),rep(Inf,length(stepzm)*(nbStates-1)+1)),ncol=2,dimnames=list(colnames(stepDM.ind),c("lower","upper")))

angleDM.ind<-matrix(0,nrow=nbStates,ncol=N*nbStates,dimnames=list(c(paste0("concentration_",1:nbStates)),c(paste0("concentration_1:(Intercept)ID",1:N),paste0("concentration_23:(Intercept)ID",1:N),paste0("concentration_2ID",1:N))))
angleDM.ind[1,1:N]<-paste0("ID",1:N)
angleDM.ind[2,N+1:N]<-angleDM.ind[3,N+1:N]<-paste0("ID",1:N)
angleDM.ind[2,2*N+1:N]<-paste0("ID",1:N)

angleworkBounds.ind<-matrix(c(rep(c(-Inf,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),-Inf),each=N),rep(Inf,2*N),rep(0,N)),ncol=2,dimnames=list(colnames(angleDM.ind),c("lower","upper")))

omegazm<-unique(hsData$ID[which(hsData$omega==0)])
omegaDM.ind<-matrix(0,nrow=nbStates*3,ncol=N*2*(nbStates-1)+length(omegazm)+1,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates)),c(paste0("shape_1:(Intercept)ID",1:N),paste0("shape2_1ID",1:N),paste0("shape_2:(Intercept)ID",1:N),paste0("shape1_2ID",1:N),paste0("zeromass_1:(Intercept)ID",omegazm),paste0("zeromass_23:(Intercept)"))))
omegaDM.ind[1,1:N]<-omegaDM.ind[4,1:N]<-paste0("ID",1:N)
omegaDM.ind[4,N+1:N]<-paste0("ID",1:N)
omegaDM.ind[2,2*N+1:N]<-omegaDM.ind[3,2*N+1:N]<-omegaDM.ind[5,2*N+1:N]<-omegaDM.ind[6,2*N+1:N]<-paste0("ID",1:N)
omegaDM.ind[2,3*N+1:N]<-omegaDM.ind[3,3*N+1:N]<-paste0("ID",1:N)
omegaDM.ind[7,4*N+1:length(omegazm)]<-paste0("ID",omegazm)
omegaDM.ind[7:9,ncol(omegaDM.ind)]<-1
omegaworkBounds.ind<-matrix(c(rep(c(-Inf,0,-Inf,0),each=N),rep(-Inf,length(omegazm)),-Inf,rep(c(Inf,Inf,Inf,Inf),each=N),rep(Inf,length(omegazm)),Inf),ncol=2,dimnames=list(colnames(omegaDM.ind),c("lower","upper")))

fixPar.ind<-list(step=c(rep(NA,N*nbStates*2),rep(NA,length(stepzm)*(nbStates-1)),boot::logit(1.e-100)),
                 omega=c(rep(NA,N*2*(nbStates-1)),rep(NA,length(omegazm)),boot::logit(1.e-100)))

Par0.ind<-getPar0(bestFit,DM=list(step=stepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),formula=~ID+0,formulaDelta=~ID+0)

# fit each track individually to get initial values for full model
registerDoParallel(cores=ncores) 
bestFit.all<-foreach(i=1:N) %dorng% {
  data.ind<-subset(hsData,ID==i)
  tmpPar1 <- Par0
  tmpPar2 <- getPar0(bestFit)
  if(!all(data.ind$step>0,na.rm=TRUE)){
    tmpstepDM<-stepDM
    tmpstepworkBounds<-stepworkBounds
    tmpstepBounds<-stepBounds
    tmpfixPar<-fixPar
  } else {
    tmpstepDM<-stepDM[1:(2*nbStates),1:(2*nbStates)]
    tmpstepworkBounds<-stepworkBounds[1:(2*nbStates),]
    tmpstepBounds<-stepBounds[1:(2*nbStates),]
    tmpPar1$step<-tmpPar1$step[1:(2*nbStates)]
    tmpPar2$Par$step<-tmpPar2$Par$step[1:(2*nbStates)]
    tmpfixPar<-fixPar
    tmpfixPar$step<-NULL
  }
  tryFits <- list()
  tryFits[[1]]<-tryCatch(fitHMM(data.ind,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar1,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,nlmPar=list(steptol=1.e-9,hessian=FALSE),prior=prior),error=function(e) e)
  tryFits[[2]]<-tryCatch(fitHMM(data.ind,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar2$Par,beta0=tmpPar2$beta,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,nlmPar=list(steptol=1.e-9,hessian=FALSE),prior=prior),error=function(e) e)
  
  if(inherits(tryFits[[1]],"error") & inherits(tryFits[[2]],"error")) stop("both sets of starting values failed for individual ",i)
  
  if(inherits(tryFits[[1]],"error")){
    fit <- tryFits[[2]]
  } else if(inherits(tryFits[[2]],"error")){
    fit <- tryFits[[1]]
  } else fit <- tryFits[[which.min(c(tryFits[[1]]$mod$minimum,tryFits[[2]]$mod$minimum))]]
  
  # fit without prior constraints
  tmpPar0 <- getPar0(fit)
  if(any(tmpPar0$delta==1)){
    del0 <- rep(1.e-100,nbStates)
    del0[which(tmpPar0$delta==1)] <- 1-1.e-100*(nbStates-1)
  } else del0 <- tmpPar0$delta
  fit <- fitHMM(data.ind,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar0$Par,beta0=tmpPar0$beta,delta0=del0,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,nlmPar=list(steptol=1.e-9,hessian=FALSE))
  
  # double check optimization with "Nelder-Mead" method
  tmpPar0 <- getPar0(fit)
  if(any(tmpPar0$delta==1)){
    del0 <- rep(1.e-100,nbStates)
    del0[which(tmpPar0$delta==1)] <- 1-1.e-100*(nbStates-1)
  } else del0 <- tmpPar0$delta
  fitHMM(data.ind,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar0$Par,beta0=tmpPar0$beta,delta0=del0,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-9,hessian=FALSE))
  
}
stopImplicitCluster()


#specify initial values for full model based on individual model fits (note this is identical to fitting each track individually but is much slower to fit and calculate hessian)
for(i in 1:N){
  bfPar<-getPar0(bestFit.all[[i]])
  Par0.ind$Par$step[match(paste0(colnames(stepDM),"ID",i),names(Par0.ind$Par$step),nomatch=0)]<-bfPar$Par$step[match(names(Par0.ind$Par$step),paste0(colnames(stepDM),"ID",i),nomatch=0)]
  Par0.ind$Par$angle[match(paste0(colnames(angleDM),"ID",i),names(Par0.ind$Par$angle),nomatch=0)]<-bfPar$Par$angle[match(names(Par0.ind$Par$angle),paste0(colnames(angleDM),"ID",i),nomatch=0)]
  Par0.ind$Par$omega[match(paste0(colnames(omegaDM),"ID",i),names(Par0.ind$Par$omega),nomatch=0)]<-bfPar$Par$omega[match(names(Par0.ind$Par$omega),paste0(colnames(omegaDM),"ID",i),nomatch=0)]
  Par0.ind$beta[paste0("ID",i),]<-bfPar$beta
  Par0.ind$delta[paste0("ID",i),]<-bestFit.all[[i]]$CIbeta$delta$est
}  
bestFit.ind<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),formula=~ID+0,formulaDelta=~ID+0,Par0=Par0.ind$Par,beta0=Par0.ind$beta,delta0=Par0.ind$delta,DM=list(step=stepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),workBounds=list(step=stepworkBounds.ind,angle=angleworkBounds.ind,omega=omegaworkBounds.ind),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),fixPar=fixPar.ind,stateNames=stateNames,optMethod="Nelder-Mead",control=list(maxit=100000,hessian=FALSE,abstol=1.e-9))

#################################################################################################################################
## compare models using AIC
#################################################################################################################################
AIC(bestFit,bestFit.sex,bestFit.ind)

# plot individual-level model
plot(bestFit.ind,ask=FALSE)

#################################################################################################################################
## Fit individual-level effects model using MIfitHMM
#################################################################################################################################

# draw realizations of position process from crwOut
miData <- MIfitHMM(crwOut,nSims=nSims,ncores=ncores,covNames=c("sex"),fit=FALSE)

# fit MIfitHMM to each track individually
registerDoParallel(cores=ncores) 
miBestFit.all<-foreach(i=1:N) %dorng% {
  data.ind<- lapply(miData$miData,function(x) subset(x,ID==i))
  
  tmpPar1 <- Par0
  tmpPar2 <- getPar0(bestFit.all[[i]])
  if(!all(data.ind$step>0,na.rm=TRUE)){
    tmpstepDM<-stepDM
    tmpstepworkBounds<-stepworkBounds
    tmpstepBounds<-stepBounds
  } else {
    tmpstepDM<-stepDM[1:(2*nbStates),1:(2*nbStates)]
    tmpstepworkBounds<-stepworkBounds[1:(2*nbStates),]
    tmpstepBounds<-stepBounds[1:(2*nbStates),]
    tmpPar1$step<-tmpPar1$step[1:(2*nbStates)]
    tmpPar2$Par$step<-tmpPar2$Par$step[1:(2*nbStates)]
  }
  tmpfixPar<-fixPar
  tmpfixPar$step<-NULL
  
  tryFits <- list()
  tryFits[[1]]<-MIfitHMM(data.ind,poolEstimates=FALSE,
                         nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar1,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                         retryFits=5,nlmPar=list(steptol=1.e-9,hessian=FALSE),prior=prior)
  tryFits[[2]]<-MIfitHMM(data.ind,poolEstimates=FALSE,
                         nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar2$Par,beta0=tmpPar2$beta,delta0=tmpPar2$delta,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                         retryFits=5,nlmPar=list(steptol=1.e-9,hessian=FALSE),prior=prior)
  
  fit <- list()
  for(j in 1:nSims){
    if(inherits(tryFits[[1]][[j]],"error") & inherits(tryFits[[2]][[j]],"error")) {
      tryFits[[1]][[j]]<-fitHMM(data.ind[[j]],
                                nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar1,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                                optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-9,hessian=FALSE),prior=prior)
      if(!inherits(tryFits[[1]][[j]],"error")){
        tmpPar1<-getPar0(tryFits[[1]][[j]])
        tryFits[[1]][[j]]<-fitHMM(data.ind[[j]],
                                  nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar1$Par,beta0=tmpPar1$beta,delta0=tmpPar1$delta,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                                  retryFits=5,nlmPar=list(steptol=1.e-9,hessian=FALSE),prior=prior)
      }
      tryFits[[2]][[j]]<-fitHMM(data.ind[[j]],
                                nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar2$Par,beta0=tmpPar2$beta,delta0=tmpPar2$delta,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                                optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-9,hessian=FALSE),prior=prior)
      if(!inherits(tryFits[[2]][[j]],"error")){
        tmpPar2<-getPar0(tryFits[[2]][[j]])
        tryFits[[2]][[j]]<-fitHMM(data.ind[[j]],
                                  nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar2$Par,beta0=tmpPar2$beta,delta0=tmpPar2$delta,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                                  retryFits=5,nlmPar=list(steptol=1.e-9,hessian=FALSE),prior=prior)
      }
    }  

    if(inherits(tryFits[[1]][[j]],"error") & inherits(tryFits[[2]][[j]],"error")) stop("both sets of starting values failed for individual ",i)
    
    if(inherits(tryFits[[1]][[j]],"error")){
      fit[[j]] <- tryFits[[2]][[j]]
    } else if(inherits(tryFits[[2]][[j]],"error")){
      fit[[j]] <- tryFits[[1]][[j]]
    } else fit[[j]] <- tryFits[[which.min(c(tryFits[[1]][[j]]$mod$minimum,tryFits[[2]][[j]]$mod$minimum))]][[j]]
  }
  
  # fit without prior constraints
  tmpPar0 <- lapply(fit,function(x) {
    tmp<-getPar0(x);
    if(any(tmp$delta==1)){
      del0 <- rep(1.e-100,nbStates)
      del0[which(tmp$delta==1)] <- 1-1.e-100*(nbStates-1)
    } else del0 <- tmp$delta
    tmp})
  fit <- MIfitHMM(data.ind,poolEstimates=FALSE,
                  nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=lapply(tmpPar0,function(x) x$Par),beta0=lapply(tmpPar0,function(x) x$beta),delta0=lapply(tmpPar0,function(x) x$delta),DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
                  nlmPar=list(steptol=1.e-9,hessian=FALSE))  
  
  # double check optimization with "Nelder-Mead" method
  tmpPar0 <- lapply(fit,function(x) {
    tmp<-getPar0(x);
    if(any(tmp$delta==1)){
      del0 <- rep(1.e-100,nbStates)
      del0[which(tmp$delta==1)] <- 1-1.e-100*(nbStates-1)
    } else del0 <- tmp$delta
    tmp})
  
  MIfitHMM(data.ind,poolEstimates=FALSE,
           nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=lapply(tmpPar0,function(x) x$Par),beta0=lapply(tmpPar0,function(x) x$beta),delta0=lapply(tmpPar0,function(x) x$delta),DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),workBounds=list(step=tmpstepworkBounds,angle=angleworkBounds,omega=omegaworkBounds),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar,stateNames=stateNames,
           optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-9))
  
}
stopImplicitCluster()

registerDoParallel(cores=ncores) 
miSum.all<-foreach(i=1:N) %dorng% {
  MIpool(miBestFit.all[[i]])
}
stopImplicitCluster()

#specify initial values for full model based on individual model fits
tmpPar0.ind <- list()
for(j in 1:nSims){
  tmpPar0.ind[[j]] <- getPar0(bestFit.ind,DM=list(step=stepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),formula=~ID+0,formulaDelt=~ID+0)
  for(i in 1:N){
    bfPar <- getPar0(miBestFit.all[[i]][[j]])
    tmpPar0.ind[[j]]$Par$step[match(paste0(colnames(stepDM),"ID",i),names(Par0.ind$Par$step),nomatch=0)]<-bfPar$Par$step[match(names(Par0.ind$Par$step),paste0(colnames(stepDM),"ID",i),nomatch=0)]
    tmpPar0.ind[[j]]$Par$angle[match(paste0(colnames(angleDM),"ID",i),names(Par0.ind$Par$angle),nomatch=0)]<-bfPar$Par$angle[match(names(Par0.ind$Par$angle),paste0(colnames(angleDM),"ID",i),nomatch=0)]
    tmpPar0.ind[[j]]$Par$omega[match(paste0(colnames(omegaDM),"ID",i),names(Par0.ind$Par$omega),nomatch=0)]<-bfPar$Par$omega[match(names(Par0.ind$Par$omega),paste0(colnames(omegaDM),"ID",i),nomatch=0)]
    tmpPar0.ind[[j]]$beta[paste0("ID",i),]<-bfPar$beta
    tmpPar0.ind[[j]]$delta[paste0("ID",i),]<-miBestFit.all[[i]][[j]]$CIbeta$delta$est
  }  
  tmpPar0.ind[[j]]$Par$step<-tmpPar0.ind[[j]]$Par$step[1:(N*2*nbStates)]
}
tmpstepDM.ind<-stepDM.ind[1:(2*nbStates),1:(N*2*nbStates)]
tmpstepworkBounds.ind<-stepworkBounds.ind[1:(N*2*nbStates),]
tmpstepBounds<-stepBounds[1:(2*nbStates),]
tmpfixPar.ind<-fixPar.ind
tmpfixPar.ind$step<-NULL

#Fit multiple imputations (warning -- this takes a looooooong time!)
miBestFit.ind<-MIfitHMM(miData$miData,nSims=nSims,ncores=ncores,poolEstimates=FALSE,
                        nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),formula=~ID+0,formulaDelta=~ID+0,Par0=lapply(tmpPar0.ind,function(x) x$Par),beta0=lapply(tmpPar0.ind,function(x) x$beta),delta0=lapply(tmpPar0.ind,function(x) x$delta),DM=list(step=tmpstepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),workBounds=list(step=tmpstepworkBounds.ind,angle=angleworkBounds.ind,omega=omegaworkBounds.ind),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),fixPar=tmpfixPar.ind,stateNames=stateNames,
                        optMethod="Nelder-Mead",control=list(maxit=100000,abstol=1.e-9))

miSum.ind<-MIpool(miBestFit.ind,ncores=ncores)

plot(miSum.ind,plotCI=TRUE,ask=FALSE)

hsActivityBudgets<-timeInStates(miBestFit.ind,by="sex",ncores=ncores)
hsActivityBudgets

save.image("harbourSealExample.RData")

setRNG::setRNG(oldRNG)
