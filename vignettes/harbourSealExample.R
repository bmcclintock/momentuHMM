library(momentuHMM)
library(sp)
library(doParallel)

nSims<-30 # number of imputations
ncores<-7 # number of cores to use in parallel

timeStep <- 2 # time step (hours)

#################################################################################################################################
## Prepare data
#################################################################################################################################

# load McClintock et al 2013 data
load(url("http://www.esapubs.org/archive/ecol/E094/072/data.RData"))

# back-calculate observation times from j_ti and ssidx
time<-numeric(length(j_ti))
ID<-factor(length(j_ti),levels=1:17)
predTimes<-list()

cumT<-c(0,cumsum(T_n)+1:N)
for(i in 1:N){
  for(t in (cumT[i]+1):(cumT[i+1]-1)){
    if((ssidx[t+1])>=(ssidx[t]+1)){
      for(s in (ssidx[t]+1):(ssidx[t+1])){
        time[s]<- (t-cumT[i]+j_ti[s])*timeStep*60*60
        ID[s]<-i
      }
    }
  }
  predTimes[[i]]<-((cumT[i]+1):(cumT[i+1]-1)-cumT[i])*timeStep*60*60
}
names(predTimes)<-unique(ID)

# project location data
data<-data.frame(ID=ID,time=time,x=x_ti,y=y_ti)
coordinates(data)<-c("x","y")
proj4string(data)<-CRS("+proj=longlat +ellps=WGS84")
data<-spTransform(data,CRS="+init=epsg:27700 +units=m")

# assume isotropic error ellipse with 50m semi-major axis (see crawl::argosDiag2Cov)
data$ln.sd.x<-rep(3.565449,nrow(data))
data$ln.sd.y<-rep(3.565449,nrow(data))
data$error.corr<-rep(0,nrow(data))

# fit crawl model to all tracks
init<-theta<-fixPar<-list()
for(i in unique(data$ID)){
  init[[i]]<-list(a=c(subset(data,ID==i)$x[1],0,subset(data,ID==i)$y[1],0),P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 5000 ^ 2, 10 * 3600 ^ 2)))
  theta[[i]]<-c(3.454, -7.945)
  fixPar[[i]]<-c(1,1,NA,NA)
}
crwOut<-crawlWrap(data,ncores=ncores,err.model=list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                  initial.state=init,
                  theta=theta,fixPar=fixPar,attempts=100,
                  predTime=predTimes,retryFits = 10)
plot(crwOut,ask=FALSE)

# replace missing dive activity data with NA
omega_t[which(NAobsvect==1)]<-NA

# create sex covariate
sex<-rep(sexind,times=T_n)
sex[sex==0]<-"M"
sex[sex==1]<-"F"
sex<-factor(sex)

# merge location and dive activity data
# 'negCons' is dummy covariate used to enforce turning angle constraints in the individual- and sex-level effects model
omega_sex<-data.frame(ID=factor(rep(1:N,times=T_n),levels=1:N),time=unlist(predTimes),omega=omega_t,sex=sex,negCons=rep(-1,length(sex)))
crwOut<-crawlMerge(crwOut,omega_sex,Time.name="time")

# prepare HMM data
hsData<-prepData(crwOut,covNames=c("sex","negCons"))

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
                 0,0,0,0,0,0,1,0,0,
                 0,0,0,0,0,0,0,1,0,
                 0,0,0,0,0,0,0,0,1),nrow=3*nbStates,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates),paste0("zeromass_",1:nbStates)),c(paste0("shape_",1:nbStates,":(Intercept)"),"scale:(Intercept)","scale_2","scale_3",paste0("zeromass_",1:nbStates,":(Intercept)"))))
stepcons<-c(1,1,1,1,2,2,1,1,1)
stepBounds<-matrix(c(0,5,
                     0,5,
                     0,5,
                     0,14400,
                     0,14400,
                     0,14400,
                     0,1,
                     0,1,
                     0,1),nrow=3*nbStates,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates),paste0("zeromass_",1:nbStates))))

angleDM<-matrix(c(1,0,0,
                  0,1,-1,
                  0,1,0),nrow=nbStates,byrow=TRUE,dimnames=list(paste0("concentration_",1:nbStates),c("concentration_1:(Intercept)","concentration_23:(Intercept)","concentration_2")))
angleBounds<-matrix(c(0,0.95,
                      0,0.95,
                      0,0.95),nrow=nbStates,byrow=TRUE,dimnames=list(paste0("concentration_",1:nbStates)))
anglecons<-c(1,2,2)
angleworkcons<-c(0,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),0)

omegaDM<-matrix(c(1,0,0,0,0,0,0,
                  0,0,1,1,0,0,0,
                  0,0,1,1,0,0,0,
                  1,1,0,0,0,0,0,
                  0,0,1,0,0,0,0,
                  0,0,1,0,0,0,0,
                  0,0,0,0,1,0,0,
                  0,0,0,0,0,1,0,
                  0,0,0,0,0,0,1),nrow=nbStates*3,byrow=TRUE,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates)),c("shape_1:(Intercept)","shape2_1","shape_2:(Intercept)","shape1_2",paste0("zeromass_",1:nbStates,":(Intercept)"))))
omegacons<-c(1,2,1,2,1,1,1)
omegaBounds<-matrix(c(1,10,
                      1,10,
                      1,10,
                      1,10,
                      1,10,
                      1,10,
                      0,1,
                      0,1,
                      0,1),nrow=nbStates*3,byrow=TRUE,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates))))


stepPar0<-c(2,2,2,1000,3000,8000,0.001,0.001,1.e-100)
anglePar0<-c(0.01,0.2,0.76)
omegaPar0<-c(1.1,9.9,9.9,9.9,1.1,1.1,c(0.5,1.e-100,1.e-100))

Par0<-getParDM(data=hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par=list(step=stepPar0,angle=anglePar0,omega=omegaPar0),DM=list(step=stepDM,angle=angleDM,omega=omegaDM),cons=list(step=stepcons,angle=anglecons,omega=omegacons),workcons=list(angle=angleworkcons),userBounds = list(step=stepBounds,angle=angleBounds,omega=omegaBounds))

fixPar<-list(step=c(rep(NA,nbStates*2),NA,NA,boot::logit(1.e-100)),
             omega=c(rep(NA,4),NA,boot::logit(c(1.e-100,1.e-100))))

bestFit<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par0,DM=list(step=stepDM,angle=angleDM,omega=omegaDM),cons=list(step=stepcons,angle=anglecons,omega=omegacons),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(angle=angleworkcons),fixPar=fixPar,stateNames=stateNames)

### Model 2: sex-level effects
# Note that factor-level covariates must be individually specified (e.g., 'sexF', 'sexM') when using pseudo-design matrix

stepDM.sex<-matrix(c("sexF",0,0,"sexM",0,0,0,0,0,0,0,0,0,0,0,0,0,
                     0,"sexF",0,0,"sexM",0,0,0,0,0,0,0,0,0,0,0,0,
                     0,0,"sexF",0,0,"sexM",0,0,0,0,0,0,0,0,0,0,0,
                     0,0,0,0,0,0,"sexF",0,0,"sexM",0,0,0,0,0,0,0,
                     0,0,0,0,0,0,"sexF","sexF",0,"sexM","sexM",0,0,0,0,0,0,
                     0,0,0,0,0,0,"sexF","sexF","sexF","sexM","sexM","sexM",0,0,0,0,0,
                     0,0,0,0,0,0,0,0,0,0,0,0,"sexF",0,0,"sexM",0,
                     0,0,0,0,0,0,0,0,0,0,0,0,0,"sexF",0,0,"sexM",
                     0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),nrow=3*nbStates,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates),paste0("zeromass_",1:nbStates)),c("shape_1:sexF","shape_2:sexF","shape_3:sexF","shape_1:sexM","shape_2:sexM","shape_3:sexM","scale_1:sexF","scale_2:sexF","scale_3:sexF","scale_1:sexM","scale_2:sexM","scale_3:sexM","zeromass_1:sexF","zeromass_2:sexF","zeromass_3:(Intercept)","zeromass_1:sexM","zeromass_2:sexM")))
stepcons.sex<-c(1,1,1,1,1,1,1,2,2,1,2,2,1,1,1,1,1)

angleDM.sex<-matrix(c("sexF",0,0,"sexM",0,0,
                      0,"sexF","sexF:negCons",0,"sexM","sexM:negCons",
                      0,"sexF",0,0,"sexM",0),nrow=nbStates,byrow=TRUE,dimnames=list(paste0("concentration_",1:nbStates),c("concentration_1:sexF","concentration_23:sexF","concentration_2:sexF","concentration_1:sexM","concentration_23:sexM","concentration_2:sexM")))
anglecons.sex<-c(1,2,2,1,2,2)
angleworkcons.sex<-c(0,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),0,0,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),0)

omegaDM.sex<-matrix(c("sexF",0,"sexM",0,0,0,0,0,0,0,0,0,
                      0,0,0,0,"sexF","sexF","sexM","sexM",0,0,0,0,
                      0,0,0,0,"sexF","sexF","sexM","sexM",0,0,0,0,
                      "sexF","sexF","sexM","sexM",0,0,0,0,0,0,0,0,
                      0,0,0,0,"sexF",0,"sexM",0,0,0,0,0,
                      0,0,0,0,"sexF",0,"sexM",0,0,0,0,0,
                      0,0,0,0,0,0,0,0,"sexF",0,0,"sexM",
                      0,0,0,0,0,0,0,0,0,1,0,0,
                      0,0,0,0,0,0,0,0,0,0,1,0),nrow=nbStates*3,byrow=TRUE,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates)),c("shape_1:sexF","shape2_1:sexF","shape_1:sexM","shape2_1:sexM","shape_2:sexF","shape1_2:sexF","shape_2:sexM","shape1_2:sexM","zeromass_1:sexF",paste0("zeromass_",2:nbStates,":(Intercept)"),"zeromass_1:sexM")))
omegacons.sex<-c(1,2,1,2,1,2,1,2,1,1,1,1)

fixPar.sex<-list(step=c(rep(NA,nbStates*2*2),NA,NA,boot::logit(1.e-100),NA,NA),
                 omega=c(rep(NA,4*2),NA,boot::logit(c(1.e-100,1.e-100)),NA))

bfPar<-getPar0(bestFit)
Par0.sex<-getPar0(bestFit,formula=~sex,DM=list(step=stepDM.sex,angle=angleDM.sex,omega=omegaDM.sex))
Par0.sex$Par$step[c("shape_1:sexF","shape_2:sexF","shape_3:sexF","shape_1:sexM","shape_2:sexM","shape_3:sexM")]<-bfPar$Par$step[c("shape_1:(Intercept)","shape_2:(Intercept)","shape_3:(Intercept)")]
Par0.sex$Par$step[c("scale_1:sexF","scale_2:sexF","scale_3:sexF","scale_1:sexM","scale_2:sexM","scale_3:sexM")]<-bfPar$Par$step[c("scale:(Intercept)","scale_2","scale_3")]

Par0.sex$Par$angle[c("concentration_1:sexF", "concentration_23:sexF",  "concentration_2:sexF",  "concentration_1:sexM", "concentration_23:sexM",  "concentration_2:sexM" )]<-bfPar$Par$angle
Par0.sex$Par$omega[c("shape_1:sexF", "shape2_1:sexF","shape_2:sexF", "shape1_2:sexF","zeromass_1:sexF")]<-bfPar$Par$omega[c("shape_1:(Intercept)","shape2_1","shape_2:(Intercept)","shape1_2", "zeromass_1:(Intercept)" )]
Par0.sex$Par$omega[c("shape_1:sexM", "shape2_1:sexM","shape_2:sexM", "shape1_2:sexM","zeromass_1:sexM")]<-bfPar$Par$omega[c("shape_1:(Intercept)","shape2_1","shape_2:(Intercept)","shape1_2", "zeromass_1:(Intercept)" )]

bestFit.sex<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par0.sex$Par,beta0=Par0.sex$beta,delta0=Par0.sex$delta,formula=~sex,DM=list(step=stepDM.sex,angle=angleDM.sex,omega=omegaDM.sex),cons=list(step=stepcons.sex,angle=anglecons.sex,omega=omegacons.sex),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(angle=angleworkcons.sex),fixPar=fixPar.sex,stateNames=stateNames)

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

stepcons.ind<-c(rep(1,3*N),rep(1,N),rep(2,2*N),rep(1,length(stepzm)*(nbStates-1)+1))
stepworkcons.ind<-c(rep(0,N*nbStates*2),rep(-boot::logit(1.e-100),length(stepzm)*(nbStates-1)),0)

angleDM.ind<-matrix(0,nrow=nbStates,ncol=N*nbStates,dimnames=list(c(paste0("concentration_",1:nbStates)),c(paste0("concentration_1:(Intercept)ID",1:N),paste0("concentration_23:(Intercept)ID",1:N),paste0("concentration_2ID",1:N))))
angleDM.ind[1,1:N]<-paste0("ID",1:N)
angleDM.ind[2,N+1:N]<-angleDM.ind[3,N+1:N]<-paste0("ID",1:N)
angleDM.ind[2,2*N+1:N]<-paste0("ID",1:N,":negCons")

anglecons.ind<-rep(c(1,2,2),each=N)
angleworkcons.ind<-rep(c(0,boot::logit((0.75 - angleBounds[3,1])/(angleBounds[3,2] - angleBounds[3,1])),0),each=N)

omegazm<-unique(hsData$ID[which(hsData$omega==0)])
omegaDM.ind<-matrix(0,nrow=nbStates*3,ncol=N*2*(nbStates-1)+length(omegazm)+1,dimnames=list(c(paste0("shape1_",1:nbStates),paste0("shape2_",1:nbStates),paste0("zeromass_",1:nbStates)),c(paste0("shape_1:(Intercept)ID",1:N),paste0("shape2_1ID",1:N),paste0("shape_2:(Intercept)ID",1:N),paste0("shape1_2ID",1:N),paste0("zeromass_1:(Intercept)ID",omegazm),paste0("zeromass_23:(Intercept)"))))
omegaDM.ind[1,1:N]<-omegaDM.ind[4,1:N]<-paste0("ID",1:N)
omegaDM.ind[4,N+1:N]<-paste0("ID",1:N)
omegaDM.ind[2,2*N+1:N]<-omegaDM.ind[3,2*N+1:N]<-omegaDM.ind[5,2*N+1:N]<-omegaDM.ind[6,2*N+1:N]<-paste0("ID",1:N)
omegaDM.ind[2,3*N+1:N]<-omegaDM.ind[3,3*N+1:N]<-paste0("ID",1:N)
omegaDM.ind[7,4*N+1:length(omegazm)]<-paste0("ID",omegazm)
omegaDM.ind[7:9,ncol(omegaDM.ind)]<-1
omegacons.ind<-c(rep(c(1,2,1,2),each=N),rep(1,length(omegazm)),1)
omegaworkcons.ind<-c(rep(0,N*2*(nbStates-1)),rep(-boot::logit(1.e-100),length(omegazm)),0)

fixPar.ind<-list(step=c(rep(NA,N*nbStates*2),rep(NA,length(stepzm)*(nbStates-1)),boot::logit(1.e-100)),
                 omega=c(rep(NA,N*2*(nbStates-1)),rep(NA,length(omegazm)),boot::logit(1.e-100)))

Par0.ind<-getPar0(bestFit,DM=list(step=stepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),formula=~ID+0)

# fit each track individually to get initial values for full model
registerDoParallel(cores=ncores) 
bestFit.all<-foreach(i=1:N) %dopar% {
  data.ind<-subset(hsData,ID==i)
  if(!all(data.ind$step>0,na.rm=TRUE)){
    fit<-fitHMM(data.ind,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=Par0,DM=list(step=stepDM,angle=angleDM,omega=omegaDM),cons=list(step=stepcons,angle=anglecons,omega=omegacons),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(angle=angleworkcons),fixPar=fixPar,stateNames=stateNames)
  }
  else {
    tmpstepDM<-stepDM[1:(2*nbStates),1:(2*nbStates)]
    tmpstepcons<-stepcons[1:(2*nbStates)]
    tmpstepBounds<-stepBounds[1:(2*nbStates),]
    tmpPar0<-Par0
    tmpPar0$step<-tmpPar0$step[1:(2*nbStates)]
    tmpfixPar<-fixPar
    tmpfixPar$step<-NULL
    fit<-fitHMM(data.ind,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar0,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),cons=list(step=tmpstepcons,angle=anglecons,omega=omegacons),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(angle=angleworkcons),fixPar=tmpfixPar,stateNames=stateNames)
  }
  fit
}
stopImplicitCluster()

#specify initial values for full model based on individual model fits
for(i in 1:N){
  bfPar<-getPar0(bestFit.all[[i]])
  Par0.ind$Par$step[match(paste0(colnames(stepDM),"ID",i),names(Par0.ind$Par$step),nomatch=0)]<-bfPar$Par$step[match(names(Par0.ind$Par$step),paste0(colnames(stepDM),"ID",i),nomatch=0)]
  Par0.ind$Par$angle[match(paste0(colnames(angleDM),"ID",i),names(Par0.ind$Par$angle),nomatch=0)]<-bfPar$Par$angle[match(names(Par0.ind$Par$angle),paste0(colnames(angleDM),"ID",i),nomatch=0)]
  Par0.ind$Par$omega[match(paste0(colnames(omegaDM),"ID",i),names(Par0.ind$Par$omega),nomatch=0)]<-bfPar$Par$omega[match(names(Par0.ind$Par$omega),paste0(colnames(omegaDM),"ID",i),nomatch=0)]
  Par0.ind$beta[paste0("ID",i),]<-bfPar$beta
}  
bestFit.ind<-fitHMM(hsData,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),formula=~ID+0,Par0=Par0.ind$Par,beta0=Par0.ind$beta,delta0=Par0.ind$delta,DM=list(step=stepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),cons=list(step=stepcons.ind,angle=anglecons.ind,omega=omegacons.ind),userBounds=list(step=stepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(step=stepworkcons.ind,angle=angleworkcons.ind,omega=omegaworkcons.ind),fixPar=fixPar.ind,stateNames=stateNames)

#################################################################################################################################
## compare models using AIC
#################################################################################################################################
AIC(bestFit,bestFit.sex,bestFit.ind)

# plot individual-level model
plot(bestFit.ind,plotCI=TRUE,ask=FALSE)

#################################################################################################################################
## Fit individual-level effects model using MIfitHMM
#################################################################################################################################

# fit MIfitHMM to each track individually
registerDoParallel(cores=ncores) 
miBestFit.all<-foreach(i=1:N) %dopar% {
  data.ind<-NULL
  data.ind$crwFits<-crwOut$crwFits[i]
  data.ind$crwPredict<-crwOut$crwPredict[which(crwOut$crwPredict$ID==i),]
  class(data.ind)<-c("crwData",class(data.ind))
  
  tmpstepDM<-stepDM[1:(2*nbStates),1:(2*nbStates)]
  tmpstepcons<-stepcons[1:(2*nbStates)]
  tmpstepBounds<-stepBounds[1:(2*nbStates),]
  tmpPar0<-getPar0(bestFit.all[[i]])
  tmpPar0$Par$step<-tmpPar0$Par$step[1:(2*nbStates)]
  tmpfixPar<-fixPar
  tmpfixPar$step<-NULL
  
  fit<-MIfitHMM(data.ind,nSims=nSims,ncores=1,poolEstimates=FALSE,
                nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),Par0=tmpPar0$Par,beta0=tmpPar0$beta,delta0=tmpPar0$delta,DM=list(step=tmpstepDM,angle=angleDM,omega=omegaDM),cons=list(step=tmpstepcons,angle=anglecons,omega=omegacons),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(angle=angleworkcons),fixPar=tmpfixPar,stateNames=stateNames,
                covNames=c("sex","negCons"))
}
stopImplicitCluster()

miSum.all<-lapply(miBestFit.all,MIpool)

tmpPar0.ind<-getPar0(bestFit.ind,DM=list(step=stepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),formula=~ID+0)

#specify initial values for full model based on individual model fits
for(i in 1:N){
  bfPar<-getPar0(miSum.all[[i]])
  tmpPar0.ind$Par$step[match(paste0(colnames(stepDM),"ID",i),names(Par0.ind$Par$step),nomatch=0)]<-bfPar$Par$step[match(names(Par0.ind$Par$step),paste0(colnames(stepDM),"ID",i),nomatch=0)]
  tmpPar0.ind$Par$angle[match(paste0(colnames(angleDM),"ID",i),names(Par0.ind$Par$angle),nomatch=0)]<-bfPar$Par$angle[match(names(Par0.ind$Par$angle),paste0(colnames(angleDM),"ID",i),nomatch=0)]
  tmpPar0.ind$Par$omega[match(paste0(colnames(omegaDM),"ID",i),names(Par0.ind$Par$omega),nomatch=0)]<-bfPar$Par$omega[match(names(Par0.ind$Par$omega),paste0(colnames(omegaDM),"ID",i),nomatch=0)]
  tmpPar0.ind$beta[paste0("ID",i),]<-bfPar$beta
}  
tmpstepDM.ind<-stepDM.ind[1:(2*nbStates),1:(N*2*nbStates)]
tmpstepcons.ind<-stepcons.ind[1:(N*2*nbStates)]
tmpstepBounds<-stepBounds[1:(2*nbStates),]
tmpPar0.ind$Par$step<-tmpPar0.ind$Par$step[1:(N*2*nbStates)]
tmpfixPar.ind<-fixPar.ind
tmpfixPar.ind$step<-NULL

#Fit multiple imputations (warning -- this takes a looooooong time!)
miBestFit.ind<-MIfitHMM(crwOut,nSims=nSims,ncores=ncores,poolEstimates=FALSE,
                        nbStates=nbStates,dist=list(step=stepDist,angle=angleDist,omega=omegaDist),formula=~ID+0,Par0=tmpPar0.ind$Par,beta0=tmpPar0.ind$beta,delta0=tmpPar0.ind$delta,DM=list(step=tmpstepDM.ind,angle=angleDM.ind,omega=omegaDM.ind),cons=list(step=tmpstepcons.ind,angle=anglecons.ind,omega=omegacons.ind),userBounds=list(step=tmpstepBounds,angle=angleBounds,omega=omegaBounds),workcons=list(angle=angleworkcons.ind,omega=omegaworkcons.ind),fixPar=tmpfixPar.ind,stateNames=stateNames,
                        covNames=c("sex","negCons"))

miSum.ind<-MIpool(miBestFit.ind,ncores=ncores)

plot(miSum.ind,plotCI=TRUE,ask=FALSE)

save.image("harbourSealExample.RData")
