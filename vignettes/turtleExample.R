library(momentuHMM)
library(sp)

# load turtle data from github
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/turtleData.RData"))

nSims <- 100 # number of imputatons
retryFits <- 250 # number attempt to re-fit based on random perturbation
ncores <- 7 # number of CPU cores

inits<-list(a=c(coordinates(turtleData)[1,1],0,
                coordinates(turtleData)[1,2],0),
            P=diag(c(100^2,100^2,100^2,100^2)))

fixPar<-c(log(1000*c(0.290,0.452,0.534,NA,NA,NA)),log(1000*c(0.122,0.239,0.301,NA,NA,NA)),NA,NA)

err.model=list(x=~lc-1,y=~lc-1)

constr=list(lower=c(rep(log(1000*0.534),3),rep(log(1000*0.301),3),rep(-Inf,2)),#c(rep(-Inf,2), log(1500),rep(-Inf,2)),#
            upper=rep(Inf,8))

predTimes <- seq(as.POSIXlt("2012-11-20 02:00:00 UTC",tz="UTC"),as.POSIXlt("2012-12-19 04:00:00 UTC",tz="UTC"),"2 hours")

crwOut<-crawlWrap(turtleData,retryFits=retryFits,Time.name="time",
                  err.model=err.model,initial.state=inits,
                  theta=c(7.730711, 8.216563, 8.505832, 7.103412, 7.245771, 7.935648, 5.371427, -10.677923),fixPar = fixPar, constr=constr,
                  predTime=predTimes) 

crwOut$crwPredict$date<-as.Date(crwOut$crwPredict$time)

miTurtleData<-MIfitHMM(crwOut,ncores=ncores,nSims=nSims,fit=FALSE,
                       spatialCovs=list(w=speedBrick,d=dirBrick,r=dirBrick),angleCovs="d")

for(j in 1:nSims){
  miTurtleData$miData[[j]]$bearing<-c(atan2(miTurtleData$miData[[j]]$y[2:nrow(miTurtleData$miData[[j]])]-miTurtleData$miData[[j]]$y[1:(nrow(miTurtleData$miData[[j]])-1)],miTurtleData$miData[[j]]$x[2:nrow(miTurtleData$miData[[j]])]-miTurtleData$miData[[j]]$x[1:(nrow(miTurtleData$miData[[j]])-1)]),NA)
  miTurtleData$miData[[j]]$angle_osc<-cos(miTurtleData$miData[[j]]$bearing-miTurtleData$miData[[j]]$r)
  miTurtleData$miData[[j]]$angle_osc[is.na(miTurtleData$miData[[j]]$angle_osc)]<-0
}

nbStates<-2
dist<-list(step="gamma",angle="wrpcauchy")
estAngleMean<-list(angle=TRUE)
circularAngleMean<-list(angle=TRUE)

DM<-list(step=list(mean=~state2(w:angle_osc),sd=~1),
         angle=list(mean=~state2(d),concentration=~1))

Par0<-list(step=c(9.132130, 9.200282, 0, 8.579123, 8.640819),
           angle=c(0, -16.0417715, -0.4749668))

turtleFits<-MIfitHMM(miTurtleData$miData,ncores=ncores,
                     nbStates=nbStates,dist=dist,Par0=Par0,DM=DM,estAngleMean=estAngleMean,circularAngleMean=circularAngleMean,retryFits=retryFits)

plot(turtleFits,plotCI=TRUE,covs=data.frame(angle_osc=1),ask=FALSE)

save.image("turtleExample.RData")
