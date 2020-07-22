library(momentuHMM)

# load turtle data from github
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/turtleData.RData"))

nSims <- 100 # number of imputatons
ncores <- 7  # number of CPU cores

fixPar<-c(log(1000*c(0.290,0.452,0.534,NA,NA,NA)),log(1000*c(0.122,0.239,0.301,NA,NA,NA)),NA,NA)

err.model=list(x=~lc-1,y=~lc-1)

constr=list(lower=c(rep(log(1000*0.534),3),rep(log(1000*0.301),3),rep(-Inf,2)),#c(rep(-Inf,2), log(1500),rep(-Inf,2)),#
            upper=rep(Inf,8))

predTimes <- seq(as.POSIXlt("2012-11-20 02:00:00 UTC",tz="UTC"),as.POSIXlt("2012-12-19 04:00:00 UTC",tz="UTC"),"2 hours")

# fit crawl model and predict locations at 2 hour time steps
crwOut<-crawlWrap(turtleData,retryFits=5,Time.name="time",
                  err.model=err.model,
                  theta=c(7.747184,   8.223590,   8.514041,   7.100828,   7.248787,   7.941202,   5.365669, -10.691232),fixPar = fixPar, constr=constr,
                  predTime=predTimes) 

# add date field to crwPredict data frame to match z values of raster bricks (speedBrick and dirBrick); see ?raster::getZ
crwOut$crwPredict$date<-as.Date(crwOut$crwPredict$time)

###################################################################################
## Fit 2-state HMM using multiple imputation
###################################################################################
nbStates<-2
dist<-list(mu="rw_mvnorm2")  # bivariate normal random walk on the positions (mu.x, mu.y)

DM <- list(mu=matrix(c("mu.x_tm1",0,"crw(mu.x_tm1,lag=1)",                    0, 0, 0,0,0,  # mu.x_tm1 is x location at previous time step; crw(mu.x_tm1,lag=1) is mu.x_tm1-mu.x_tm2 at lag 1 (previous x velocity)
                       "mu.x_tm1",0,                    0,"crw(mu.x_tm1,lag=1)","u",0,0,0,
                       0,"mu.y_tm1","crw(mu.y_tm1,lag=1)",                    0, 0, 0,0,0,  # mu.y_tm1 is y location at previous time step; crw(mu.y_tm1,lag=1) is mu.y_tm1-mu.y_tm2 at lag 1 (previous y velocity)
                       0,"mu.y_tm1",                    0,"crw(mu.y_tm1,lag=1)","v",0,0,0,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,1,1,0,
                       0,0,0,0,0,0,0,1,
                       0,0,0,0,0,0,0,1,
                       0,0,0,0,0,1,0,0,
                       0,0,0,0,0,1,1,0),10,8,byrow=TRUE,dimnames=list(NULL,c("x:x_tm1","y:y_tm1","xy_1:vel","xy_2:vel","xy:uv","sigma_12:(Intercept)","sigma_2","sigma.xy_12:(Intercept)"))))

Par0=list(mu=c(1, 1, 0.5, 0.7, 0.3, 15, 0.5, 0))
beta0 <- matrix(c(-3.5,-4),1,2)

fixPar <- list(mu=c(1,1,NA,NA,NA,NA,NA,0))

## full model using DM formulas (don't need to include mu.x_tm1 or mu.y_tm1 intercept terms when using formulas)
#DM <- list(mu=list(mean.x=~crw(mu.x_tm1,lag=1)+state2(u),mean.y=~crw(mu.y_tm1,lag=1)+state2(v),sigma.x=~1,sigma.xy=~1,sigma.y=~1))
#Par0 <- list(mu=   c(1, 0.7, 1, 0.9, 0.3, 1, 0.7, 1, 0.7, 0.4, 14, 15, 0, 0, 14, 15))
#fixPar <- list(mu=c(1,  NA, 1,  NA,  NA, 1,  NA, 1,  NA,  NA, NA, NA, 0, 0, NA, NA))

bestData <- prepData(crwOut,spatialCovs=list(u=uBrick,v=vBrick),altCoordNames="mu")     
checkPar0(bestData,nbStates=nbStates,dist=dist,DM=DM,Par0=Par0,beta0=beta0,fixPar=fixPar)

set.seed(1,kind="L'Ecuyer-CMRG",normal.kind="Inversion")
miData <- MIfitHMM(crwOut,ncores=ncores,nSims=nSims,
                   fit=FALSE,mvnCoords = "mu",
                   spatialCovs=list(u=uBrick,v=vBrick),altCoordNames = "mu")


# do initial optimization using penalization to help find decent starting values for each imputation
set.seed(1,kind="L'Ecuyer-CMRG",normal.kind="Inversion")
turtleFits<-MIfitHMM(miData,ncores=ncores,
                     nbStates=nbStates,dist=dist,Par0=Par0,beta0=beta0,DM=DM,fixPar=fixPar,
                     retryFits=10,mvnCoords = "mu",
                     spatialCovs=list(u=uBrick,v=vBrick),altCoordNames = "mu",prior=function(par) sum(dnorm(par[c(3,4,5,7,length(par)-2:1)],0,c(rep(1,4),rep(2,2)),log=TRUE)))

fPar0 <- lapply(turtleFits$HMMfits,function(x) getPar0(x)$Par)
fbeta0 <- lapply(turtleFits$HMMfits,function(x) getPar0(x)$beta)
fdelta0 <- lapply(turtleFits$HMMfits,function(x) getPar0(x)$delta)

# redo optimization without penalization
set.seed(1,kind="L'Ecuyer-CMRG",normal.kind="Inversion")
turtleFits<-MIfitHMM(miData,ncores=ncores,
                     nbStates=nbStates,dist=dist,Par0=fPar0,beta0=fbeta0,delta0=fdelta0,DM=DM,fixPar=fixPar,
                     retryFits=50,mvnCoords = "mu",
                     spatialCovs=list(u=uBrick,v=vBrick),altCoordNames = "mu")

###################################################################################

plot(turtleFits,plotCI=TRUE,ask=FALSE)
plotPR(turtleFits,ncores=ncores)

# plot relative to ocean surface current speed on 2 December 2012
plotSpatialCov(turtleFits,spatialCov=speedBrick$X2012.12.02)

# simulate from the model, raster spatial extent is too small for simulation so many will fail
set.seed(1,kind="L'Ecuyer-CMRG",normal.kind="Inversion")
simTurtle <- simData(model=turtleFits,spatialCovs=list(u=uBrick,v=vBrick),covs=turtleFits$miSum$data["date"],initialPosition = c(turtleFits$miSum$data$mu.x[1],turtleFits$miSum$data$mu.y[1]),retrySims = 100,states=TRUE)
plotSpatialCov(simTurtle,spatialCov=speedBrick$X2012.12.02,states=simTurtle$states)
plot(simTurtle,dataNames=c("mu.x","mu.y"))

save.image("turtleExample_rw_mvnorm2.RData")
