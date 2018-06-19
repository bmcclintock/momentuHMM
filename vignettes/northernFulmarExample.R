# northern fulmar tracks modelled with biased movement relative to colony and fishing vessels
# this example uses data from Pirotta et al (2018; https://doi.org/10.1111/1365-2656.12830) and fits a very similar (but not identical) model
# 6 states: ARS at sea (1="seaARS"), transit at sea (2="seaTr"), ARS toward boat (3="boatARS"), transit toward boat (4="boatTR"), ARS toward colony (5="colonyARS"), and transit toward colony (6="colonyTr")
library(momentuHMM)
library(sp)

nbStates <- 6
stateNames <- c("seaARS","seaTr","boatARS","boatTr","colonyARS","colonyTr")

dist <- list(step="weibull",angle="wrpcauchy",d="lnorm") # distributions for step length ('step'), turning angle ('angle'), and distance to nearest boat ('d')

# load data provided by Pirotta et al
raw_data <- read.csv(url("https://datadryad.org/bitstream/handle/10255/dryad.174482/Fulmar_trackingData.csv?sequence=1"),stringsAsFactors = FALSE)

raw_data$ID <- raw_data$tripID
raw_data$Date <- as.POSIXct(raw_data$Date,tz="UTC",format="%d/%m/%Y %H:%M") # not sure these are actually UTC, but doesn't matter for this analysis

# project data
oldProj <- CRS("+proj=longlat +datum=WGS84")
newProj <- CRS("+init=epsg:27700")
coordinates(raw_data) <- c("Longitude","Latitude")
proj4string(raw_data) <- oldProj
raw_data <- as.data.frame(spTransform(raw_data,newProj))
coordinates(raw_data) <- c("Boat_Longitude","Boat_Latitude")
proj4string(raw_data) <- oldProj
raw_data <- as.data.frame(spTransform(raw_data,newProj))

# first use prepData to calculate colony distance covariate (which is needed to calculate at sea covariate below)
colony <- data.frame(x=-3.1,y=59.12) # colony coordinates as reported by Pirotta et al
coordinates(colony) <- c("x","y")
proj4string(colony) <- oldProj
colony <- as.matrix(as.data.frame(spTransform(colony,newProj)))
rownames(colony) <- "colony"
colony_dist <- prepData(raw_data,coordNames=c("Longitude","Latitude"),centers=colony)

# calculate at sea covariate
sea.angle <- NULL
for(id in unique(colony_dist$ID)) {
  idat <- subset(colony_dist,ID==id)
  nbSubObs <- length(which(colony_dist$ID==id))
  max_dist <- as.numeric(idat[which.max(idat$colony.dist),c("x","y")])
  max_angle <- momentuHMM:::distAngle(colony,colony,max_dist)[2]
  sea.angle <- c(sea.angle, rep(max_angle,nbSubObs))
}
raw_data$sea.angle<-sea.angle

# calculate time since left colony covariate ('time')
time <- aInd <- NULL
for(id in unique(raw_data$ID)) {
  idInd <- which(raw_data$ID==id)
  aInd <- c(aInd,idInd[1])
  nbSubObs <- length(idInd)
  time <- c(time, (1:nbSubObs)/nbSubObs)
}
raw_data$time <- time

# get boat data into centroids argument format
boat_data <- list(boat=data.frame(Date=raw_data$Date,x=raw_data$Boat_Longitude,y=raw_data$Boat_Latitude))

# format and merge all data and covariates for analysis
fulmar_data <- prepData(raw_data,coordNames=c("Longitude","Latitude"),centers=colony,centroids=boat_data,covNames="time",angleCovs = "sea.angle")
fulmar_data$d <- fulmar_data$boat.dist # momentuHMM doesn't like data streams and covariates to have same name, so create identical data column with different name

# standarize boat.dist covariate
fulmar_data$boat.dist <- scale(fulmar_data$boat.dist)

# specify similar data stream probability distribution parameter constraints as Pirotta et al using pseudo-design matrices
stepDM <- matrix(c(1,0,0,0,
                   0,1,0,0,
                   1,0,0,0,
                   0,1,0,0,
                   1,0,0,0,
                   0,1,0,0,
                   0,0,1,0,
                   0,0,1,1,
                   0,0,1,0,
                   0,0,1,1,
                   0,0,1,0,
                   0,0,1,1),2*nbStates,4,byrow=TRUE,dimnames=list(c(paste0("shape_",1:nbStates),paste0("scale_",1:nbStates)),
                                                                  c("shape:ARS","shape:Tr","scale:(Intercept)","scale:Tr")))

# constrain scale parameters such that Tr > ARS 
stepworkBounds <- matrix(c(-Inf,Inf,
                           -Inf,Inf,
                           -Inf,Inf,
                           0,Inf),ncol(stepDM),2,byrow=TRUE)


nbTrips <- length(unique(fulmar_data$ID))
angleDM <- matrix(c("sea.angle",0,0,0,0,rep(0,2*nbTrips),
                    "sea.angle",0,0,0,0,rep(0,2*nbTrips),
                    0,"boat.angle",0,0,0,rep(0,2*nbTrips),
                    0,"boat.angle",0,0,0,rep(0,2*nbTrips),
                    0,0,"colony.angle",0,0,rep(0,2*nbTrips),
                    0,0,"colony.angle",0,0,rep(0,2*nbTrips),
                    0,0,0,1,0,paste0("ID",1:nbTrips),rep(0,nbTrips),         # allow seaARS angle concentration parameter to vary by trip
                    0,0,0,1,1,paste0("ID",1:nbTrips),paste0("ID",1:nbTrips), # allow seaTr angle concentration parameter to vary by trip
                    0,0,0,1,0,rep(0,2*nbTrips),
                    0,0,0,1,1,rep(0,2*nbTrips),
                    0,0,0,1,0,rep(0,2*nbTrips),
                    0,0,0,1,1,rep(0,2*nbTrips)),2*nbStates,3+2+2*nbTrips,byrow=TRUE,dimnames=list(c(paste0("mean_",1:nbStates),paste0("concentration_",1:nbStates)),
                                                                  c("mean:sea","mean:boat","mean:colony","concentration:(Intercept)","concentration:Tr",paste0("concentration:ID",1:nbTrips,":(Intercept)"),paste0("concentration:ID",1:nbTrips,":Tr"))))

# constrain concentration parameters such that Tr > ARS 
angleworkBounds <- matrix(c(-Inf,Inf,
                            -Inf,Inf,
                            -Inf,Inf,
                            -Inf,Inf,
                            0,Inf,
                            rep(c(-Inf,Inf),nbTrips),
                            rep(c(0,Inf),nbTrips)),ncol(angleDM),2,byrow=TRUE)

dDM <- matrix(c(1,1,0,0,
                1,1,0,0,
                1,0,0,0,
                1,0,0,0,
                1,1,0,0,
                1,1,0,0,
                0,0,1,1,
                0,0,1,1,
                0,0,1,0,
                0,0,1,0,
                0,0,1,1,
                0,0,1,1),2*nbStates,4,byrow=TRUE,dimnames=list(c(paste0("location_",1:nbStates),paste0("scale_",1:nbStates)),
                                                                  c("location:(Intercept)","location:noboat","scale:(Intercept)","scale:noboat")))

# constrain location and scale parameters such that sea and colony > boat
dworkBounds <- matrix(c(-Inf,Inf,
                        0,Inf,
                        -Inf,Inf,
                        0,Inf),ncol(dDM),2,byrow=TRUE)

DM<-list(step=stepDM,angle=angleDM,d=dDM)
workBounds <- list(step=stepworkBounds,angle=angleworkBounds,d=dworkBounds)

# state transition formula similar to Pirotta et al
formula <- ~ toState3(boat.dist) + toState4(boat.dist) + toState5(time) + toState6(time)

# specify knownStates
knownStates <- rep(NA,nrow(fulmar_data))
knownStates[aInd] <- 2 # Priotta et al assumed all animals start in state 2 ('seaTr')

# fix delta_2 = 1 because assuming initial state is known for each track
fixPar <- list(delta=c(100,rep(0,nbStates-2)))
fixPar$delta <- exp(c(0,fixPar$delta))/sum(exp(c(0,fixPar$delta)))
fixPar$angle <- c(rep(1.e+7,3),rep(NA,2+2*nbTrips)) # fix angle mean parameters to large value and therefore constrain model to BRW (instead of BCRW)

# to speed things up, use starting values obtained after exploring the likelihood surface using retryFits
stepPar0 <- c(-0.04,0.93,6.2,0.85)
anglePar0 <- c(0,0,0,-0.63,0.16,-0.3,0.65,-0.48,1.26,-0.79,-0.18,-12.43,-13.74,-1.31,-1.29,-15.96,-13.59,-9.71,2.62)
dPar0 <- c(8.01,0.8,-0.36,-5)
beta0 <- matrix(c(-2.59,0,0,-29.49,-54.49,0,-12.29,0.02,0,-8.05,0,5.47,-6.54,0,4.48,-2.19,0,0,-11.6,-0.54,0,-49.99,-95.78,0,-42.07,0,80.76,-20.89,0,39.1,-4.22,0,0,-6.53,0,0,-0.29,2.96,0,-5.35,0,0.53,-9.46,0,-4.31,-3.89,0,0,-2.56,0,0,-18.37,-27.52,0,-8.59,0,-3.75,-2.37,0,2.01,-2.41,0,0,-4.47,0,0,-11.47,0.4,0,-11.43,0.43,0,-4.93,0,3.46,-4.97,0,0,-10.44,0,0,-11.78,-0.15,0,-38.36,-72.02,0,0.48,0,-3.7),3,nbStates*(nbStates-1))
Par0 <- list(Par=list(step=stepPar0,angle=anglePar0,d=dPar0),beta0=beta0)

#############################################################################################
# fit 6 state model with trip effects on sea states concentration parameter
#############################################################################################
m2<-fitHMM(fulmar_data,nbStates,dist,Par0=Par0$Par,beta0=Par0$beta0,formula=formula,estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),
           DM=DM,
           workBounds=workBounds,fixPar=fixPar,knownStates=knownStates,stateNames=stateNames,
           prior=function(par) {sum(dnorm(par[27+1:90],0,100,log=TRUE))}) # use priors on t.p.m to help prevent numerical issues 
#############################################################################################

plot(m2,covs=data.frame(sea.angle=0,boat.angle=0,colony.angle=0),ask=FALSE)

#activity budgets
#Pirotta et al report 0.277 (seaARS), 0.203 (seaTr), 0.182 (boatARS), 0.06 (boatTr), 0.115 (colonyARS), and 0.163 (colonyTr)
timeInStates(m2)

#activity budgets by individual bird
timeInStates(m2, by="birdID")

plotSat(m2,zoom=7,shape=c(17,1,17,1,17,1),size=2,col=rep(c("#E69F00", "#56B4E9", "#009E73"),each=2),stateNames=c("sea ARS","sea Transit","boat ARS","boat Transit","colony ARS","colony Transit"),projargs=newProj,ask=FALSE)

# plot stationary state probabilities as a function of covariates
plotStationary(m2,covs=data.frame(boat.dist=-0.5),legend.pos="topright")
