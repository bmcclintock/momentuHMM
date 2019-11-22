# northern fulmar tracks modelled with biased movement relative to colony and fishing vessels
# this example uses data from Pirotta et al (2018; https://doi.org/10.1111/1365-2656.12830) and fits a very similar (but not identical) model
# 6 states: ARS at sea (1="seaARS"), transit at sea (2="seaTr"), ARS toward boat (3="boatARS"), transit toward boat (4="boatTR"), ARS toward colony (5="colonyARS"), and transit toward colony (6="colonyTr")
library(momentuHMM)
library(sp)

nbStates <- 6
stateNames <- c("seaARS","seaTr","boatARS","boatTr","colonyARS","colonyTr")

dist <- list(step="weibull",angle="wrpcauchy",d="lnorm") # distributions for step length ('step'), turning angle ('angle'), and distance to nearest boat ('d')

# load data provided by Pirotta et al
raw_data <- read.csv(url("https://datadryad.org/stash/downloads/file_stream/45899"),stringsAsFactors = FALSE)

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
                           0,Inf),ncol(stepDM),2,byrow=TRUE,
                         dimnames=list(colnames(stepDM),c("lower","upper")))


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
                            rep(c(0,Inf),nbTrips)),ncol(angleDM),2,byrow=TRUE,
                          dimnames=list(colnames(angleDM),c("lower","upper")))

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
                        0,Inf),ncol(dDM),2,byrow=TRUE,
                      dimnames=list(colnames(dDM),c("lower","upper")))

DM<-list(step=stepDM,angle=angleDM,d=dDM)
workBounds <- list(step=stepworkBounds,angle=angleworkBounds,d=dworkBounds)

# state transition formula similar to Pirotta et al
formula <- ~ toState3(boat.dist) + toState4(boat.dist) + toState5(time) + toState6(time)

# use betaCons to constrain: 1) ARS and Tr intercepts; and 2) boat and time effects to be the same for each sea, boat, and colony state
betaCons <- matrix(1:(3*nbStates*(nbStates-1)),3,nbStates*(nbStates-1),
                   dimnames=list(c("(Intercept)","boat.dist","time"),
                                 paste(rep(1:nbStates,each=nbStates),"->",rep(1:nbStates,nbStates))[(1:(nbStates*nbStates))[-diag(matrix(1:(nbStates*nbStates),nbStates,nbStates))]]))
betaCons["(Intercept)",c("1 -> 4","1 -> 6","3 -> 2","3 -> 4","3 -> 6","5 -> 2","5 -> 4","5 -> 6")] <- betaCons["(Intercept)","1 -> 2"] # constrain ARS -> Tr intercept
betaCons["(Intercept)",c("1 -> 5","3 -> 1","3 -> 5","5 -> 1","5 -> 3")] <- betaCons["(Intercept)","1 -> 3"] # constrain ARS -> ARS intercept
betaCons["(Intercept)",c("2 -> 3","2 -> 5","4 -> 1","4 -> 3","4 -> 5","6 -> 1","6 -> 3","6 -> 5")] <- betaCons["(Intercept)","2 -> 1"] # constrain Tr -> ARS intercept
betaCons["(Intercept)",c("2 -> 6","4 -> 2","4 -> 6","6 -> 2","6 -> 4")] <- betaCons["(Intercept)","2 -> 4"] # constrain Tr -> Tr intercept
betaCons["boat.dist",c("1 -> 4","2 -> 3","2 -> 4")] <- betaCons["boat.dist","1 -> 3"] # constrain boat.dist 1 -> 3 = 1 -> 4 = 2 -> 3 = 2 -> 4
betaCons["boat.dist","4 -> 3"] <- betaCons["boat.dist","3 -> 4"] # constrain boat.dist 3 -> 4 = 4 -> 3
betaCons["boat.dist",c("5 -> 4","6 -> 3","6 -> 4")] <- betaCons["boat.dist","5 -> 3"] # constrain boat.dist 5 -> 3 = 5 -> 4 = 6 -> 3 = 6 -> 4
betaCons["time",c("1 -> 6","2 -> 5","2 -> 6")] <- betaCons["time","1 -> 5"] # constrain time 1 -> 5 = 1 -> 6 = 2 -> 5 = 2 -> 6
betaCons["time",c("3 -> 6", "4 -> 5", "4 -> 6")] <- betaCons["time","3 -> 5"] # constrain time 3 -> 5 = 3 -> 6 = 4 -> 5 = 4 -> 6
betaCons["time","6 -> 5"] <- betaCons["time","5 -> 6"] # constrain time 5 -> 6 = 6 -> 5

# specify knownStates
knownStates <- rep(NA,nrow(fulmar_data))
knownStates[aInd] <- 2 # Priotta et al assumed all animals start in state 2 ('seaTr')

# fix delta_2 = 1 because assuming initial state is known for each track
fixPar <- list(delta=c(100,rep(0,nbStates-2)))
fixPar$delta <- exp(c(0,fixPar$delta))/sum(exp(c(0,fixPar$delta)))
fixPar$angle <- c(rep(1.e+7,3),rep(NA,2+2*nbTrips)) # fix angle mean parameters to large value and therefore constrain model to BRW (instead of BCRW)

# to speed things up, use starting values obtained after exploring the likelihood surface using retryFits
stepPar0 <- c(-0.06,0.92,6.24,0.83)
anglePar0 <- c(0,0,0,-0.78,0.29,-0.22,0.8,0.08,1.29,-0.88,0.21,-11.21,-13.76,-10.27,-58.07,-16.23,-14.69,-9.49,2.52)
dPar0 <- c(8.03,0.79,-0.34,-35.15)
beta0 <- matrix(c(-4.5,0,0,-5.08,-3.63,0,-4.5,-3.63,0,-5.08,0,1.66,-4.5,0,1.66,-3.8,0,0,-3.8,-3.63,0,-4.11,-3.63,0,-3.8,0,1.66,-4.11,0,1.66,-5.08,0,0,-4.5,0,0,-4.5,-3.73,0,-5.08,0,1.09,-4.5,0,1.09,-3.8,0,0,-4.11,0,0,-3.8,-3.73,0,-3.8,0,1.09,-4.11,0,1.09,-5.08,0,0,-4.5,0,0,-5.08,-1.93,0,-4.5,-1.93,0,-4.5,0,1.68,-3.8,0,0,-4.11,0,0,-3.8,-1.93,0,-4.11,-1.93,0,-3.8,0,1.68),3,nbStates*(nbStates-1))
Par0 <- list(Par=list(step=stepPar0,angle=anglePar0,d=dPar0),beta0=beta0)

#############################################################################################
# fit 6 state model with trip effects on sea states concentration parameter
#############################################################################################
m2<-fitHMM(fulmar_data,nbStates,dist,Par0=Par0$Par,beta0=Par0$beta0,formula=formula,estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),
           DM=DM,
           workBounds=workBounds,betaCons=betaCons,fixPar=fixPar,knownStates=knownStates,stateNames=stateNames)
#############################################################################################

plot(m2,covs=data.frame(sea.angle=0,boat.angle=0,colony.angle=0),ask=FALSE)

#activity budgets
#Pirotta et al report 0.277 (seaARS), 0.203 (seaTr), 0.182 (boatARS), 0.06 (boatTr), 0.115 (colonyARS), and 0.163 (colonyTr)
timeInStates(m2)

#activity budgets by individual bird
timeInStates(m2, by="birdID")

plotSat(m2,zoom=7,shape=c(17,1,17,1,17,1),size=2,col=rep(c("#E69F00", "#56B4E9", "#009E73"),each=2),stateNames=c("sea ARS","sea Transit","boat ARS","boat Transit","colony ARS","colony Transit"),projargs=newProj,ask=FALSE)

# plot stationary state probabilities as a function of covariates
plotStationary(m2,covs=data.frame(boat.dist=-0.5),plotCI=TRUE,legend.pos="topright")
