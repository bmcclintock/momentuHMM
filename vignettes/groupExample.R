##This example is based on simulation Scenario A for the group dynamic model of Langrock et al (2014; Methods in Ecology and Evolution 5: 190-199)
library(momentuHMM)

oldRNG<-setRNG::setRNG()

setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=4)

#########################################################################
## Simulate group centroid path as BCRW relative to origin ##############
#########################################################################
dist <- list(step="gamma", angle="vm")
nbObs <- 250

Parc <- list(step=c(15,10), angle = c(0.15,log(1)))
DMc <- list(angle=list(mean=~center1.angle,concentration=~1))

centroidData <- simData(nbStates=1,dist=dist,Par=Parc,DM=DMc,circularAngleMean=list(angle=TRUE),centers=matrix(0,1,2),obsPerAnimal = nbObs)
plot(centroidData,ask=FALSE)
#########################################################################

#########################################################################
## Simulate individual paths with state1 as BRW relative to centroid ####
#########################################################################
nbAnimals <- 20
nbStates<-2
stateNames <- c("group","solitary")

Par <- list(step=c(30,50,15,25), angle = c(1.e+7,log(2.5),log(5)))
DM <- list(angle=list(mean=~state1(centroid.angle),concentration=~1))
  
beta <- matrix(c(-2.944439,-1.734601),1,nbStates)

# calculate stationary distribution
gamma <- diag(nbStates)
gamma[!gamma] <- exp(beta)
gamma <- t(gamma)
gamma <- gamma/apply(gamma,1,sum)
delta <- solve(diag(nbStates) - t(gamma) + 1, rep(1, nbStates))

# draw random initial locations for each individual
initialPositions <- vector("list")
for (i in 1:nbAnimals) {
  initialPositions[[i]] <- runif(2, -10, 10)
}

groupData <- simData(nbAnimals=nbAnimals,nbStates=nbStates,dist=dist,Par=Par,beta=beta,delta=delta,DM=DM,circularAngleMean=list(angle=TRUE),centroids=list(centroid=data.frame(x=centroidData$x,y=centroidData$y)),obsPerAnimal=nbObs,initialPosition=initialPositions,states=TRUE,stateNames=stateNames)
plot(groupData,compact=TRUE,ask=FALSE)
#########################################################################

#########################################################################
## Fit group dynamic model to simulated paths ###########################
Par0 <- list(step=c(30,50,15,25), angle = c(0,log(2.5),log(5)))
groupFit <- fitHMM(groupData,nbStates=nbStates,dist=dist,Par=Par0,DM=DM,stationary=TRUE,estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),stateNames=stateNames)
plot(groupFit,ask=FALSE)
#########################################################################

#########################################################################
## Simulate individual paths with state1 as BCRW relative to centroid ###
## using angleStrength and vmConsensus                                ###
#########################################################################
dist2 <- list(step="gamma", angle="vmConsensus")

# using these parameter values, minimum and maximum concentrations for state 1 are abs(1-min(centroid.dist)*0.01)*0.5 and (1+max(centroid.dist)*0.01)*0.5, respectively
Par2 <- list(step=c(30,50,15,25), angle = c(0.01,log(0.5),log(5)))
DM2 <- list(angle=list(mean=~state1(angleStrength(centroid.angle,strength=centroid.dist)),kappa=~1))

groupData2 <- simData(nbAnimals=nbAnimals,nbStates=nbStates,dist=dist2,Par=Par2,beta=beta,delta=delta,DM=DM2,circularAngleMean=list(angle=TRUE),centroids=list(centroid=data.frame(x=centroidData$x,y=centroidData$y)),obsPerAnimal=nbObs,initialPosition=initialPositions,states=TRUE,stateNames=stateNames)
plot(groupData2,compact=TRUE,ask=FALSE)
#########################################################################

#########################################################################
## Fit group dynamic model to simulated paths ###########################
groupFit2 <- fitHMM(groupData2,nbStates=nbStates,dist=dist2,Par=Par0,DM=DM2,stationary=TRUE,estAngleMean=list(angle=TRUE),circularAngleMean=list(angle=TRUE),stateNames=stateNames)
plot(groupFit2,plotCI=TRUE,covs=data.frame(centroid.angle=0),ask=FALSE)
#########################################################################

setRNG::setRNG(oldRNG)