library(momentuHMM)
library(data.tree)

# load the data from Adam et al (available from https://doi.org/10.1111/2041-210X.13241)
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/horn_shark_data_set.RData"))

# coarse-scale data
data <- data.frame(ID=unlist(mapply(function(x) rep(paste0("seg",x),nrow(steps[[x]])-1),1:length(steps))),
                   level="1",
                   steps=unlist(lapply(steps,function(x) x$steps[-nrow(x)])),
                   stepCat=unlist(lapply(steps,function(x) x$cats[-nrow(x)])),
                   odbas=NA)
data$stepCat[which(is.na(data$steps))] <- NA

### add extra rows for fine-scale data
# level=1  covariate indicates when coarse scale behavior switching can occur (i.e., coarse scale t.p.m. used when level=1)
# level=2i covariate indicates start of each fine scale interval (i.e., fine scale initial distribution used when level=2i)
# level=2  otherwise (i.e., fine scale t.p.m used when level=2)
odbas <- unlist(odbas)
sharkData <- NULL
for(i in 1:nrow(data)){
  fineInd <- data.frame(ID=data$ID[i],level="2",steps=NA,stepCat=NA,odbas=odbas[(i-1)*50+1:50])
  tmp <- rbind(data[i,,drop=FALSE],data.frame(ID=data$ID[i],level="2i",steps=NA,stepCat=NA,odbas=NA),fineInd)
  sharkData <- rbind(sharkData,tmp)
}

# prepare hierarchical data
sharkData <- prepData(sharkData,coordNames=NULL,hierLevels=c("1","2i","2"))

# data summary
summary(sharkData,dataNames=names(sharkData)[-1])

### define hierarchical HMM: states 1-3 = coarse state 1 (high activity); states 4-6 = coarse state 2 (resting); states 7-9 = coarse state 3 (travelling)
hierStates <- data.tree::Node$new("shark HHMM states")
hierStates$AddChild("activity")   # zero distance travelled, high activity
hierStates$activity$AddChild("a1", state=1)
hierStates$activity$AddChild("a2", state=2)
hierStates$activity$AddChild("a3", state=3)
hierStates$AddChild("resting")  # zero distance travelled, low activity
hierStates$resting$AddChild("r1", state=4)
hierStates$resting$AddChild("r2", state=5)
hierStates$resting$AddChild("r3", state=6)
hierStates$AddChild("transit")    # travelling
hierStates$transit$AddChild("t1", state=7)
hierStates$transit$AddChild("t2", state=8)
hierStates$transit$AddChild("t3", state=9)

print(hierStates,"state")

nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))
nCat <- 8 # number of stepCat categories

# data stream distributions: level 1 = coarse level (stepCat="cat8"); level 2 = fine level (odbas="gamma")
hierDist <- data.tree::Node$new("shark HHMM dist")
hierDist$AddChild("level1")
hierDist$level1$AddChild("stepCat", dist=paste0("cat",nCat)) # categorical distribution with 8 categories (1, 2, 3, 4, 5, 6, 7, or 8)
hierDist$AddChild("level2")
hierDist$level2$AddChild("odbas", dist="gamma")
# equivalent
#hierDist <- data.tree::as.Node(list(name="shark HHMM dist",
#                           level1=list(stepCat=list(dist="cat8")),
#                           level2=list(odbas=list(dist="gamma"))))

print(hierDist,"dist")

### defining start values based on Adam et al
mu0 <- c(0.191, 0.323, 0.721, 0.084, 0.15, 0.228, 0.094, 0.191, 0.39) # real scale
sd0 <- c(0.047, 0.051, 0.248, 0.021, 0.025, 0.033, 0.026, 0.039, 0.159) # real scale
probs0 <- c(1e+10, -1e+10, 0.958, 5.064, 3.261, 4.021, 0.473, 2.568) # working scale

# constrain coarse-scale data stream parameters for corresponding fine-scale states
DM <- list(stepCat=matrix(cbind(c(rep(c(1,1,1),2),rep(0,(nCat-1)*nbStates-6)),kronecker(diag(nCat-1),c(0,0,0,0,0,0,1,1,1))),
                       nrow=(nCat-1)*nbStates,
                       ncol=8,
                       dimnames=list(paste0(rep(paste0("prob",1:(nCat-1),"_"),each=nbStates),1:nbStates),
                                     c(paste0(c("prob1_12",paste0("prob",1:(nCat-1),"_3")),":(Intercept)")))))

Par0 <- list(stepCat=probs0,
             odbas=c(mu0,sd0))
fixPar <- list(stepCat=c(1.e+10,-1.e+10,rep(NA,6)))

hierBeta <- data.tree::Node$new("shark beta")
hierBeta$AddChild("level1",beta=matrix(c(-0.651, -0.169, -1.884, -0.369, -1.596, -0.737),1,length(hierStates$children)*(length(hierStates$children)-1)))
hierBeta$AddChild("level2")
hierBeta$level2$AddChild("activity",beta=matrix(c(-2.902, -14.032, 3.059, -1.37, 8.098, 11.66),1,length(hierStates$activity$children)*(length(hierStates$activity$children)-1)))
hierBeta$level2$AddChild("resting",beta=matrix(c(-3.264, -14.279, 3.107, 0.252, 13.861, 16.468),1,length(hierStates$resting$children)*(length(hierStates$resting$children)-1)))
hierBeta$level2$AddChild("transit",beta=matrix(c(-3.21, -21.32, 3.463, -0.598, 14.636, 17.811),1,length(hierStates$transit$children)*(length(hierStates$transit$children)-1)))

hierDelta <- data.tree::Node$new("shark delta")
hierDelta$AddChild("level1",delta=matrix(c(0.582, 2.894),1,length(hierStates$children)-1))
hierDelta$AddChild("level2")
hierDelta$level2$AddChild("activity",delta=matrix(c(-0.001, -1.1),1,length(hierStates$activity$children)-1))
hierDelta$level2$AddChild("resting",delta=matrix(c(-0.103, -0.105),1,length(hierStates$resting$children)-1))
hierDelta$level2$AddChild("transit",delta=matrix(c(0.24, -0.777),1,length(hierStates$transit$children)-1))

# check hierarchical model specification and parameters
checkPar0(sharkData,hierStates=hierStates,hierDist=hierDist,Par0=Par0,DM=DM,hierBeta=hierBeta,hierDelta=hierDelta,fixPar=fixPar)

# compute hierarchical state transition probabilities based on initial values
iTrProbs <- getTrProbs(sharkData,hierStates=hierStates,hierBeta=hierBeta,hierDist=hierDist)
iTrProbs$level1$gamma[,,1] # tpm at first time step for level1
lapply(iTrProbs$level2,function(x) x$gamma[,,1]) # tpm at first time step for level2

hhmm <- fitHMM(sharkData,hierStates=hierStates,hierDist=hierDist,Par0=Par0,hierBeta=hierBeta,hierDelta=hierDelta,DM=DM,fixPar=fixPar,
               nlmPar=list(print.level=2,iterlim=1000,print.level=2,stepmax=25,hessian=FALSE))
hhmm
plot(hhmm,ask=FALSE)
plotPR(hhmm)
plotStates(hhmm,animals=1:5,ask=FALSE)

### most likely state sequence #######################
states <- viterbi(hhmm)
hStates <- viterbi(hhmm,hierarchical=TRUE)
# coarse scale sequence
coarseStates <- hStates$level1
# fine scale sequence
fineStates <- hStates$level2
######################################################

# transition probabilities
trProbs12 <- getTrProbs(hhmm,covIndex=c(1,3))

# stationary distributions
stats12 <- stationary(hhmm,covIndex=c(1,3))

### coarse scale #####################################
# t.p.m.
round(trProbs12$level1$gamma[,,1],3)

# stationary distribution
round(stats12[[1]]$level1[1,],3)
######################################################

### fine scale #######################################
# t.p.m.
lapply(trProbs12$level2,function(x) round(x$gamma[,,1],3))

# stationary distribution
lapply(stats12[[1]]$level2,function(x) round(x[1,],3))
######################################################

### Simulate from fitted model ##############################
obsPerLevel <- data.tree::Node$new("simHierData")
obsPerLevel$AddChild("level1",obs=7) # number of level 1 observations
obsPerLevel$AddChild("level2",obs=10)  # number of level 2 observations that follow each level 1 observation

simHHMM <- simHierData(model=hhmm, nbAnimals = length(unique(sharkData$ID)), obsPerLevel = obsPerLevel, states = TRUE)
plot(simHHMM,dataNames=hierDist$Get("name",filterFun=isLeaf),animals=1:5,ask=FALSE)
#############################################################

save.image("hornSharkExample.RData")
