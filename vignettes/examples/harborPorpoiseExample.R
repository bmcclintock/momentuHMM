# hierarchical HMM harbor porpoise example from Leos-Barajas et al (https://doi.org/10.1007/s13253-017-0282-9)
# this implementation is virtually identical except that the coarse- and fine-scale initial distributions are not assumed to be the corresponding stationary distributions

library(momentuHMM)
library(data.tree)

# load harbor porpoise data from Leos-Barajas et al
load(url("https://static-content.springer.com/esm/art%3A10.1007%2Fs13253-017-0282-9/MediaObjects/13253_2017_282_MOESM2_ESM.rdata"))

# convert date_time to POSIX
data <- lapply(data,function(x) {x$date_time <- as.POSIXct(x$date_time,tz="UTC"); x})

### add 2 extra rows for each time step where coarse scale behavior switches occur (hourly in this case)
# level=1  covariate indicates when coarse scale behavior switching can occur (i.e., coarse scale t.p.m. used when level=1)
# level=2i covariate indicates start of each fine scale interval (i.e., fine scale initial distribution used when level=2i)
# level=2  otherwise (i.e., fine scale t.p.m used when level=2)
porpoiseData <- NULL
for(i in 1:length(data)){
  coarseInd <- data.frame(date_time=as.POSIXct(format(data[[i]]$date_time[1],format="%Y-%m-%d %H:%M"),tz="UTC"),level=c("1","2i"),dive_duration=NA,maximum_depth=NA,dive_wiggliness=NA)
  tmp <- rbind(coarseInd,data.frame(data[[i]],level="2"))
  porpoiseData <- rbind(porpoiseData,tmp)
}

# prepare hierarchical data
porpoiseData <- prepData(porpoiseData,coordNames=NULL,hierLevels=c("1","2i","2"))

# data summary
summary(porpoiseData,dataNames=names(porpoiseData)[-1])

### define hierarchical HMM: states 1-3 = coarse state 1 (nonforaging); states 4-6 = coarse state 2 (foraging)
hierStates <- Node$new("harbor porpoise HHMM states")
hierStates$AddChild("nonforaging")
hierStates$nonforaging$AddChild("nf1", state=1)
hierStates$nonforaging$AddChild("nf2", state=2)
hierStates$nonforaging$AddChild("nf3", state=3)
hierStates$AddChild("foraging")
hierStates$foraging$AddChild("f1", state=4)
hierStates$foraging$AddChild("f2", state=5)
hierStates$foraging$AddChild("f3", state=6)
# equivalent
#hierStates <- as.Node(list(name="harbor porpoise HHMM states",
#                           nonforaging=list(nf1=list(state=1),nf2=list(state=2),nf3=list(state=3)),
#                           foraging=list(f1=list(state=4),f2=list(state=5),f3=list(state=6))))

print(hierStates,"state")

nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))

# data stream distributions: level 1 = coarse level (no data streams); level 2 = fine level (dive_duration="gamma", maximum_depth="gamma", dive_wiggliness="gamma")
hierDist <- Node$new("harbor porpoise HHMM dist")
hierDist$AddChild("level1")
hierDist$AddChild("level2")
hierDist$level2$AddChild("dive_duration", dist="gamma")
hierDist$level2$AddChild("maximum_depth", dist="gamma")
hierDist$level2$AddChild("dive_wiggliness", dist="gamma")
# equivalent
#hierDist <- as.Node(list(name="harbor porpoise HHMM dist",
#                           level1=list(),
#                           level2=list(dive_duration=list(dist="gamma"),maximum_depth=list(dist="gamma"),dive_wiggliness=list(dist="gamma"))))

print(hierDist,"dist")

# defining start values
dd.mu0 = rep(c(5,30,100),hierStates$count)
dd.sigma0 = rep(c(5,15,40),hierStates$count)
md.mu0 = rep(c(5,15,40),hierStates$count)
md.sigma0 = rep(c(2,5,20),hierStates$count)
dw.mu0 = rep(c(2,10,40),hierStates$count)
dw.sigma0 = rep(c(2,10,20),hierStates$count)
dw.pi0 = rep(c(0.2,0.01,0.01),hierStates$count)

Par0 <- list(dive_duration=c(dd.mu0,dd.sigma0),
             maximum_depth=c(md.mu0,md.sigma0),
             dive_wiggliness=c(dw.mu0,dw.sigma0,dw.pi0))

# constrain fine-scale data stream distributions to be same for coarse-scale states
dw_DM <- matrix(cbind(kronecker(c(1,1,0,0,0,0),diag(3)),
                      kronecker(c(0,0,1,1,0,0),diag(3)),
                      kronecker(c(0,0,0,0,1,1),diag(3))),
                nrow=nbStates*3,
                ncol=9,
                dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates),paste0("zeromass_",1:nbStates)),
                                paste0(rep(c("mean","sd","zeromass"),each=3),c("_14:(Intercept)","_25:(Intercept)","_36:(Intercept)"))))

DM <- list(dive_duration=dw_DM[1:(2*nbStates),1:6],
           maximum_depth=dw_DM[1:(2*nbStates),1:6],
           dive_wiggliness=dw_DM)

# get initial parameter values for data stream probability distributions on the working scale
Par <- getParDM(porpoiseData,hierStates=hierStates,hierDist=hierDist,Par=Par0,DM=DM)

# define hierarchical t.p.m. formula(s)
hierFormula <- Node$new("harbor porpoise HHMM formula")
hierFormula$AddChild("level1", formula=~1)
hierFormula$AddChild("level2", formula=~1)
# equivalent
#hierFormula <- as.Node(list(name="harbor porpoise HHMM formula",
#                           level1=list(formula=~1),
#                           level2=list(formula=~1)))

# define hierarchical initial distribution formula(s)
hierFormulaDelta <- Node$new("harbor porpoise HHMM formulaDelta")
hierFormulaDelta$AddChild("level1", formulaDelta=~1)
hierFormulaDelta$AddChild("level2", formulaDelta=~1)


## set initial values for beta and delta to nudge coarse state 1 to "nonforaging" and coarse state 2 to "foraging"

# initial values ('beta') for t.p.m. at each level of hierarchy
hierBeta <- Node$new("harbor porpoise beta")
hierBeta$AddChild("level1",beta=matrix(c(-1, -1),1))
hierBeta$AddChild("level2")
hierBeta$level2$AddChild("nonforaging",beta=matrix(c(0,-1,1,0,1,1),1))
hierBeta$level2$AddChild("foraging",beta=matrix(c(-1,1,1,2,0,3),1))
print(hierBeta,"beta")

# initial values ('delta') for initial distribution at each level of hierarchy
hierDelta <- Node$new("harbor porpoise delta")
hierDelta$AddChild("level1",delta=matrix(0,1))
hierDelta$AddChild("level2")
hierDelta$level2$AddChild("nonforaging",delta=matrix(c(30, 30),1))
hierDelta$level2$AddChild("foraging",delta=matrix(c(-6, 2),1))
print(hierDelta,"delta")

# check hierarchical model specification and parameters
checkPar0(porpoiseData,hierStates=hierStates,hierDist=hierDist,Par0=Par,hierFormula=hierFormula,hierFormulaDelta=hierFormulaDelta,DM=DM,hierBeta=hierBeta,hierDelta=hierDelta)

# compute hierarchical state transition probabilities based on initial values
iTrProbs <- getTrProbs(porpoiseData,hierStates=hierStates,hierBeta=hierBeta,hierFormula=hierFormula,hierDist=hierDist)
iTrProbs$level1$gamma[,,1] # tpm at first time step for level1
lapply(iTrProbs$level2,function(x) x$gamma[,,1]) # tpm at first time step for level2

# fit hierarchical HMM
hhmm <- fitHMM(porpoiseData,hierStates=hierStates,hierDist=hierDist,hierFormula=hierFormula,hierFormulaDelta=hierFormulaDelta,Par0=Par,hierBeta=hierBeta,hierDelta=hierDelta,DM=DM,nlmPar=list(print.level=2))
hhmm
plot(hhmm,ask=FALSE)
plotPR(hhmm)
plotStates(hhmm, ask = FALSE)

### most likely state sequence #######################
states <- viterbi(hhmm)
hStates <- viterbi(hhmm,hierarchical=TRUE)
# coarse scale sequence
coarseStates <- hStates$level1
# fine scale sequence
fineStates <- hStates$level2
######################################################

# stationary distributions
stats <- stationary(hhmm)

### coarse scale #####################################
# stationary distribution
stats[[1]]$level1[1,]
######################################################

### fine scale #######################################
# stationary distribution
lapply(stats[[1]]$level2,function(x) x[1,])
######################################################

### Simulate from fitted model ##############################
obsPerLevel<-Node$new("simHierData")
obsPerLevel$AddChild("level1",obs=100) # number of level 1 observations
obsPerLevel$AddChild("level2",obs=25)  # number of level 2 observations that follow each level 1 observation

simHHMM <- simHierData(model=hhmm, obsPerLevel = obsPerLevel, states = TRUE)
plot(simHHMM,dataNames=hierDist$Get("name",filterFun=isLeaf),ask=FALSE)
#############################################################

save.image("harborPorpoiseExample.RData")
