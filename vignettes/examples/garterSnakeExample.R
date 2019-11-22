# hierarchical HMM garter snake example from Leos-Barajas et al (https://doi.org/10.1007/s13253-017-0282-9)

library(momentuHMM)
library(data.tree)

# load garter snake data from Leos-Barajas et al
load(url("https://static-content.springer.com/esm/art%3A10.1007%2Fs13253-017-0282-9/MediaObjects/13253_2017_282_MOESM1_ESM.rdata"))

W <- dim(dataAr)[3] # number of individuals
M <- dim(dataAr)[2] # number of time series per individual

### add 2 extra rows for each time step where coarse scale behavior switches occur (when every time series segment starts in this case)
# level=1  indicates when coarse scale behavior switching can occur (i.e., coarse scale t.p.m. used when level=1)
# level=2i indicates start of each fine scale interval (i.e., fine scale initial distribution used when level=2i)
# level=2  otherwise (i.e., fine scale t.p.m used when level=2)
snakeData <- NULL
for(w in 1:W){
  coarseInd <- data.frame(ID=w,level=c("1","2i"),step=NA)
  for(m in 1:M){
    tmp <- rbind(coarseInd,data.frame(ID=w,level="2",step=sqrt(dataAr[,m,w])))
    snakeData <- rbind(snakeData,tmp)
  }
}

# prepare hierarchical data
snakeData <- prepData(snakeData,coordNames=NULL,hierLevels=c("1","2i","2"))
plot(snakeData,dataNames="step",ask=FALSE)

# number of mixtures
mixtures <- 1

### define hierarchical HMM: states 1-3 = coarse state 1; states 4-6 = coarse state 2; states 7-9 = coarse state 3
hierStates <- Node$new("garter snake HHMM states")
hierStates$AddChild("internalState1")
hierStates$internalState1$AddChild("mo1", state=1) # motionless
hierStates$internalState1$AddChild("ex1", state=2) # slow exploratory
hierStates$internalState1$AddChild("es1", state=3) # rapid escape
hierStates$AddChild("internalState2")
hierStates$internalState2$AddChild("mo2", state=4) # motionless
hierStates$internalState2$AddChild("ex2", state=5) # slow exploratory
hierStates$internalState2$AddChild("es2", state=6) # rapid escape
hierStates$AddChild("internalState3")
hierStates$internalState3$AddChild("mo3", state=7) # motionless
hierStates$internalState3$AddChild("ex3", state=8) # slow exploratory
hierStates$internalState3$AddChild("es3", state=9) # rapid escape
# equivalent
#hierStates <- as.Node(list(name="garter snake HHMM states",
#                           internalState1=list(mo1=list(state=1),ex1=list(state=2),es1=list(state=3)),
#                           internalState2=list(mo2=list(state=4),ex2=list(state=5),es2=list(state=6)),
#                           internalState3=list(mo3=list(state=7),ex3=list(state=8),es3=list(state=9))))

print(hierStates,"state")

nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))

# data stream distributions: level 1 = coarse level (no data streams); level 2 = fine level (step="gamma")
hierDist <- Node$new("garter snake HHMM dist")
hierDist$AddChild("level1")
hierDist$AddChild("level2")
hierDist$level2$AddChild("step", dist="gamma")

# equivalent
#hierDist <- as.Node(list(name="garter snake HHMM dist",
#                           level1=list(),
#                           level2=list(step=list(dist="gamma"))))

print(hierDist,"dist")

# defining start values for step data stream
mu0 <- c(0.121,0.678,1.375)
sd0 <- c(0.06,0.321,0.4875)
Par0 <- list(step=c(rep(mu0,hierStates$count),rep(sd0,hierStates$count)))

# constrain fine-scale data stream distributions to be same for coarse-scale states
DM <- list(step=matrix(cbind(kronecker(c(1,1,1,0,0,0),diag(3)),
                             kronecker(c(0,0,0,1,1,1),diag(3))),
                       nrow=nbStates*2,
                       ncol=6,
                       dimnames=list(c(paste0("mean_",1:nbStates),paste0("sd_",1:nbStates)),
                                       paste0(rep(c("mean","sd"),each=3),c("_147:(Intercept)","_258:(Intercept)","_369:(Intercept)")))))

# get initial parameter values for data stream probability distributions on the working scale
Par <- getParDM(snakeData,hierStates=hierStates,hierDist=hierDist,Par=Par0,DM=DM)

# define hierarchical t.p.m. formulas (same as default)
hierFormula <- Node$new("garter snake formula")
hierFormula$AddChild("level1", formula=~1)
hierFormula$AddChild("level2", formula=~1)
# equivalent
#hierFormula <- as.Node(list(name="garter snake formula",
#                           level1=list(formula=~1),
#                           level2=list(formula=~1)))

# define hierarchical initial distribution formulas (same as default)
hierFormulaDelta <- Node$new("garter snake formulaDelta")
hierFormulaDelta$AddChild("level1", formulaDelta=~1)
hierFormulaDelta$AddChild("level2", formulaDelta=~1)
# equivalent
#hierFormulaDelta <- as.Node(list(name="garter snake formulaDelta",
#                           level1=list(formula=~1),
#                           level2=list(formula=~1)))

# initial values ('beta') and constraints ('betaCons') for t.p.m. at each level of hierarchy
# when mixtures>1, this will only include mixture at top level (via the betaCons argument)
level1states <- hierStates$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2) # reference states for level1
level1covNames <- colnames(model.matrix(hierFormula$level1$formula,snakeData))                        # covariate names for level1
level1dimNames <- list(paste0(level1covNames,"_mix",rep(1:mixtures,each=length(level1covNames))),
                       c(sapply(level1states,function(x) paste(rep(x,each=hierStates$count-1),"->",level1states[-which(level1states==x)]))))
hierBeta <- Node$new("garter snake beta")
hierBeta$AddChild("level1")
hierBeta$AddChild("level2")
hierBeta$level2$AddChild("internalState1")
hierBeta$level2$AddChild("internalState2")
hierBeta$level2$AddChild("internalState3")

betaCons <- data.tree::Clone(hierBeta)
betaCons$name <- "betaCons"
hierDelta <- data.tree::Clone(hierBeta)
hierDelta$name <- "hierDelta"
deltaCons <- data.tree::Clone(hierBeta)
deltaCons$name <- "deltaCons"

hierBeta$level1$beta <- matrix(c(1.24, 0.44, 1.1, -0.87, -1.40, -1.11),
                                        nrow=length(level1covNames)*mixtures,
                                        ncol=hierStates$count*(hierStates$count-1),byrow=TRUE,
                                        dimnames=level1dimNames)
betaCons$level1$betaCons <- matrix(seq(1,hierStates$count*(hierStates$count-1)*mixtures),
                                        nrow=length(level1covNames)*mixtures,
                                        ncol=hierStates$count*(hierStates$count-1),
                                        dimnames=level1dimNames)

beta0_level2 <- list()
beta0_level2$internalState1 <- c(-2.99, -5.06, 3.93, 1.25, 34.84, 35.97)
beta0_level2$internalState2 <- c(-1.72, -2.78, 3.54, 2.83, 34.72, 36.21)
beta0_level2$internalState3 <- c(-5.10, -17.69, 5.76, -11.34, 33.39, 37.37)
covNames <- colnames(model.matrix(hierFormula$level2$formula,snakeData)) # covariate names for level2
for(jj in 1:hierStates$count){
  j <-  names(hierStates$children)[jj]
  ref <- hierStates[[j]]$Get(function(x) Aggregate(x,"state",min),filterFun=function(x) x$level==2) # reference states for internalState j
  states <- hierStates[[j]]$Get("state",filterFun = isLeaf)                                         # states for internalState j
  dimNames <- list(paste0(covNames,"_mix",rep(1:mixtures,each=length(covNames))),paste0(rep(states,each=hierStates[[j]]$count-1)," -> ",states[-which(states==ref)]))
  hierBeta$level2[[j]]$beta <- matrix(beta0_level2[[j]],
                                          nrow=length(covNames)*mixtures,
                                          ncol=hierStates[[j]]$count*(hierStates[[j]]$count-1),byrow=TRUE,
                                          dimnames=dimNames)
  betaCons$level2[[j]]$betaCons <- matrix(rep(seq(1,hierStates$count*(hierStates$count-1)*mixtures,mixtures),each=length(covNames)*mixtures),
                                          nrow=length(covNames)*mixtures,
                                          ncol=hierStates[[j]]$count*(hierStates[[j]]$count-1),
                                          dimnames=dimNames)
}




# initial values ('delta') and constraints ('deltaCons') for initial distribution at each level of hierarchy
# when mixtures>1, this will only include mixture at top level (via the deltaCons field)
level1DeltaCovNames <- colnames(model.matrix(hierFormulaDelta$level1$formulaDelta,snakeData)) # delta covariate names for level1
level1DeltaDimNames <- list(paste0(level1DeltaCovNames,"_mix",rep(1:mixtures,each=length(level1DeltaCovNames))),paste0("state ",level1states[-1]))
hierDelta$level1$delta <- matrix(rep(c(-2.5,-3.5),hierStates$count-1),
                                        nrow=length(level1DeltaCovNames)*mixtures,
                                        ncol=(hierStates$count-1),byrow=TRUE,
                                        dimnames=level1DeltaDimNames)
deltaCons$level1$deltaCons <- matrix(seq(1,(hierStates$count-1)*mixtures),
                                          nrow=length(level1DeltaCovNames)*mixtures,
                                          ncol=(hierStates$count-1),
                                          dimnames=level1DeltaDimNames)

delta0_level2 <- list()
delta0_level2$internalState1 <- c(-1.39, 0.16)
delta0_level2$internalState2 <- c(-1.27, 2.33) 
delta0_level2$internalState3 <- c(0.34, -0.26)
level2DeltaCovNames <- colnames(model.matrix(hierFormulaDelta$level2$formulaDelta,snakeData)) # delta covariate names for level2
for(jj in 1:hierStates$count){
  j <- paste0("internalState",jj)
  states <- hierStates[[j]]$Get("state",filterFun = isLeaf) # states for internalState j
  deltaDimNames <- list(paste0(level2DeltaCovNames,"_mix",rep(1:mixtures,each=length(level2DeltaCovNames))),names(states)[-1])
  hierDelta$level2[[j]]$delta <- matrix(delta0_level2[[j]],
                                          nrow=length(level2DeltaCovNames)*mixtures,
                                          ncol=(hierStates[[j]]$count-1),byrow=TRUE,
                                          dimnames=deltaDimNames)
  deltaCons$level2[[j]]$deltaCons <- matrix(rep(seq(1,(hierStates$count-1)*mixtures,mixtures),each=length(level2DeltaCovNames)*mixtures),
                                            nrow=length(level2DeltaCovNames)*mixtures,
                                            ncol=(hierStates[[j]]$count-1),
                                            dimnames=deltaDimNames)
}

# hierBeta must be a list when mixtures>1
if(mixtures>1) hierBeta <- list(beta=hierBeta, pi=rep(1/mixtures,mixtures))

# check hierarchical model specification and parameters
checkPar0(snakeData,hierStates=hierStates,hierDist=hierDist,Par0=Par,hierFormula=hierFormula,hierFormulaDelta=hierFormulaDelta,mixtures=mixtures,DM=DM,hierBeta=hierBeta,hierDelta=hierDelta,betaCons=betaCons,deltaCons=deltaCons)

# compute hierarchical state transition probabilities based on initial values
iTrProbs <- getTrProbs(snakeData,hierStates=hierStates,hierBeta=hierBeta,hierFormula=hierFormula,hierDist=hierDist,mixtures=mixtures)
iTrProbs$level1$gamma[,,1] # tpm at first time step for level1
lapply(iTrProbs$level2,function(x) x$gamma[,,1]) # tpm at first time step for level2

# fit hierarchical HMM
hhmm <- fitHMM(snakeData,hierStates=hierStates,hierDist=hierDist,hierFormula=hierFormula,hierFormulaDelta=hierFormulaDelta,mixtures=mixtures,Par0=Par,DM=DM,hierBeta=hierBeta,hierDelta=hierDelta,betaCons=betaCons,deltaCons=deltaCons,nlmPar=list(print.level=2))
hhmm
plot(hhmm,ask=FALSE)
plotPR(hhmm)
plotStates(hhmm, animals="18", ask = FALSE)

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

### fine scale ###################################
# stationary distribution
lapply(stats[[1]]$level2,function(x) x[1,])
######################################################

### Simulate from fitted model ##############################
obsPerLevel<-Node$new("simHierData")
obsPerLevel$AddChild("level1",obs=M)               # number of level 1 observations
obsPerLevel$AddChild("level2",obs=dim(dataAr)[1])  # number of level 2 observations that follow each level 1 observation

simHHMM <- simHierData(nbAnimals=W,model=hhmm, obsPerLevel = obsPerLevel, states = TRUE)
plot(simHHMM,dataNames=hierDist$Get("name",filterFun=isLeaf),ask=FALSE)
#############################################################

save.image("garterSnakeExample.RData")