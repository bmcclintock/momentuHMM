library(momentuHMM)
library(data.tree)

# load the data from Adam et al (available from https://doi.org/10.1111/2041-210X.13241)
load(url("https://raw.github.com/bmcclintock/momentuHMM/develop/vignettes/Atlantic_cod_data_set.RData"))

# coarse-scale data
data <- data.frame(level="1",step=steps,angle=angles,vertical=NA,time=0)

### add extra rows for fine-scale data
# level=1  covariate indicates when coarse scale behavior switching can occur (i.e., coarse scale t.p.m. used when level=1)
# level=2i covariate indicates start of each fine scale interval (i.e., fine scale initial distribution used when level=2i)
# level=2  otherwise (i.e., fine scale t.p.m used when level=2)
codData <- NULL
timeSeq <- seq(from=0,to=23+5/6,length=144) # time of day covariate
for(i in 1:nrow(data)){
  fineInd <- data.frame(level="2",step=NA,angle=NA,vertical=verticals[[i]],time=timeSeq)
  tmp <- rbind(data[i,,drop=FALSE],data.frame(level="2i",step=NA,angle=NA,vertical=NA,time=0),fineInd)
  codData <- rbind(codData,tmp)
}

# prepare hierarchical data
codData <- prepData(codData,coordNames=NULL,covNames="time",hierLevels=c("1","2i","2"))

# data summary
summary(codData,dataNames=names(codData)[-1])

### define hierarchical HMM: states 1-3 = coarse state 1 (resident/foraging); states 4-6 = coarse state 2 (mobile/foraging); states 7-9 = coarse state 3 (travelling/migrating)
hierStates <- data.tree::Node$new("cod HHMM states")
hierStates$AddChild("resForage")   # resident/foraging
hierStates$resForage$AddChild("rF1", state=1)
hierStates$resForage$AddChild("rF2", state=2)
hierStates$resForage$AddChild("rF3", state=3)
hierStates$AddChild("mobForage")  # mobile/foraging
hierStates$mobForage$AddChild("mF1", state=4)
hierStates$mobForage$AddChild("mF2", state=5)
hierStates$mobForage$AddChild("mF3", state=6)
hierStates$AddChild("transit")    # travelling/migrating
hierStates$transit$AddChild("t1", state=7)
hierStates$transit$AddChild("t2", state=8)
hierStates$transit$AddChild("t3", state=9)

print(hierStates,"state")

nbStates <- length(hierStates$Get("state",filterFun=data.tree::isLeaf))

# data stream distributions: level 1 = coarse level (step="gamma", angle="vm"); level 2 = fine level (vertical="gamma")
hierDist <- data.tree::Node$new("cod HHMM dist")
hierDist$AddChild("level1")
hierDist$level1$AddChild("step", dist="gamma")
hierDist$level1$AddChild("angle", dist="vm")
hierDist$AddChild("level2")
hierDist$level2$AddChild("vertical", dist="gamma")
# equivalent
#hierDist <- data.tree::as.Node(list(name="cod HHMM dist",
#                           level1=list(step=list(dist="gamma"),angle=list(dist="vm")),
#                           level2=list(vertical=list(dist="gamma"))))

print(hierDist,"dist")

### defining start values based on those reported by Adam et al
hm.mu0 <- c(5.482, 6.786, 14.914)
hm.sigma0 <- c(4.27, 4.714, 11.242)
ha.mu0 <- c(0.011, -0.299, 0.044)
ha.kappa0 <- c(1.571, 1.426, 2.15)
vm.mu0 <- vm.sigma0 <- vm.pi0 <- list()
vm.mu0[[1]] <- c(0.116, 0.303, 0.691)
vm.mu0[[2]] <- c(0.109, 0.056, 0.351)
vm.mu0[[3]] <- c(0.125, 0.514, 1.987)
vm.sigma0[[1]] <- c(0.096, 0.261, 0.636)
vm.sigma0[[2]] <- c(0.043, 0.047, 0.342)
vm.sigma0[[3]] <- c(0.109, 0.462, 1.878)
vm.pi0[[1]] <- c(0.014, 0.003, 1.050e-06)
vm.pi0[[2]] <- c(1.791e-08, 0.035, 0.002)
vm.pi0[[3]] <- c(0.012, 2.933e-04, 3.462e-09)


Par0 <- list(step=c(rep(hm.mu0,each=3),rep(hm.sigma0,each=3)),
             angle=c(rep(ha.mu0,each=3),rep(ha.kappa0,each=3)),
             vertical=c(unlist(vm.mu0),unlist(vm.sigma0),unlist(vm.pi0)))

# constrain coarse-scale data stream parameters for corresponding fine-scale states
DM <- list(step=matrix(kronecker(diag(6),c(1,1,1)),
                       nrow=2*nbStates,
                       ncol=6,
                       dimnames=list(paste0(rep(c("mean_","sd_"),each=nbStates),1:nbStates),
                                     c(paste0(rep(c("mean_","sd_"),each=3),1:length(hierStates$children),":(Intercept)")))))
DM$angle <- DM$step
dimnames(DM$angle) <- list(paste0(rep(c("mean_","concentration_"),each=nbStates),1:nbStates),
                           c(paste0(rep(c("mean_","concentration_"),each=3),1:length(hierStates$children),":(Intercept)")))

# get initial parameter values for data stream probability distributions on the working scale
Par <- getParDM(codData,hierStates=hierStates,hierDist=hierDist,Par=Par0,DM=DM,estAngleMean = list(angle=TRUE))

# define hierarchical t.p.m. formula(s)
hierFormula <- data.tree::Node$new("cod HHMM formula")
hierFormula$AddChild("level1", formula=~1)
hierFormula$AddChild("level2", formula=~cosinor(time, period=24))

hierBeta <- data.tree::Node$new("cod beta")
hierBeta$AddChild("level1",beta=matrix(c(-18.585, -2.86, -2.551, -1.641, -2.169, -2.415),1,length(hierStates$children)*(length(hierStates$children)-1)))
hierBeta$AddChild("level2")
hierBeta$level2$AddChild("resForage",beta=matrix(c(-2.562, -3.403,  2.765, -1.607,  2.273,  4.842, 
                                                   -0.665,  -0.26, -0.681, -0.149, -2.728, -2.798, 
                                                   -0.027,   0.26,  0.191,  0.667,  0.123, -0.262),3,length(hierStates$resForage$children)*(length(hierStates$resForage$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("mobForage",beta=matrix(c(-2.156, -3.662,   3.01,  0.597, -0.313,  2.897, 
                                                    0.067,  -1.22, -0.799, -0.797,   0.15,  0.379, 
                                                   -0.112, -0.195, -0.269, -0.215,  1.539,  0.728),3,length(hierStates$mobForage$children)*(length(hierStates$mobForage$children)-1),byrow=TRUE))
hierBeta$level2$AddChild("transit",  beta=matrix(c( -2.53, -4.279,  2.507, -0.228, 10.803, 12.873, 
                                                    -0.04,  1.221, -0.301,  0.284, -0.106, -0.077, 
                                                    0.629, -0.226, -0.253, -0.303,  0.011,  0.036),3,length(hierStates$transit$children)*(length(hierStates$transit$children)-1),byrow=TRUE))


hierDelta <- data.tree::Node$new("cod delta")
hierDelta$AddChild("level1",delta=matrix(c(15.776, 4.78),1))
hierDelta$AddChild("level2")
hierDelta$level2$AddChild("resForage",delta=matrix(c(-0.643, -2.416),1))
hierDelta$level2$AddChild("mobForage",delta=matrix(c(1.181, 0.46),1))
hierDelta$level2$AddChild("transit",delta=matrix(c(-0.357, -0.624),1))


# check hierarchical model specification and parameters
checkPar0(codData,hierStates=hierStates,hierDist=hierDist,Par0=Par,hierFormula=hierFormula,DM=DM,hierBeta=hierBeta,hierDelta=hierDelta,estAngleMean = list(angle=TRUE))

# compute hierarchical state transition probabilities based on initial values
iTrProbs <- getTrProbs(codData,hierStates=hierStates,hierBeta=hierBeta,hierFormula=hierFormula,hierDist=hierDist)
iTrProbs$level1$gamma[,,1] # tpm at first time step for level1
lapply(iTrProbs$level2,function(x) x$gamma[,,1]) # tpm at first time step for level2

hhmm <- fitHMM(codData,hierStates=hierStates,hierDist=hierDist,hierFormula=hierFormula,Par0=Par,hierBeta=hierBeta,hierDelta=hierDelta,DM=DM,estAngleMean = list(angle=TRUE),
               nlmPar=list(print.level=2,iterlim=2500,steptol=1e-06,stepmax=150))
hhmm
plot(hhmm,plotCI=TRUE,ask=FALSE)
plotPR(hhmm)
plotStates(hhmm,ask=FALSE)

### most likely state sequence #######################
states <- viterbi(hhmm)
hStates <- viterbi(hhmm,hierarchical=TRUE)
# coarse scale sequence
coarseStates <- hStates$level1
# fine scale sequence
fineStates <- hStates$level2
######################################################

# stationary distributions
stats <- stationary(hhmm, covs=data.frame(time=timeSeq))

### coarse scale #####################################
# stationary distribution
stats[[1]]$level1[1,]
######################################################

### fine scale #######################################
# stationary distribution
stats[[1]]$level2

# plot stationary distribution as function of time of day
plotStationary(hhmm, plotCI=TRUE)
######################################################

save.image("codExample.RData")