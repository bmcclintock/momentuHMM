library(momentuHMM)

# load ses data from github
tracks <- read.csv(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/sesData.csv"))
tracks <- tracks[,c(1,2,3)]
colnames(tracks) <- c("ID","x","y")

# derive steps and angles, and distance/angle to center
center <- matrix(c(70,-49),nrow=1,dimnames=list("colony"))
data <- prepData(data=tracks, type="LL", centers=center)

# standardize distances
data$colony.dist <- scale(data$colony.dist)

# time spent since left colony
time <- NULL
for(id in unique(data$ID)) {
  nbSubObs <- length(which(data$ID==id))
  
  # approximately in months for interval = 9.6h
  time <- c(time, (1:nbSubObs)/75)
}

data$time <- time

# include time since departure on 2->4 (column 6) and distance to colony on 1->2 (column 1) in t.p.m. formula
#formula <- ~ colony.dist + time
formula <- ~ betaCol1(colony.dist) + betaCol6(time)

stateNames <- c("outbound","search","forage","inbound")

#############################################################################################
# Fit a covariate-free model 
#############################################################################################
# initial parameters
stepPar0 <- c(25,5,1,25,10,5,3,10)
anglePar0 <- c(15,5,2,15)

# constrain transition probabilities
fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100),nrow=1)

m1 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), Par0=list(step=stepPar0, angle=anglePar0),
             fixPar=list(beta=fixbeta), stateNames = stateNames)
#############################################################################################



#############################################################################################
# biased random walk model of Michelot et al 2017 
#############################################################################################
# constrain transition probabilities
#fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100,
#                    NA,   0,   0,   0,NA,NA,   0,NA,   0,   0,   0,   0,
#                    NA,   0,   0,   0,NA,NA,   0,NA,   0,   0,   0,   0),
#                  nrow=3,byrow=TRUE)
fixbeta <- matrix(c(NA,-100,-100,-100,NA,NA,-100,NA,-100,-100,-100,-100,
                    NA,   0,   0,   0, 0, 0,   0, 0,   0,   0,   0,   0,
                     0,   0,   0,   0, 0,NA,   0, 0,   0,   0,   0,   0),
                  nrow=3,byrow=TRUE)

angleFormula <- ~ state1(colony.angle) + state4(colony.angle)

# constrain states 1 and 4 to biased random walks (outbound=repulsion, inbound=attraction)
fixPar <- list(angle=c(-100,100,NA,NA,NA,NA),beta=fixbeta)

Par0 <- getPar0(model=m1, nbStates=4, DM=list(angle=list(mean=angleFormula, concentration=~1)), 
                estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=TRUE), formula=formula)

m2 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
             beta0=Par0$beta, fixPar=fixPar, formula=formula, DM=list(angle=list(mean=angleFormula, concentration=~1)), 
             estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=TRUE), stateNames = stateNames)
#############################################################################################


#############################################################################################
# biased correlated random walk model including covariates on step length and turning angle concentration parameters
#############################################################################################
distFormula <- ~ state1(colony.dist) + state4(colony.dist)
stepDM <- list(mean=distFormula, sd=distFormula)
angleDM <- list(mean=angleFormula, concentration=distFormula)

fixPar <- list(beta=fixbeta)

Par0 <- getPar0(model=m2, nbStates=4, DM = list(step=stepDM, angle=angleDM), estAngleMean=list(angle=TRUE), 
                circularAngleMean=list(angle=TRUE), formula=formula)
Par0$Par$angle[c("mean_1:(colony.angle)","mean_4:(colony.angle)")] <- 0

m3 <- fitHMM(data=data, nbStates=4, dist=list(step="gamma",angle="vm"), Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
             beta0=Par0$beta, fixPar=fixPar, formula=formula, DM = list(step=stepDM, angle=angleDM), 
             estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=TRUE), stateNames = stateNames)
#############################################################################################

AIC(m1,m2,m3)

states <- viterbi(m3)

save.image("sesExample.RData")

#########################
## Plot decoded tracks ##
#########################
library(maps) # for map plots
library(mapdata) # for map plots
library(marmap) # to plot bathymetry
library(sp) # for degAxis

pal <- c("#78c679","#F0E442","#E31A1C","#88419d")

bathdat <- getNOAA.bathy(-20, 160, -72, -45, keep=TRUE)

pdf(file="plot_sesResults1.pdf",width=12,height=6)
plot(data$x, data$y, axes=FALSE, xlab="longitude", ylab="latitude", col="white")
degAxis(1)
degAxis(2)
plot(bathdat, add=TRUE,
     image=TRUE,
     step = 10000,
     deepest=-5000,
     shallowest=2000,
     bpal=colorRampPalette(c("dodgerblue3","lightsteelblue2"))(100),
     col=grey(0,0.1),
     lwd=0.1)
map('worldHires', add=TRUE, fill=TRUE, col='white')

for(id in unique(data$ID)) {
    ind <- which(data$ID==id)
    segments(x0 = data$x[ind[-length(ind)]], y0 = data$y[ind[-length(ind)]], x1 = data$x[ind[-1]], y1 = data$y[ind[-1]], 
             col = pal[states[ind[-length(ind)]]], lwd=2)
}

legend("topleft", legend = stateNames,
       col = pal, lwd=2, bg="white")
dev.off()