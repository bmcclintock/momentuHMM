## pilot whale example from Isojunno et al (2017; https://doi.org/10.1002/ecs2.2044)
## this example uses discrete individual-level random effects on state transition probabilities

library(momentuHMM)

## load data 
load(url("https://raw.github.com/bmcclintock/momentuHMM/master/vignettes/pilotWhaleData.RData"))

## prepare data
names(data.keep)[1] <- "ID"
pilotData <- prepData(data.keep,coordNames = NULL)

## 11 data streams
dist <- list(dive.dur = "weibull",
             dive.depth = "gamma",
             GR.speed2 = "gamma",
             dive.pitchvar2 = "beta",
             breath.headchange = "vm",
             GR.size = "pois",
             GR.tight = "bern",
             dive.CS.pres = "bern",
             dive.SS.pres = "bern",
             presurf = "bern",
             postsurf = "bern")

## initial values
Par0 <- list(dive.dur = c(1.9, 2.72, 1.64, 4.21, 1.3, 8.14, 1.53, 0.79),
             dive.depth = c(10.23, 315.91, 10.74, 5.51, 4.93, 233.12, 6.32, 1.91),
             GR.speed2 = c(1.15, 1.32, 1.36, 1.57, 0.66, 0.51, 0.77, 0.76),
             dive.pitchvar2 = c(2.21, 2.88, 1.82, 3.18, 17.94, 6.04, 16.06, 55.36),
             breath.headchange = c(3.07, 5.65, 2.64, 18.02),
             GR.size = c(6, 7.39, 20.39, 9.52),
             GR.tight = c(0.89, 0.66, 0.76, 0.81),
             dive.CS.pres = c(0.76, 0.99, 0.41, 0.47),
             dive.SS.pres = c(0.72, 0.98, 0.39, 0.41),
             presurf = c(0.76, 0.99, 0.71, 0.71),
             postsurf = c(0.81, 0.96, 0.72, 0.66))

beta0 <- matrix(c(-2.38, -3.86, -1.22, 0.21, -1.6, -0.47, -3.56, -4.15, -2.29, -1.17, -3.05, -2.54),nrow=1)

stateNames <- c("exploratory","foraging","crowded","directed")

## fit model with single mixture on TPM
fitmix1 <- fitHMM(pilotData,nbStates=4,dist=dist,Par0=Par0,beta0=beta0,stationary=TRUE,mixtures=1,stateNames=stateNames,nlmPar=list(print.level=2))

## fit model with 2 mixtures on TPM
Par0_mix2 <- getPar0(fitmix1,mixtures=2)
Par0_mix2$beta$beta[1,] <- c(-2.26, -3.93, -0.58, 0.03, -2.25, -0.26, -3.38, -4.79, -2.82, -1.06, -3.3, -3.43)
Par0_mix2$beta$beta[2,] <- c(-2.51, -3.32, -2.63, 0.03, -1.26, -0.12, -96.8, -3.62, -1.75, -1.76, -2.14, -1.38)
Par0_mix2$beta$pi <- c(0.73, 0.27)

# check model specification
checkPar0(pilotData,nbStates=4,dist=dist,Par0=Par0_mix2$Par,beta0=Par0_mix2$beta,stationary=TRUE,mixtures=2,stateNames=stateNames)

# fit 2 mixture model
fitmix2 <- fitHMM(pilotData,nbStates=4,dist=dist,Par0=Par0_mix2$Par,beta0=Par0_mix2$beta,stationary=TRUE,mixtures=2,stateNames=stateNames,nlmPar=list(print.level=2))

## fit model with 3 mixtures on TPM
Par0_mix3 <- getPar0(fitmix2,mixtures=3)
Par0_mix3$beta$beta[1,] <- c(-2.15, -4.31, -1.09, 0.28, -1.88, -0.3, -3.5, -4.71, -3.11, -0.68, -2.49, -2.6)
Par0_mix3$beta$beta[2,] <- c(-2.5, -2.47, 0.63, -17.22, -13.18, 0.59, -3.92, -13.96, -2.27, -1.25, -3.57, -3.75)
Par0_mix3$beta$beta[3,] <- c(-2.71, -3.48, -3.01, -0.35, -1.12, -0.1, -96.8, -2.98, -1.53, -2.29, -2.07, -1.55)
Par0_mix3$beta$pi <- c(0.4, 0.4, 0.2)

# fit 3 mixture model
fitmix3 <- fitHMM(pilotData,nbStates=4,dist=dist,Par0=Par0_mix3$Par,beta0=Par0_mix3$beta,stationary=TRUE,mixtures=3,stateNames=stateNames,nlmPar=list(print.level=2))

## fit model with 4 mixtures on TPM
Par0_mix4 <- getPar0(fitmix3,mixtures=4)
Par0_mix4$beta$beta[1,] <- c(-2.28, -4.09, -1.21, 0.05, -1.74, -0.26, -3.32, -5.43, -2.99, -0.76, -2.31, -2.25)
Par0_mix4$beta$beta[2,] <- c(-2.47, -2.45, 0.68, -17.22, -13.18, 0.59, -3.97, -13.97, -2.27, -1.28, -3.57, -3.75)
Par0_mix4$beta$beta[3,] <- c(-2.73, -3.47, -3, -0.37, -1.14, -0.09, -96.8, -2.97, -1.53, -2.32, -2.07, -1.54)
Par0_mix4$beta$beta[4,] <- c(-1.65, -25.12, -0.36, 2.11, -13.24, -16.48, -2.9, -3.45, -3.56, -0.39, -3.71, -31.97)
Par0_mix4$beta$pi <- c(0.33, 0.4, 0.2, 0.07)

# fit 4 mixture model
fitmix4 <- fitHMM(pilotData,nbStates=4,dist=dist,Par0=Par0_mix4$Par,beta0=Par0_mix4$beta,stationary=TRUE,mixtures=4,stateNames=stateNames,nlmPar=list(print.level=2))

## fixed effects model
Par0_fix <- getPar0(fitmix4,formula=~0+ID,mixtures=1)
fitfix <- fitHMM(pilotData,nbStates=4,dist=dist,formula=~0+ID,stationary=TRUE,Par0=Par0_fix$Par,stateNames=stateNames,nlmPar=list(print.level=2))

## calculate AIC and AIC weights
AIC(fitmix1,fitmix2,fitmix3,fitmix4,fitfix)
AICweights(fitmix1,fitmix2,fitmix3,fitmix4,fitfix)

## calculate mixture probabilities for each individual
mixProbs <- mixtureProbs(fitmix3,getCI=TRUE)
round(mixProbs$est,3)

## calculate state transition probabilities for each mixture
trProbs <- getTrProbs(fitmix3, covIndex=1)
# mixture 1
round(trProbs[[1]][,,1],2)
# mixture 2
round(trProbs[[2]][,,1],2)
# mixture 3
round(trProbs[[3]][,,1],2)

save.image("pilotWhaleExample.RData")
