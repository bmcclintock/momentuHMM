## pilot whale example from Isojunno et al (2017; https://doi.org/10.1002/ecs2.2044)
## this example uses discrete individual-level random effects on state transition probabilities

## load data 
load(url("https://raw.github.com/bmcclintock/momentuHMM/develop/vignettes/pilotWhaleData.RData"))

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
delta0 <- c(0.32, 0.32, 0.15, 0.21)

stateNames <- c("exploratory","foraging","crowded","directed")

## fit model with single mixture on TPM
fitmix1 <- fitHMM(pilotData,nbStates=4,dist=dist,Par0=Par0,beta0=beta0,stationary=TRUE,mixtures=1,stateNames=stateNames,nlmPar=list(print.level=2))

## fit model with 2 mixtures on TPM
Par0_mix2 <- getPar0(fitmix1,mixtures=2)
Par0_mix2$beta$beta[1,] <- c(-2.26, -3.93, -0.58, 0.03, -2.25, -0.26, -3.38, -4.79, -2.82, -1.06, -3.3, -3.43)
Par0_mix2$beta$beta[2,] <- c(-2.51, -3.32, -2.63, 0.03, -1.26, -0.12, -96.8, -3.62, -1.75, -1.76, -2.14, -1.38)
Par0_mix2$beta$pi <- c(0.73, 0.27)
fitmix2 <- fitHMM(pilotData,nbStates=4,dist=dist,Par0=Par0_mix2$Par,beta0=Par0_mix2$beta,stationary=TRUE,mixtures=2,stateNames=stateNames,nlmPar=list(print.level=2))

## calculate AIC and AIC weights
AIC(fitmix1,fitmix2)
AICweights(fitmix1,fitmix2)

## calculate mixture probabilities for each individual
mixtureProbs(fitmix2)

## calculate state transition probabilities for each mixture
trProbs <- getTrProbs(fitmix2)
# mixture 1
round(trProbs[[1]][,,1],2)
# mixture 2
round(trProbs[[2]][,,1],2)

save.image("pilotWhaleExample.RData")
