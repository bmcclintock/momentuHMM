library(momentuHMM)
load("elephantData.RData")

inits<-list(a=c(elephantData$x[1],0,elephantData$y[1],0),P = diag(c(5000 ^ 2,10 * 3600 ^ 2, 5000 ^ 2, 10 * 3600 ^ 2)))

crwOut<-crawlWrap(elephantData,ncores=1,timeStep="hour",initial.state=inits,theta=c(2.761119, -8.195757),fixPar=c(NA,NA))
plot(crwOut,ask=FALSE)

### create momentuHMMData object from crwData object
dataHMM <- prepData(crwOut,covNames="temp")
### add cosinor covariate based on hour of day
dataHMM$hour <- as.integer(strftime(dataHMM$time, format = "%H", tz="GMT"))

stateNames<-c("encamped","exploratory")

DM <- list(step = list(mean = ~ temp * cosinor(hour, period = 24),
                       sd = ~ temp * cosinor(hour, period = 24)),
           angle = list(concentration = ~ temp))
formula <- ~ temp * cosinor(hour, period = 24)

elephantFit <- fitHMM(data = dataHMM, nbStates = 2, 
                      dist = list(step = "gamma", angle = "wrpcauchy"),
                      Par0 = Par0$Par, beta0 = Par0$beta,
                      DM = DM, formula = formula,
                      estAngleMean = list(angle=FALSE), stateNames = stateNames)

plot(elephantFit, plotCI = TRUE, covs = data.frame(hour=12))

save.image("elephantExample.RData")
