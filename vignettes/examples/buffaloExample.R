library(momentuHMM)
library(raster)
library(ctmcmove)
library(sp)

## download buffalo data
load(url("https://github.com/henryrscharf/Hooten_et_al_EL_2018/raw/master/data/buffalo/buffalo_Cilla.RData"))

## download distance to water covariate raster
load(url("https://github.com/henryrscharf/Hooten_et_al_EL_2018/raw/master/data/buffalo/dist2sabie.RData"))

# Store the original data
original_values <- values(dist2sabie)
original_extent <- extent(dist2sabie)
original_res <- res(dist2sabie)

# Create new raster with proper CRS using the correct syntax
dist2sabie <- raster(original_extent, 
                     resolution = original_res,
                     crs = CRS("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))

# Set values and name
values(dist2sabie) <- original_values
names(dist2sabie) <- "dist2sabie"

## standardize dist2sabie based on slope of gradient
dist2sabie_scaled <- dist2sabie / mean(values(raster::terrain(dist2sabie, opt = "slope")), na.rm = T)

# calculate gradient
D_scaled <- ctmcmove::rast.grad(dist2sabie_scaled)

## W (recharge function covariates)
# near_sabie = indicator for <500m from water
intercept <- raster(dist2sabie)
values(intercept) <- 1
W <- stack(list("intercept" = intercept,
                "near_sabie" = dist2sabie < 0.5e3))
W_names <- names(W)

## orthogonalize W based on locations ----
W_ortho <- W
W_path <- extract(x = W, y = matrix(buffalo_proj@coords, ncol = 2))
obstimes <- as.numeric(buffalo_proj$POSIX) / 3600 # convert time to numeric hours
W_tilde <- apply(W_path * c(0, diff(obstimes)), 2, cumsum)
W_tilde_svd <- svd(W_tilde)
W_tilde_proj_mat <- W_tilde_svd$v %*% diag(W_tilde_svd$d^(-1))
W_mat <- as.matrix(W)
W_mat_proj <- W_mat %*% W_tilde_proj_mat
for(layer in 1:ncol(W_mat)){
  values(W_ortho[[layer]]) <- W_mat_proj[, layer]
  names(W_ortho[[layer]]) <- paste0("svd", layer)
}

lnError <- crawl::argosDiag2Cov(50,50,0) # assume 50m isotropic error ellipse
buffaloData <- data.frame(ID=1,
                          time=obstimes,
                          x=buffalo_proj@coords[, 1],
                          y=buffalo_proj@coords[, 2],
                          ln.sd.x=lnError$ln.sd.x, 
                          ln.sd.y = lnError$ln.sd.y, 
                          error.corr= lnError$error.corr)

# plot observed track and distance to water raster
plotSpatialCov(buffaloData,dist2sabie)

crwOut <- crawlWrap(buffaloData,theta=c(6.5,-.1),fixPar=c(1,1,NA,NA),
                    err.model = list(x=~ln.sd.x-1,y=~ln.sd.y-1,rho=~error.corr),
                    timeStep=0.25, # predict at 15 min time steps
                    attempts=10)

spatialCovs <- list(W_intercept=W_ortho$svd1,
                    W_near_sabie=W_ortho$svd2,
                    dist2sabie=dist2sabie,
                    D.x=D_scaled$rast.grad.x,
                    D.y=D_scaled$rast.grad.y)

hmmData <- prepData(crwOut, spatialCovs = spatialCovs, altCoordNames = "mu")

nbStates <- 2
stateNames <- c("charged","discharged")
dist <- list(mu="rw_mvnorm2") # bivariate normal random walk (mu.x, mu.y)

DM <- list(mu=matrix(c("mu.x_tm1",         0,    0,0,0,0,
                       "mu.x_tm1",         0,"D.x",0,0,0,
                                0,"mu.y_tm1",    0,0,0,0,
                                0,"mu.y_tm1","D.y",0,0,0,
                                0,         0,    0,1,0,0,
                                0,         0,    0,0,1,0,
                                0,         0,    0,0,0,1,
                                0,         0,    0,0,0,1,
                                0,         0,    0,1,0,0,
                                0,         0,    0,0,1,0),
                     5*nbStates,
                     6,byrow=TRUE,
                     dimnames=list(c(paste0("mean.",rep(c("x_","y_"),each=nbStates),1:nbStates),
                                     paste0("sigma.",rep(c("x_","xy_","y_"),each=nbStates),1:nbStates)),
                                   c("x:x_tm1",
                                     "y:y_tm1",
                                     "xy:D",
                                     "sigma_1:(Intercept)",
                                     "sigma_2:(Intercept)",
                                     "sigma_12:(Intercept)"))))

# starting values
Par0=list(mu=c(1, 1, 0, log(85872.66), log(37753.53), 0))
g0 <- 0 # recharge function at time 0
theta <- c(0,0,0) # recharge function parameters

## specify recharge formula
# note that formula for theta requires an 'intercept' term, so we fix it to zero using fixPar
formula <- ~ recharge(g0=~1,theta=~W_intercept+W_near_sabie)

## remove Markov property
betaRef <- c(1,1) # make state 1 the reference state
betaCons <- matrix(c(1,2),2,2) # 1 -> 1 = 2 -> 1 and 1 -> 2 = 2 -> 2

## set fixed parameters
# no intercept or estimated coefficients for t.p.m based entirely on recharge model
# recharge coefficient set to -1 because charged state (state 1) is the reference state
# initial distribution doesn't affect likelihood when Markov property is removed
fixPar <- list(mu = c(Par0$mu[1:2],NA,NA,NA,Par0$mu[6]),
               beta = matrix(c(0,-1,0,-1),2,2),
               delta = c(0.5,0.5),
               theta = c(0,NA,NA)) # fix extra 'intercept' term to zero

# check recharge model specification
checkPar0(hmmData,nbStates=nbStates,dist=dist,formula=formula,
          Par0=Par0,beta0=list(beta=fixPar$beta,g0=g0,theta=theta),delta0=fixPar$delta,
          fixPar=fixPar,DM=DM,betaRef=betaRef,betaCons=betaCons,stateNames=stateNames)

# fit to best predicted path from crawl to get starting values for multiple imputation
buffaloFit <- fitHMM(hmmData,nbStates=nbStates,dist=dist,formula=formula,
                     Par0=Par0,beta0=list(g0=g0,theta=theta),
                     fixPar=fixPar,DM=DM,betaRef=betaRef,betaCons=betaCons,stateNames=stateNames,
                     mvnCoords="mu",optMethod="Nelder-Mead",control=list(maxit=1000))

# multiple imputation fits
bestPar <- getPar(buffaloFit)
set.seed(1,kind="Mersenne-Twister",normal.kind="Inversion")
buffaloFits <- MIfitHMM(crwOut, nSims=28, ncores=4,
                        spatialCovs = spatialCovs, 
                        mvnCoords="mu", altCoordNames = "mu",
                        nbStates=nbStates, dist=dist, formula=formula,
                        Par0=bestPar$Par, beta0=bestPar$beta,
                        fixPar=fixPar, DM=DM, betaRef=betaRef, betaCons=betaCons, stateNames=stateNames,
                        retryFits=3, retrySD=list(mu=c(0,0,3,0,0,0),g0=1,theta=c(0,1,1)),
                        optMethod="Nelder-Mead",control=list(maxit=100000))
buffaloFits
plot(buffaloFits,plotCI=TRUE,ask=FALSE)
plotSpatialCov(buffaloFits,dist2sabie)

trProbs <- getTrProbs(buffaloFits, getCI=TRUE)
# plot estimates and CIs for Pr(discharged) at each time step
plot(trProbs$est[1,2,],type="l", 
     ylim=c(0,1), ylab="Pr(discharged)", xlab="t", col=c("#E69F00", "#56B4E9")[buffaloFits$miSum$Par$states])
arrows(1:dim(trProbs$est)[3],
       trProbs$lower[1,2,],
       1:dim(trProbs$est)[3],
       trProbs$upper[1,2,],
       length=0.025, angle=90, code=3, col=c("#E69F00", "#56B4E9")[buffaloFits$miSum$Par$states], lwd=1.3)
abline(h=0.5,lty=2)

# proportion of entire time series spent in each state
buffaloFits$miSum$Par$timeInStates

# histograms of distance to water by state
par(mfrow=c(2,1))
hist(buffaloFits$miSum$data$dist2sabie[which(buffaloFits$miSum$Par$states==1)],main=stateNames[1],xlab="distance to water (m)")
hist(buffaloFits$miSum$data$dist2sabie[which(buffaloFits$miSum$Par$states==2)],main=stateNames[2],xlab="distance to water (m)")

save.image("buffaloExample.RData")
