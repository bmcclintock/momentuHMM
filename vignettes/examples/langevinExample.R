library(momentuHMM)
library(raster)
library(viridis)
library(ggplot2)
library(patchwork)

#####################
## Load tracks ##
#####################
tracks <- read.csv(url("https://raw.github.com/bmcclintock/momentuHMM/develop/vignettes/SSLpreddat.csv"))
tracks$time <- as.POSIXct(tracks$time,format="%Y-%m-%d %H:%M:%S",tz="UTC")
tracks$time <- as.numeric(tracks$time)
tracks$time <- (tracks$time-min(tracks$time))/3600

#####################
## Load covariates ##
#####################
download.file(url="https://raw.github.com/bmcclintock/momentuHMM/develop/vignettes/aleut_habitat.gri?raw=true",destfile = "aleut_habitat.gri")
download.file(url="https://raw.github.com/bmcclintock/momentuHMM/develop/vignettes/aleut_habitat.grd?raw=true",destfile = "aleut_habitat.grd")
hbfull <- brick("aleut_habitat.grd", values=TRUE)
covlist0 <- list(bathy = hbfull$bathy,
                 slope = hbfull$slope,
                 d2site = hbfull$d2site)

# Convert to km
for(i in 1:length(covlist0)) {
  extent(covlist0[[i]]) <- extent(c(xmin(covlist0[[i]]), xmax(covlist0[[i]]), 
                                    ymin(covlist0[[i]]), ymax(covlist0[[i]]))/1000)
  projection(covlist0[[i]]) <- gsub("units=m", "units=km", projection(covlist0[[i]]))
}

ncov <- length(covlist0)
# Resample covariates to the same grid
for(i in 2:ncov)
  covlist0[[i]] <- resample(covlist0[[i]],covlist0[[1]])

# convert to km
tracks$x <- tracks$x/1000
tracks$y <- tracks$y/1000

# Crop covariates to area of interest
border <- 30
covlist0$bathy <- covlist0$bathy/1000
covlist0$d2site <- covlist0$d2site/1000
lim <- c(min(tracks$x)-border,max(tracks$x)+border,min(tracks$y)-border,max(tracks$y)+border)
covlist0 <- lapply(covlist0, crop, y=extent(lim))

# prepare data and calculate habitat gradients
langData <- prepData(tracks,
                     CT=TRUE, # identifies the data as being in continuous time
                     coordNames = c("x","y"),
                     altCoordNames = "mu",
                     spatialCovs = covlist0,
                     gradient=TRUE) # calculates gradients for spatial covariates

# specify design matrix
DM <- list(mu=matrix(c("mu.x_tm1","langevin(bathy.x)","langevin(slope.x)","langevin(d2site.x)",0,0,
                       "mu.y_tm1","langevin(bathy.y)","langevin(slope.y)","langevin(d2site.y)",0,0,
                       0,                  0,                  0,                   0,1,0,
                       0,                  0,                  0,                   0,1,0,
                       0,                  0,                  0,                   0,0,1),
                     nrow=5,byrow=TRUE,dimnames = list(c("mean.x","mean.y","sd.x","sd.y","corr.xy"),
                                                       c("mean:mu_tm1","bathy","slope","d2site","sd:(Intercept)","corr.xy:(Intercept)"))))

# fit model
fitLangevin <- fitCTHMM(langData,
                        nbStates=1,
                        dist=list(mu="rw_mvnorm2"),
                        DM=DM,
                        Par0=list(mu=c(1,0,0,0 ,1.25,0)),
                        fixPar=list(mu=c(1,rep(NA,3),NA,0)),
                        mvnCoords = "mu",
                        nlmPar=list(print.level=2))

# calculate utilization distribution
logUD <- fitLangevin$CIbeta$mu$est[2]*covlist0$bathy+fitLangevin$CIbeta$mu$est[3]*covlist0$slope + fitLangevin$CIbeta$mu$est[4]*covlist0$d2site
UD <- exp(logUD)/sum(exp(values(logUD)))

UDmat <- data.frame(coordinates(UD),val=values(UD))
ggtheme <- theme(axis.title = element_text(size=10), axis.text = element_text(size=10),
                 legend.title = element_text(size=12), legend.text = element_text(size=10),
                 title = element_text(size=10))

# Plot utilization distribution
p1 <- ggplot(UDmat,aes(x,y)) + geom_raster(aes(fill=val)) +
  coord_equal() + scale_fill_viridis(name=expression(pi)) +
  xlab("Easting (km)") + ylab("Northing (km)") + ggtheme +
  geom_point(aes(x,y), data=tracks, size=0.5)

# Plot log-utilization distribution
p2 <- ggplot(UDmat,aes(x,y)) + geom_raster(aes(fill=log(val))) +
  coord_equal() + scale_fill_viridis(name=expression(log(pi))) +
  xlab("Easting (km)") + ylab("Northing (km)") + ggtheme +
  geom_point(aes(x,y), data=tracks, size=0.5)

p1 + p2 + plot_layout(1,2)

# simulate from the fitted model
set.seed(18,"Mersenne-Twister","Inversion")
simLangevin <- simCTHMM(model=fitLangevin,
                        spatialCovs=covlist0,
                        matchModelObs = TRUE,
                        retrySims = 5)
p3 <- plotSpatialCov(simLangevin,log(UD),return=TRUE)
p3 + scale_fill_viridis(name=expression(log(pi)))
