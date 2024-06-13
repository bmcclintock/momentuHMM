# load packages
library(momentuHMM)

# create raster based on forest 
boundary <- forest
boundary[boundary>0] <- NA
boundary <- raster::distance(boundary)
names(boundary) <- "boundary"
raster::plot(boundary)

# specify bivariate normal biased correlated random walk with potential function based on distance to "water"
# mu.x_tm1 is x location at previous time step; crw(mu.x_tm1,lag=1) is mu.x_tm1-mu.x_tm2 at lag 1 (previous x velocity)
# mu.y_tm1 is y location at previous time step; crw(mu.y_tm1,lag=1) is mu.y_tm1-mu.y_tm2 at lag 1 (previous y velocity)
dist <- list(mu="rw_mvnorm2") # bivariate normal random walk
DM <- list(mu=list(mean.x=~mu.x_tm1+crw(mu.x_tm1,lag=1)+boundary.x, 
                   mean.y=~mu.y_tm1+crw(mu.y_tm1,lag=1)+boundary.y, 
                   sd.x=~1,
                   sd.y=~1,
                   corr.xy=~1))

# specify parameters; negative coefficients for gradients results in attraction to "water"
Par <- list(mu=c(1,0.75,-1500,1,0.75,-1500,log(sqrt(100000)),log(sqrt(100000)),0)) # link scale
names(Par$mu) <- c("mu.x_tm1","crw(mu.x_tm1,lag=1)","boundary.x","mu.y_tm1","crw(mu.y_tm1,lag=1)","boundary.y","sd.x","sd.y","corr.xy")

# simulate and plot
set.seed(2,kind="Mersenne-Twister",normal.kind="Inversion")
simBound <- simData(nbStates=1, obsPerAnimal = 10000, dist=dist, Par=Par, DM=DM, spatialCovs=list(boundary=boundary), gradient=TRUE, mvnCoords="mu",initialPosition=c(25000,75000))
plot(simBound,dataNames=c("mu.x","mu.y"),ask=FALSE)
raster::plot(boundary)
points(simBound$mu.x,simBound$mu.y,type="l")
