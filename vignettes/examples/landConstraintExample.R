# load packages
library(momentuHMM)
library(ctmcmove)

# create raster based on forest 
boundary <- forest
boundary[boundary>0] <- NA
boundary <- raster::distance(boundary)
names(boundary) <- "boundary"
proj4string(boundary) <- sp::CRS("+init=epsg:4326") # boundary needs to have a CRS for ctmcmove::rast.grad
raster::plot(boundary)

# compute gradients in x- and y- directions
grad <- ctmcmove::rast.grad(boundary)
grad.x <- grad$rast.grad.x
grad.y <- grad$rast.grad.y 

# specify bivariate normal biased correlated random walk with potential function based on distance to "water"
# mu.x_tm1 is x location at previous time step; crw(mu.x_tm1,lag=1) is mu.x_tm1-mu.x_tm2 at lag 1 (previous x velocity)
# mu.y_tm1 is y location at previous time step; crw(mu.y_tm1,lag=1) is mu.y_tm1-mu.y_tm2 at lag 1 (previous y velocity)
dist <- list(mu="rw_mvnorm2") # bivariate normal random walk
DM <- list(mu=list(mean.x=~mu.x_tm1+crw(mu.x_tm1,lag=1)+grad.x, 
                   mean.y=~mu.y_tm1+crw(mu.y_tm1,lag=1)+grad.y, 
                   sigma.x=~1,
                   sigma.xy=~1,
                   sigma.y=~1))

# specify parameters; negative coefficients for gradients results in attraction to "water"
Par <- list(mu=c(1,0.75,-750000,1,0.75,-1500,log(100000),0,log(100000))) # link scale
names(Par$mu) <- c("mu.x_tm1","crw(mu.x_tm1,lag=1)","grad.x","mu.y_tm1","crw(mu.y_tm1,lag=1)","grad.y","sigma.x","sigma.xy","sigma.y")

# simulate and plot
set.seed(4,kind="Mersenne-Twister",normal.kind="Inversion")
simBound <- simData(nbStates=1, obsPerAnimal = 10000, dist=dist, Par=Par, DM=DM, spatialCovs=list(grad.x=grad.x,grad.y=grad.y), mvnCoords="mu",initialPosition=c(25000,75000))
plot(simBound,dataNames=c("mu.x","mu.y"),ask=FALSE)
raster::plot(boundary$boundary)
points(simBound$mu.x,simBound$mu.y,type="l")
