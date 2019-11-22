
#' Example data simulation
#'
#' Generate data used in other functions' examples and unit tests.
#' 
# @importFrom gstat gstat vgm
# @importFrom raster raster
# @importFrom sp gridded spplot
# @importFrom stats predict

oldRNG<-setRNG::setRNG()

setRNG::setRNG(kind="Mersenne-Twister",normal.kind="Inversion",seed=10)

# simulate data
nbAnimals <- 2
nbStates <- 2
nbCovs <- 2
mu<-c(15,150)
sigma<-c(10,20)
angleMean <- c(0,0)
kappa <- c(0.7,1.5)
stepPar <- c(mu,sigma)
anglePar <- c(angleMean,kappa)
stepDist <- "gamma"
angleDist <- "vm"
zeroInflation <- FALSE
obsPerAnimal <- 100
beta <- matrix(c(-1.5,-1.5,0.25,-0.25,-0.1,0.1),3,2,byrow=TRUE)

simPar <- list(nbAnimals=nbAnimals,nbStates=nbStates,angleMean=angleMean,dist=list(step=stepDist,angle=angleDist),zeroInflation=list(step=zeroInflation,angle=FALSE))

data <- simData(nbAnimals=nbAnimals,nbStates=nbStates,dist=list(step=stepDist,angle=angleDist),
                Par=list(step=stepPar,angle=anglePar),formula=~cov1+cov2,beta=beta,nbCovs=nbCovs,zeroInflation=list(step=zeroInflation),
                obsPerAnimal=obsPerAnimal,states=TRUE)

# estimate model
mu0 <- c(20,70)
sigma0 <- c(10,30)
kappa0 <- c(1,1)
stepPar0 <- c(mu0,sigma0)
anglePar0 <- c(-pi/2,pi/2,kappa0)
formula <- ~cov1+cos(cov2)
nbCovs <- length(attr(terms(formula), "term.labels"))

beta0 <- matrix(c(rep(-1.5,nbStates*(nbStates-1)),rep(0,nbStates*(nbStates-1)*nbCovs)),
                nrow=nbCovs+1,byrow=TRUE)
delta0 <- rep(1,nbStates)/nbStates

par0 <- list(Par=list(step=stepPar0,angle=anglePar0),formula=formula,nbCovs=nbCovs,beta0=beta0,
             delta0=delta0)

m <- fitHMM(data=data,nbStates=nbStates,Par0=list(step=stepPar0,angle=anglePar0),beta0=beta0,
              delta0=delta0,formula=formula,dist=list(step=stepDist,angle=angleDist),estAngleMean=list(angle=TRUE))

example <- list(m=m,simPar=simPar,par0=par0)

err.model <- list(x= ~ ln.sd.x - 1, y =  ~ ln.sd.y - 1, rho =  ~ error.corr)

obsData<-simObsData(data,lambda=2,errorEllipse=list(M=c(0,50),m=c(0,50),r=0))

crwOut <- crawlWrap(obsData,theta=c(4,0),fixPar=c(1,1,NA,NA),err.model=err.model)

bestData<-prepData(crwOut,covNames=c("cov1","cov2"))
bestFit<-fitHMM(bestData,nbStates=nbStates,Par0=list(step=stepPar0,angle=anglePar0),beta0=beta0,
                  delta0=delta0,formula=formula,dist=list(step=stepDist,angle=angleDist),estAngleMean=list(angle=TRUE))

bPar<-getPar(bestFit)

miFits<-MIfitHMM(crwOut,nSims=4,nbStates=nbStates,Par0=bPar$Par,beta0=bPar$beta,
                  formula=formula,dist=list(step=stepDist,angle=angleDist),estAngleMean=list(angle=TRUE),
                  covNames=c("cov1","cov2"))

miExample <- list(obsData=obsData,bPar=bPar)

#set.seed(3)
#buffer<-100000        # grid half-width
#res<-1000             # grid cell resolution
#xy <- expand.grid(seq(-buffer,buffer,res), seq(-buffer,buffer,res))
#names(xy) <- c("x","y")
#g <- gstat::gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=gstat::vgm(psill=0.025,model="Exp",range=buffer), nmax=20)
#forest <- predict(g, newdata=xy, nsim=1)
#forest$sim1<-scale(forest$sim1,center=mean(forest$sim1))
#forest$sim1<-forest$sim1/max(forest$sim1)
#forest$sim1[which(forest$sim1<0)]<-0
#sp::gridded(forest) = ~x+y
#sp::spplot(forest[1])
#forest<-raster::raster(forest)
#names(forest)<-"forest"

#usethis::use_data(example,miExample,forest,compress="xz",overwrite=TRUE)
usethis::use_data(example,miExample,compress="xz",overwrite=TRUE)

setRNG::setRNG(oldRNG)

