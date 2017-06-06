
library(moveHMM)					
set.seed(1122334455)					
## Simulate covariate values 					
# (slopes in degrees	 ranging between 0 and 40	 and temperatures	 fluctuating 		
# around 10 degrees - the latter won't affect the state switching in the model below).					
# The model will be such that the haggises are most likely to be in the exploratory 					
# state (state 2) when at slopes of around 20 degrees (the slope that perfectly matches 					
# the difference in their leg lengths). For slopes close to 0 and slopes close to 40	 				
# the animals become essentially immobile (due to the differences in their leg lengths) 					
# and hence need to slowly crawl back (state 1) to slope levels better suited for them.					
# We are using a quadratic predictor to achieve this setup.					
slopes <- NULL					
for (haggis in 1:15) {					
  # for each of the 15 haggises	 simulate 400 slope values				
  arsim <- rep(NA,	400)				
  arsim[1] <- runif(1,	0.5	,1)			
  for (k in 2:400)					
    arsim[k] <- 0.9*(arsim[k-1]-0.4)+0.4+rnorm(1,	0,	0.7)			
  
  slope <- 40*plogis(arsim)					
  slopes <- c(slopes,	slope)				
}					
temps <- NULL					
for (haggis in 1:15) {					
  # for each of the 15 haggises	 simulate 400 temperature values				
  arsim2 <- rep(NA,	400)				
  arsim2[1] <- rnorm(1,	10,	3)			
  for (k in 2:400)					
    arsim2[k] <- 0.9*(arsim2[k-1]-10)+10+rnorm(1,	0	,2)			
  
  temps <- c(temps,	arsim2)				
}					
# data frame of covariate values					
covs <- data.frame(slope=slopes,	slope2=slopes^2,	temp=temps)			
# specify parameters of the step length and turning angle distributions					
stepPar <- c(1,	5,	0.5,	3) # step distribution parameters		
anglePar <- c(pi,	-0.3,	1,	8) # angle distribution parameters 		
# (angle mean of -0.3 in state 2	 to mimic the haggises' tendency to run 				
# clockwise around hills when active)					
# specify regression coefficients for the transition probabilities					
beta <- matrix(c(-3.5,	0.8,	 # intercept			
                 0.35,	-0.4,	 # slope			
                 -0.01,	0.01,	 # slope^2			
                 0,	0),	 # temp (no effect)			
               nrow=4,	byrow=TRUE)				
# simulate wild haggis movement data (15 haggises)					
dataraw <- moveHMM::simData(nbAnimals=15,	nbStates=2,	stepDist="gamma",	angleDist="vm"	,	
                   stepPar=stepPar,	anglePar=anglePar,	beta=beta,			
                   covs=covs,	obsPerAnimal=400,	states=TRUE)			
# only keep relevant columns (animals' ID	 x/y coordinates	 and covariates)			
rawHaggis <- dataraw[,	c(1,	4,	5,	6,	8)]
# save simulated data into csv file					
write.csv(rawHaggis,	file="rawHaggises.csv",	row.names=FALSE)			
