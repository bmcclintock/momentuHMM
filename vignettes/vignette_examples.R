library(sp)
library(doParallel)

example_wd <- ("~/Documents/Dropbox/current projects/moveHMM extension/momentuHMM/vignette examples/")

append.RData <- function(x, file) {
  old.objects <- load(file)
  save(list = c(old.objects, deparse(substitute(x))), file = file, version=2) # version=2 covers R 1.40 to 3.5.0; if version>2 then produces CRAN warning about depenency on R >= 3.5.0
}

###################################################
### Elephant example
###################################################
#source(paste0(getwd(),"/examples/elephantExample.R"))
load(paste0(example_wd,"elephantExample.RData"))

activityBudgets<-table(viterbi(m3))/nrow(m3$data)

save(activityBudgets,file=paste0(getwd(),"/vignette_inputs.RData"))

#latlongDat<-as.data.frame(elephantData)
#sp::coordinates(latlongDat)<-c("x","y")
#sp::proj4string(latlongDat) = sp::CRS("+proj=utm +zone=30N +datum=WGS84" )
#latlongDat<-sp::spTransform(latlongDat,sp::CRS("+proj=longlat +datum=WGS84"))
#satPlotStates<-plotSat(as.data.frame(latlongDat),zoom=8,location=c(median(rawData$lon),median(rawData$lat)),ask=FALSE,return=TRUE,states=viterbi(m3),col=c("#009E73", "#F0E442"),stateNames=c("encamped","exploratory"))

acfLag<-24*14
pdf(file=paste0(getwd(),"/plot_elephantResults%03d.pdf"),onefile=FALSE,width=5,height=5)
plot(m3,ask=FALSE,plotCI=TRUE,covs=data.frame(hour=12))
acf(pr[["stepRes"]], lag.max = acfLag, na.action = na.pass, xlab="Lag (hours)", main = "")
acf(elephantData$step,na.action=na.pass,lag=acfLag,main="",xlab="Lag (hours)")
dev.off()

for(plt in seq(1,17)[-c(1,2,9,10,12,13,15,16,17)])
  unlink(paste0("plot_elephantResults0",ifelse(plt>9,"","0"),plt,".pdf"))

png(filename="elephant_plotSat.png",width=5,height=5,units="in",res=80)
plotSat(data.frame(x=rawData$lon,y=rawData$lat),zoom=8,location=c(median(rawData$lon),median(rawData$lat)),ask=FALSE)
dev.off()

rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### Northern fur seal example
###################################################
#source(paste0(getwd(),"/examples/nfsExample.R"))
load(paste0(example_wd,"nfsExample.RData"))

nfsTimeInStates<-nfsFits$miSum$Par$timeInStates

append.RData(nfsTimeInStates,file=paste0(getwd(),"/vignette_inputs.RData"))

pdf(file=paste0(getwd(),"/plot_nfsResults.pdf"),width=5,height=5)
plot(nfsFits,ask=FALSE,legend.pos="topright")
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### Turtle example
###################################################
#source(paste0(getwd(),"/examples/turtleExample.R"))
load(paste0(example_wd,"turtleExample.RData"))

turtle_miSum<-turtleFits$miSum[c("Par","data")]

append.RData(turtle_miSum,file=paste0(getwd(),"/vignette_inputs.RData"))

library(TeachingDemos)
library(raster)

speedrast<-speedBrick[["X2012.12.02"]]
speedrast<-projectRaster(speedrast,crs="+init=epsg:26718 +proj=utm +zone=18 +datum=NAD27 +units=km +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat")
dirrast<-dirBrick[["X2012.12.02"]]
dirrast<-projectRaster(dirrast,crs="+init=epsg:26718 +proj=utm +zone=18 +datum=NAD27 +units=km +no_defs +ellps=clrk66 +nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat")

rastxy<-coordinates(speedrast)
rastx0<-unique(rastxy[,1])
rasty0<-unique(rastxy[,2])
rastx<-rep(rastx0,length(rasty0))
rasty<-rep(rasty0,each=length(rastx0))

x<-turtleFits$miSum$data$x/1000
y<-turtleFits$miSum$data$y/1000

pdf(file=paste0(getwd(),"/plot_turtleResults%03d.pdf"),onefile=FALSE)
plot(turtleFits,plotCI=TRUE,covs=data.frame(angle_osc=1),ask=FALSE)
dev.off()

for(plt in seq(1,6)[-c(1,2,4)])
  unlink(paste0("plot_turtleResults0",ifelse(plt>9,"","0"),plt,".pdf"))


pdf(file=paste0(getwd(),"/plot_turtleResults2.pdf"),width=7,height=7*50/73)
par(mar=c(4,4,0,1))
raster::plot(speedrast,col=gray.colors(20, start=1, end=0.3),xlab="easting (km)",ylab="northing (km)",legend.args=list(text='       speed (m/s)', side=3, line=1),box=FALSE,bty="n")
my.symbols(rastx,rasty,ms.arrows, xsize=20,ysize=20,add=TRUE,angle=getValues(dirrast), r=getValues(speedrast), length=.015)
points(x,y,type="o",pch=20,col=c("#E69F00", "#56B4E9")[turtleFits$miSum$Par$states],cex=.5)
segments(x0=x[-length(x)],y0=y[-length(y)],x1=x[-1],y1=y[-1],
         col=c("#E69F00", "#56B4E9")[turtleFits$miSum$Par$states][-length(turtleFits$miSum$Par$states)],lwd=1.3)
#points(coordinates(turtleData)/1000,cex=0.15)
for(i in 1:length(turtleFits$miSum$errorEllipse))
  lines(turtleFits$miSum$errorEllipse[[i]]/1000,col=adjustcolor(c("#E69F00", "#56B4E9")[turtleFits$miSum$Par$states][i],alpha.f=0.25),cex=0.6)
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### Grey seal example
###################################################
#source(paste0(getwd(),"/examples/greySealExample.R"))
load(paste0(example_wd,"greySealResults.RData"))

greySealTimeInStates<-greySealPool$Par$timeInStates

append.RData(greySealTimeInStates,file=paste0(getwd(),"/vignette_inputs.RData"))

pdf(file=paste0(getwd(),"/plot_greySealResults%03d.pdf"),onefile=FALSE)
plot(greySealPool,plotCI=TRUE,ask=FALSE)
dev.off()

for(plt in seq(1,18)[-c(2,6,9,13)])
  unlink(paste0("plot_greySealResults0",ifelse(plt>9,"","0"),plt,".pdf"))

load(paste0(getwd(),"/greySealData.RData"))
load(paste0(getwd(),"/coastUTM.RData"))

colvect<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

png(file=paste0(getwd(),"/plot_greySealResults1.png"),width=7.25,height=5,units="in",res=90)
plot(greySealPool$data$x,greySealPool$data$y,type="o", asp=1, pch=20,ylab="Latitude",xlab="Longitude",ylim=c(6060000,6289000),xlim=c(513000,925000),xaxt="n",yaxt="n",xaxs="i",yaxs="i")
axis(1,labels=c(-2,0,2,4),at=c(565312,696098,(434680-312928+696098),(565330-312928+696098)))
axis(2,labels=c(55,56,57),at=c(6094820,6206100,6317500))

points(greySealData$x,greySealData$y,pch=20)
cx=1.6

for(j in 1:5) {
  points(greySealPool$data$x[greySealPool$Par$states==j],greySealPool$data$y[greySealPool$Par$states==j],col=colvect[j],asp=1, pch=20,cex=cx)
}
points(greySealPool$data$x[length(greySealPool$data$x)],greySealPool$data$y[length(greySealPool$data$y)],col=colvect[greySealPool$Par$states[length(greySealPool$Par$states)]],asp=1, pch=20,cex=cx)

for(i in 1:(length(greySealPool$data$x)-1)) {
  polygon(greySealPool$errorEllipse[[i]],lty=2,border=rgb(red=col2rgb(colvect,alpha=FALSE)[1,greySealPool$Par$states[i]]/255,green=col2rgb(colvect,alpha=FALSE)[2,greySealPool$Par$states[i]]/255,blue=col2rgb(colvect,alpha=FALSE)[3,greySealPool$Par$states[i]]/255,alpha=0.25))
}
polygon(greySealPool$errorEllipse[[length(greySealPool$errorEllipse)]],lty=2,border=rgb(red=col2rgb(colvect,alpha=FALSE)[1,greySealPool$Par$states[length(greySealPool$Par$states)]]/255,green=col2rgb(colvect,alpha=FALSE)[2,greySealPool$Par$states[length(greySealPool$Par$states)]]/255,blue=col2rgb(colvect,alpha=FALSE)[3,greySealPool$Par$states[length(greySealPool$Par$states)]]/255,alpha=0.25))

#for(i in 1:nrow(greySealPool$Par$stateProbs$est)){
#  tmp<-sort(greySealPool$Par$stateProbs$est[i,])
#  for(j in 1:(5-1)){
#    if(tmp[j]>0.05){
#      ind<-order(greySealPool$Par$stateProbs$est[2,])[j]
#      points(greySealPool$data$x[i],greySealPool$data$y[i],col=colvect[ind],asp=1, pch=20,cex=cx/(j+1))
#    }
#  }
#}

leg=character(5)
leg[1]="Abertay haul-out state"
leg[2]="Farne Islands haul-out state"
leg[3]="Dogger Bank foraging state"
leg[4]="Low-speed exploratory state"
leg[5]="High-speed exploratory state"

legend("topright",leg,pch=rep(20,length(leg)),col=colvect,bty="n")

lines(coastUTM$Easting, coastUTM$Northing)

#add estimated locations for centres of attraction
points(centers,col="#D55E00",pch=20,cex=cx)
dev.off()

png(file=paste0(getwd(),"/plot_greySealResults2.png"),width=7.25,height=5,units="in",res=90)
plot(greySealSim$x,greySealSim$y,type="o", asp=1, pch=20,ylab="Latitude",xlab="Longitude",ylim=c(6060000,6289000),xlim=c(513000,925000),xaxt="n",yaxt="n",xaxs="i",yaxs="i")
points(greySealSim$x,greySealSim$y,col=colvect[greySealSim$states],pch=20)
axis(1,labels=c(-2,0,2,4),at=c(565312,696098,(434680-312928+696098),(565330-312928+696098)))
axis(2,labels=c(55,56,57),at=c(6094820,6206100,6317500))
legend("topright",leg,pch=rep(20,length(leg)),col=colvect,bty="n")
lines(coastUTM$Easting, coastUTM$Northing)
points(centers,col="#D55E00",pch=20,cex=1.6)
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### Elephant seal example
###################################################
#source(paste0(getwd(),"/examples/sesExample.R"))
load(paste0(example_wd,"sesExample.RData"))

sesCIbeta<-m3$CIbeta

append.RData(tracks,file=paste0(getwd(),"/vignette_inputs.RData"))
append.RData(sesCIbeta,file=paste0(getwd(),"/vignette_inputs.RData"))

pdf(file=paste0(getwd(),"/plot_sesResults%03d.pdf"),onefile=FALSE)
plot(m3,plotCI=TRUE,ask=FALSE)
dev.off()

for(plt in seq(1,27)[-c(1,3,7,9)])
  unlink(paste0("plot_sesResults0",ifelse(plt>9,"","0"),plt,".pdf"))


library(maps) # for map plots
library(mapdata) # for map plots
library(marmap) # to plot bathymetry
library(sp) # for degAxis
png(file=paste0(getwd(),"/plot_sesResults2.png"),width=12,height=6,units="in",res=80)
plot(data$x, data$y, axes=FALSE, xlab="longitude", ylab="latitude", col="white")
degAxis(1)
degAxis(2)
map('worldHires', add=TRUE, fill=TRUE, col='white')

pal <- c("#78c679","#F0E442","#E31A1C","#88419d")

for(id in unique(data$ID)) {
  ind <- which(data$ID==id)
  segments(x0 = data$x[ind[-length(ind)]], y0 = data$y[ind[-length(ind)]], x1 = data$x[ind[-1]], y1 = data$y[ind[-1]], 
           col = pal[states[ind[-length(ind)]]], lwd=2)
}

legend("topleft", legend = stateNames,
       col = pal, lwd=2, bg="white")
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### harbour seal example
###################################################
#source(paste0(getwd(),"/examples/harbourSealExample.R"))
load(paste0(example_wd,"harbourSealExample.RData"))

append.RData(hsActivityBudgets,file=paste0(getwd(),"/vignette_inputs.RData"))

options(scipen=10)
cx=0.9
tmpData<-as.data.frame(hsData)
coordinates(tmpData)<-c("x","y")
proj4string(tmpData)<-CRS("+init=epsg:27700 +units=m")
tmpData<-spTransform(tmpData,CRS="+proj=longlat +ellps=WGS84")

png(file=paste0(getwd(),"/plot_harbourSeal.png"),width=7.25,height=5,units="in",res=84)
par(mfrow=c(1,2))
for(id in c(8,1)){
  freqs<-miSum.ind$Par$states[hsData$ID==id]
  x<-tmpData$x[which(hsData$ID==id)]
  y<-tmpData$y[which(hsData$ID==id)]
  errorEllipse<-miSum.ind$errorEllipse[which(hsData$ID==id)]
  errorEllipse<-lapply(errorEllipse,function(x){x<-as.data.frame(x);
                                                coordinates(x)<-c("x","y");
                                                proj4string(x)<-CRS("+init=epsg:27700 +units=m");
                                                x<-spTransform(x,CRS="+proj=longlat +ellps=WGS84");
                                                as.data.frame(x)})
  
  #png(file=paste0(getwd(),"/plot_harbourSeal",id,".png"),width=5,height=7.25,units="in",res=90)
  plot(x,y,type="o",pch=20,cex=.5,xlab="longitude",ylab="latitude",cex.lab=0.8,cex.axis=.7)
  points(x_ti[which(data$ID==id)],y_ti[which(data$ID==id)],pch=20,type="o",cex=0.5,col=gray(.5))
  map('worldHires', c('UK', 'Ireland', 'Isle of Man','Isle of Wight'), col=gray(0.85),fill=T,add=T)
  map.axes(cex.axis=.7)
  points(x,y,type="o",pch=20,col=c("#E69F00", "#56B4E9", "#009E73")[freqs],cex=.5)
  segments(x0=x[-length(x)],y0=y[-length(y)],x1=x[-1],y1=y[-1],
           col=c("#E69F00", "#56B4E9", "#009E73")[freqs][-length(freqs)],lwd=1.3)
  #points(coordinates(turtleData)/1000,cex=0.15)
  for(i in 1:length(errorEllipse))
    lines(errorEllipse[[i]],col=adjustcolor(c("#E69F00", "#56B4E9", "#009E73")[freqs][i],alpha.f=0.25),cex=0.6)
}
legend("topright",c(stateNames,"observed"),col=c("#E69F00", "#56B4E9", "#009E73",gray(.5)),pch=20,cex=.7)
dev.off()

pdf(file=paste0(getwd(),"/plot_harbourSealResults%03d.pdf"),onefile=FALSE)
plot(miSum.ind,plotCI=TRUE,ask=FALSE)
dev.off()

for(plt in seq(1,37)[-c(20)])
  unlink(paste0("plot_harbourSealResults0",ifelse(plt>9,"","0"),plt,".pdf"))

rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### northern fulmar example
###################################################
#source(paste0(getwd(),"/examples/northernFulmarExample.R"))
load(paste0(example_wd,"northernFulmarExample.RData"))
fulmarElapsedTime <- round(m2$mod$elapsedTime/60,0)
timeIn1 <- timeInStates(m2)
timeIn2 <- timeInStates(m2,by="birdID")
append.RData(fulmarElapsedTime,file=paste0(getwd(),"/vignette_inputs.RData"))
append.RData(timeIn1,file=paste0(getwd(),"/vignette_inputs.RData"))
append.RData(timeIn2,file=paste0(getwd(),"/vignette_inputs.RData"))
#append.RData(m2,file=paste0(getwd(),"/vignette_inputs.RData"))
#append.RData(newProj,file=paste0(getwd(),"/vignette_inputs.RData"))
png(file=paste0(getwd(),"/plot_northernFulmarExample.png"),width=7.25,height=5,units="in",res=80)
plotSat(m2,zoom=7,shape=c(17,1,17,1,17,1),size=2,col=rep(c("#E69F00", "#56B4E9", "#009E73"),each=2),stateNames=c("sea ARS","sea Transit","boat ARS","boat Transit","colony ARS","colony Transit"),projargs=newProj,ask=FALSE)
dev.off()

rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### group dynamic example
###################################################
#source(paste0(getwd(),"/examples/groupExample.R"))
load(paste0(example_wd,"groupExample.RData"))
pdf(file=paste0(getwd(),"/plot_groupExampleCentroid%03d.pdf"),onefile=FALSE)
plot(centroidData,ask=FALSE)
dev.off()
unlink(paste0("plot_groupExampleCentroid002.pdf"))

pdf(file=paste0(getwd(),"/plot_groupExample%03d.pdf"),onefile=FALSE)
plot(groupData,compact=TRUE,ask=FALSE)
dev.off()

for(plt in seq(1,nbAnimals+1)[-1])
  unlink(paste0("plot_groupExample0",ifelse(plt>9,"","0"),plt,".pdf"))

pdf(file=paste0(getwd(),"/plot_groupExampleResults%03d.pdf"),onefile=FALSE)
plot(groupFit,ask=FALSE)
dev.off()

for(plt in seq(1,nbAnimals+3)[-c(7,8)])
  unlink(paste0("plot_groupExampleResults0",ifelse(plt>9,"","0"),plt,".pdf"))

###################################################
### pilot whale example
###################################################
#source(paste0(getwd(),"/examples/pilotWhaleExample.R"))
load(paste0(example_wd,"pilotWhaleExample.RData"))
fitmix1_Par <- getPar(fitmix1)
append.RData(fitmix1_Par,file=paste0(getwd(),"/vignette_inputs.RData"))
Par0_mix2 <- getPar0(fitmix1,mixtures=2)
Par0_mix2$beta$beta[1,] <- c(-2.26, -3.93, -0.58, 
                              0.03, -2.25, -0.26, 
                             -3.38, -4.79, -2.82, 
                             -1.06,  -3.3, -3.43)
Par0_mix2$beta$beta[2,] <- c(-2.51, -3.32, -2.63, 
                              0.03, -1.26, -0.12, 
                             -96.8, -3.62, -1.75, 
                             -1.76, -2.14, -1.38)
Par0_mix2$beta$pi <- c(0.73, 0.27)
append.RData(Par0_mix2,file=paste0(getwd(),"/vignette_inputs.RData"))
fitmix2_Par <- getPar(fitmix2)
append.RData(fitmix2_Par,file=paste0(getwd(),"/vignette_inputs.RData"))
fitmix3_Par <- getPar(fitmix3)
append.RData(fitmix3_Par,file=paste0(getwd(),"/vignette_inputs.RData"))
fitmix4_Par <- getPar(fitmix4)
append.RData(fitmix4_Par,file=paste0(getwd(),"/vignette_inputs.RData"))
fitfix_Par <- getPar(fitfix)
append.RData(fitfix_Par,file=paste0(getwd(),"/vignette_inputs.RData"))
#pilotWhaleAIC <- AIC(fitmix1,fitmix2,fitmix3,fitfix)
#append.RData(pilotWhaleAIC,file=paste0(getwd(),"/vignette_inputs.RData"))
#pilotWhaleAICweights <- AICweights(fitmix1,fitmix2,fitmix3,fitfix)
#append.RData(pilotWhaleAICweights,file=paste0(getwd(),"/vignette_inputs.RData"))
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### harbor porpoise HHMM example
###################################################
#source(paste0(getwd(),"/examples/harborPorpoiseExample.R"))
load(paste0(example_wd,"harborPorpoiseExample.RData"))
hhmmPar <- getPar(hhmm)
append.RData(hhmmPar,file=paste0(getwd(),"/vignette_inputs.RData"))
pdf(file=paste0(getwd(),"/plot_harborPorpoiseStates%03d.pdf"),width=8,height=11,onefile=FALSE)
plotStates(hhmm,ask=FALSE)
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### garter snake HHMM example
###################################################
#source(paste0(getwd(),"/examples/garterSnakeExample.R"))
load(paste0(example_wd,"garterSnakeExample.RData"))
hhmm2Par <- getPar(hhmm)
append.RData(hhmm2Par,file=paste0(getwd(),"/vignette_inputs.RData"))
pdf(file=paste0(getwd(),"/plot_garterSnakePR.pdf"),width=5,height=2.5)
plotPR(hhmm)
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### atlantic cod HHMM example
###################################################
#source(paste0(getwd(),"/examples/codExample.R"))
load(paste0(example_wd,"codExample.RData"))
hhmm3Par <- getPar(hhmm)
append.RData(hhmm3Par,file=paste0(getwd(),"/vignette_inputs.RData"))
pdf(file=paste0(getwd(),"/plot_codExample%03d.pdf"),width=5,height=5,onefile = FALSE)
plot(hhmm, ask=FALSE)
dev.off()
pdf(file=paste0(getwd(),"/plot_codStationary%03d.pdf"),width=5,height=5,onefile = FALSE)
plotStationary(hhmm, plotCI=TRUE)
dev.off()

for(plt in c(3,4,5,6))
  unlink(paste0("plot_codExample00",plt,".pdf"))

rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### horn shark HHMM example
###################################################
#source(paste0(getwd(),"/examples/hornSharkExample.R"))
load(paste0(example_wd,"hornSharkExample.RData"))
hhmm4Par <- getPar(hhmm)
append.RData(hhmm4Par,file=paste0(getwd(),"/vignette_inputs.RData"))
trProbs12 <- getTrProbs(hhmm,covIndex=c(1,3))
append.RData(trProbs12,file=paste0(getwd(),"/vignette_inputs.RData"))
stats12 <- stationary(hhmm,covIndex=c(1,3))
append.RData(stats12,file=paste0(getwd(),"/vignette_inputs.RData"))
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### buffalo recharge example
###################################################
#source(paste0(getwd(),"/examples/buffaloExample.R"))
load(paste0(example_wd,"buffaloExample.RData"))
bufPar <- buffaloFits$miSum$Par$beta[c("mu","g0","theta")]
bufPar$timeInStates <- buffaloFits$miSum$Par$timeInStates
append.RData(bufPar,file=paste0(getwd(),"/vignette_inputs.RData"))

pdf(file=paste0(getwd(),"/plot_buffaloExample%03d.pdf"),width=8,height=3,onefile = FALSE)
plot(buffaloFits,plotCI=TRUE,legend.pos="bottom",ask=FALSE)
dev.off()

pdf(file="plot_buffaloStates.pdf",width=7.5,height=5)
plotSpatialCov(buffaloFits,dist2sabie)
dev.off()

for(plt in c(1:10,12))
  unlink(paste0("plot_buffaloExample0",ifelse(plt>9,"","0"),plt,".pdf"))

#trProbs <- getTrProbs(buffaloFits, getCI=TRUE)
# plot estimates and CIs for Pr(discharged) at each time step
pdf(file=paste0(getwd(),"/plot_buffaloResults.pdf"),width=8,height=3)
par(mar=c(5,4,4,2)-c(0,0,2,1)) # bottom, left, top, right
plot(trProbs$est[1,2,],type="l", 
     ylim=c(0,1), ylab="Pr(discharged)", xlab="t", col=c("#E69F00", "#56B4E9")[buffaloFits$miSum$Par$states])
arrows(1:dim(trProbs$est)[3],
       trProbs$lower[1,2,],
       1:dim(trProbs$est)[3],
       trProbs$upper[1,2,],
       length=0.025, angle=90, code=3, col=c("#E69F00", "#56B4E9")[buffaloFits$miSum$Par$states], lwd=1.3)
abline(h=0.5,lty=2)
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### land constraint example
###################################################
#source(paste0(getwd(),"/examples/landConstraintExample.R"))
load(paste0(example_wd,"landConstraintExample.RData"))
png(file=paste0(getwd(),"/plot_landConstraintExample.png"),width=6,height=6,units="in",res=80)
raster::plot(boundary$boundary,legend.width=1, legend.shrink=0.75,legend.args=list(text='            Distance to water', font=2, line=1, cex=0.8))
points(simBound$mu.x,simBound$mu.y,type="l")
dev.off()
rm(list=ls()[-which(ls()=="example_wd")])


# reduce size of png files
system(paste("pngquant -f --ext .png --speed=1 *.png"))

# reduce size of pdf files
system(paste0("$HOME/Documents/Dropbox/'current projects'/'moveHMM extension'/momentuHMM/momentuHMM/vignettes/./script.sh"))
