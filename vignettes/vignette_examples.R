library(sp)

append.RData <- function(x, file) {
  old.objects <- load(file)
  save(list = c(old.objects, deparse(substitute(x))), file = file)
}

###################################################
### Elephant example
###################################################
#source(paste0(getwd(),"/elephantExample.R"))
load(paste0(example_wd,"elephantExample.RData"))

activityBudgets<-table(viterbi(m3))/nrow(m3$data)

save(activityBudgets,file=paste0(getwd(),"/vignette_inputs.RData"))

#latlongDat<-as.data.frame(elephantData)
#sp::coordinates(latlongDat)<-c("x","y")
#sp::proj4string(latlongDat) = sp::CRS("+proj=utm +zone=30N +datum=WGS84" )
#latlongDat<-sp::spTransform(latlongDat,sp::CRS("+proj=longlat +datum=WGS84"))
#satPlotStates<-plotSat(as.data.frame(latlongDat),zoom=8,location=c(median(rawData$lon),median(rawData$lat)),ask=FALSE,return=TRUE,states=viterbi(m3),col=c("#009E73", "#F0E442"),stateNames=c("encamped","exploratory"))

acfLag<-24*14
pdf(file=paste0(getwd(),"/plot_elephantResults%03d.pdf"),onefile=FALSE)
plot(m3,ask=FALSE,plotCI=TRUE,covs=data.frame(hour=12))
acf(pr[["stepRes"]], lag.max = acfLag, na.action = na.pass, xlab="Lag (hours)", main = "")
acf(elephantData$step,na.action=na.pass,lag=acfLag,main="",xlab="Lag (hours)")
dev.off()

for(plt in seq(1,17)[-c(1,2,9,10,12,13,15,16,17)])
  unlink(paste0("plot_elephantResults0",ifelse(plt>9,"","0"),plt,".pdf"))

png(filename="elephant_plotSat.png",width=6,height=6,units="in",res=90)
plotSat(data.frame(x=rawData$lon,y=rawData$lat),zoom=8,location=c(median(rawData$lon),median(rawData$lat)),ask=FALSE)
dev.off()

rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### Northern fur seal example
###################################################
#source(paste0(getwd(),"/nfsExample.R"))
load(paste0(example_wd,"nfsExample.RData"))

nfsTimeInStates<-nfsFits$miSum$Par$timeInStates

append.RData(nfsTimeInStates,file=paste0(getwd(),"/vignette_inputs.RData"))

pdf(file=paste0(getwd(),"/plot_nfsResults.pdf"))
plot(nfsFits,ask=FALSE)
dev.off()
rm(list=ls()[-which(ls()=="example_wd" | ls()=="append.RData")])

###################################################
### Turtle example
###################################################
#source(paste0(getwd(),"/turtleExample.R"))
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
plot(speedrast,col=gray.colors(20, start=1, end=0.3),xlab="easting (km)",ylab="northing (km)",legend.args=list(text='       speed (m/s)', side=3, line=1),box=FALSE,bty="n")
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
#source(paste0(getwd(),"/greySealExample_TPM.R"))
load(paste0(example_wd,"greySealResults_TPM.RData"))

greySealTimeInStates<-greySealPool$Par$timeInStates

append.RData(greySealTimeInStates,file=paste0(getwd(),"/vignette_inputs.RData"))

pdf(file=paste0(getwd(),"/plot_greySealResults%03d.pdf"),onefile=FALSE)
plot(greySealPool,plotCI=TRUE,ask=FALSE)
dev.off()

for(plt in seq(1,18)[-c(2,6,9,13)])
  unlink(paste0("plot_greySealResults0",ifelse(plt>9,"","0"),plt,".pdf"))

load(paste0(getwd(),"/greySealData_TPM.RData"))
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
#source(paste0(getwd(),"/sesExample.R"))
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
pdf(file=paste0(getwd(),"/plot_sesResults2.pdf"),width=12,height=6)
plot(data$x, data$y, axes=FALSE, xlab="longitude", ylab="latitude", col="white")
degAxis(1)
degAxis(2)
map('worldHires', add=TRUE, fill=TRUE, col='white')

for(id in unique(data$ID)) {
  ind <- which(data$ID==id)
  segments(x0 = data$x[ind[-length(ind)]], y0 = data$y[ind[-length(ind)]], x1 = data$x[ind[-1]], y1 = data$y[ind[-1]], 
           col = pal[states[ind[-length(ind)]]], lwd=2)
}

legend("topleft", legend = stateNames,
       col = pal, lwd=2, bg="white")
dev.off()
rm(list=ls()[-which(ls()=="example_wd")])




