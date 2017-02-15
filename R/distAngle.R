distAngle<-function(x,y,z){
  c(sqrt((y[1]-z[1])^2+(y[2]-z[2])^2),turnAngle(x,y,z))
}