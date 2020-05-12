#Plot of variables in the data (prism)
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
library(ggplot2)
library(usmap)
library(precintcon)
library(BMS)
library(adaptMCMC)
library(geoR) #for variogram
library(housingData) #for returning county centroid
library(gstat)
library(spatstat)
#1. Tmax/Tmin distribution for all counties !
load("Data_prism")
first10yr_index<-which(Data_prism$year<=10)
last10yr_index<-which(Data_prism$year>=23)
par(mar=c(5,5,3,2))
plot(density(Data_prism$Tmax[first10yr_index]),col="blue",xlab="Daily maximum temperature (°C)",ylab="Density",main="Daily maximum temperature 
     distribution in eastern US (Mar-Aug)",cex.axis=1.5,cex.lab=1.5)
lines(density(Data_prism$Tmax[last10yr_index]),col="red")
abline(v=29,lty=2)
legend("topleft",col=c("blue","red"),lty = c(1,1),legend = c("1981-1990","2003-2012"),cex=1.5, bty = "n")


a=CDF(density(Data_prism$Tmax[first10yr_index]))
b=CDF(density(Data_prism$Tmax[last10yr_index]))
time=seq(29,33.6,by=0.05)
par(mar=c(5,5,3,2))
plot(time,1-a(time),type="l",lwd=2,log="y",col="blue",xlab="Daily maximum temperature (°C)",ylab="Survival function",cex.axis=1.5,cex.lab=1.5)
lines(time,1-b(time),col="red",lwd=2)
legend(x=32,y=0.3,col=c("blue","red"),lty = c(1,1),lwd=c(2,2),legend = c("1981-1990","2003-2012"),cex=1.5, bty = "n")


time=seq(13,34,by=0.1)
par(mar=c(5,5,3,2))
plot(time,a(time),type="l",lwd=2,col="blue",xlab="Daily maximum temperature (°C)",ylab="CDF",cex.axis=1.5,cex.lab=1.5)
lines(time,b(time),col="red",lwd=2)
legend("topleft",col=c("blue","red"),lty = c(1,1),lwd=c(2,2),legend = c("1981-1990","2003-2012"),cex=1.5, bty = "n")
abline(v=29,lty=2)


par(mar=c(5,5,3,2))
plot(density(Data_prism$Tmin[first10yr_index]),col="blue",xlab="Daily minimum temperature (°C)",ylab="Density",main="Daily minimum temperature 
     distribution in eastern US (Mar-Aug)",cex.axis=1.5,cex.lab=1.5)
lines(density(Data_prism$Tmin[last10yr_index]),col="red")
abline(v=0,lty=2)
legend("topleft",col=c("blue","red"),lty = c(1,1),legend = c("1981-1990","2003-2012"),cex=1.5, bty = "n")


a(29)

b(29)