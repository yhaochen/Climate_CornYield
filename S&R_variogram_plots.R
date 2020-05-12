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
library(fields)
load("S&Rimplementation/Prism_S&Rmodel")
#variogram using fields
meanres<-rep(NA,length(levels(Data_prism$fips)))
lat<-rep(NA,length(levels(Data_prism$fips)))
lon<-rep(NA,length(levels(Data_prism$fips)))
GDD<-rep(NA,length(levels(Data_prism$fips)))
EDD<-rep(NA,length(levels(Data_prism$fips)))
Pr<-rep(NA,length(levels(Data_prism$fips)))
fip<-rep(NA,length(levels(Data_prism$fips)))
for (i in 1:length(levels(Data_prism$fips))){
  index<-which(Data_prism$fips==levels(Data_prism$fips)[i])
  if (length(index)!=0){
    meanres[i]<-mean(Data_prism$resfix[index])
    lat[i]<-Data_prism$lat[index[1]]
    lon[i]<-Data_prism$lon[index[1]]
    GDD[i]<-mean(Data_prism$GDD[index])
    EDD[i]<-mean(Data_prism$EDD[index])
    Pr[i]<-mean(Data_prism$Pr[index])
    fip[i]<-levels(Data_prism$fips)[i]
  }
}
countymeanres<-data.frame(res=meanres,lon=lon,lat=lat,fips=fip)
countymeanres$londis<-111*(countymeanres$lon+82)*cos(countymeanres$lat/180*pi)
countymeanres$latdis<-111*(countymeanres$lat-37)
countymeanres<-countymeanres[complete.cases(countymeanres), ]

#fields package
variog<-vgram(loc = matrix(data=c(countymeanres$londis,countymeanres$latdis),nrow=1821,ncol=2), y = countymeanres$res, 
              dmax = 300, N =60,type=c("correlogram"))
plot(variog, breaks = pretty(variog$d, 60, eps.correct = 1), add=FALSE,pch=20,lwd=2,cex.lab=1.4,cex.axis=1.4,xlab="Distance (km)",ylab="Semivariance",main="")
boxplotVGram(variog, N=30, breaks = pretty(variog$d, 50, eps.correct = 1), plot=TRUE,pch=20,cex.lab=1.4,cex.axis=1.4,xlab="Distance (km)",ylab="Semivariance",main="")

variog<-vgram(loc = matrix(data=c(countymeanres$londis,countymeanres$latdis),nrow=1821,ncol=2), y = countymeanres$res,
              dmax = 300, N =30,type=c("cross-correlogram"))
plot(variog, breaks = pretty(variog$d, 30, eps.correct = 1), add=FALSE,pch=20,lwd=2,cex.lab=1.4,cex.axis=1.4,
     xlab="Distance (km)",ylab="Correlation coefficient",main="")
boxplotVGram(variog, N=60, breaks = pretty(variog$d, 60, eps.correct = 1), plot=TRUE)



#range in each direction

coordinates(countymeanres)<-~londis+latdis
range<-rep(NA,60)
for (i in 1:60){
  angle<-i*3
  V<-variogram(res~1,locations=countymeanres$londis+countymeanres$latdis,data=countymeanres,width=10,cutoff=500,alpha=angle,tol.hor=15)
  Vfit=fit.variogram(V, vgm(c("Sph")), fit.sills = TRUE, fit.ranges = TRUE)
  range[i]<-Vfit$range[2]
  #plot(V,Vfit,pch=20,col="black",xlab="Distance (km)",ylab="Semivariance",cex.lab=1.4,cex.axis=1.4)
}
plot(seq(3,180,3),range,type="l",xlab="Degrees (0=north,90=east)",ylab="Range of fitted variogram (km)",cex.lab=1.4,cex.axis=1.4)

#scatter plot
countyavg<-data.frame(res=meanres,GDD=GDD,EDD=EDD,Pr=Pr,fips=fip,lon=lon,lat=lat)
countyavg<-countyavg[complete.cases(countyavg), ]
plot(countyavg$Pr,countyavg$res,pch=20,xlab="County Precipitation mean (mm)",ylab="County residual mean (bush/acre)",cex.lab=1.4,cex.axis=1.4)


#each year's variogram
Yeartoget=1
index<-which(Data_prism$year==Yeartoget)
Datatoget<-Data_prism[index,]
Datatoget$londis<-111*(Datatoget$lon+82)*cos(Datatoget$lat/180*pi)
Datatoget$latdis<-111*(Datatoget$lat-37)
coordinates(Datatoget)<-~londis+latdis
V<-variogram(resfix~1,locations=Datatoget$londis+Datatoget$latdis,data=Datatoget,width=5,cutoff=300)
Vfit=fit.variogram(V, vgm(c("Sph")), fit.sills = TRUE, fit.ranges = TRUE)
plot(V,Vfit,pch=20,col="black",xlab="Distance (km)",ylab="Semivariance",cex.lab=1.4,cex.axis=1.4)
