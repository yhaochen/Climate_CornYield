#climate proj

#This file reads GDD and EDD data in PRISM grids

#remove all data/variables and plots
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
source("latlong2county.R")
source("GDDEDD.R")

# 2086 to 2090 MACAv2 METDATA

#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_GFDL-ESM2M_r1i1p1_rcp45_2081_2085_CONUS_daily.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon-360,each=dim_lat), y = rep(rev(grid_lat),dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

#INPUT 1: which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)
GDD<-matrix(NA,nrow=countynum,ncol=5) #want to get countywise GDD and EDD for each year
EDD<-matrix(NA,nrow=countynum,ncol=5)
Precip<-matrix(NA,nrow=countynum,ncol=5)
tmax<-matrix(NA,nrow=countynum,ncol=5) #max min temperature
tmin<-matrix(NA,nrow=countynum,ncol=5)



  Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_GFDL-ESM2M_r1i1p1_rcp45_2081_2085_CONUS_daily.nc")
  Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmin_GFDL-ESM2M_r1i1p1_rcp45_2081_2085_CONUS_daily.nc")
  mettmax<-nc_open(Tmax_file)
  mettmin<-nc_open(Tmin_file)
for (i in 1:5){
  k=60
  if (i>=4){ #when exist Feb29
    k=61
  }
  Tmax<-ncvar_get(mettmax,varid = "air_temperature",start = c(596,1,k+365*(i-1)),count = c(dim_lon-595,dim_lat,184)) 
  Tmin<-ncvar_get(mettmin,varid = "air_temperature",start = c(596,1,k+365*(i-1)),count = c(dim_lon-595,dim_lat,184))
  #data size is too large, start from 596(lon), which is exactly 100 deg W
  for (n in 1:countynum){ 
    #determine which grids to get based on which counties we need, and calculate county average Tmax,Tmin,GDD,EDD,Spring/fall frost
    gridstoget<-which(countyname==countytoget[n],arr.ind = T) # （lat,lon）
    gridnum<-dim(gridstoget)[1]
    Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=184)
    Tmingrid<-matrix(NA,nrow=gridnum,ncol=184)
    for (j in 1:gridnum){
      #Tmaxgrid is the daily data of Tmax in given county given year
      Tmaxgrid[j, ]<-Tmax[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ] #dimension=(lon,lat,time), reverse order of lat because the first row is the northest lat
      Tmingrid[j, ]<-Tmin[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ]
    }
    Tmaxgrid<-Tmaxgrid[complete.cases(Tmaxgrid), ]
    Tmingrid<-Tmingrid[complete.cases(Tmingrid), ]
    gridnum<-dim(Tmaxgrid)[1]
    tmax[n,i]<-mean(Tmaxgrid,na.rm=TRUE)-273.15
    tmin[n,i]<-mean(Tmingrid,na.rm=TRUE)-273.15
    #GDD[n,i]<-GDDEDD(Tmaxgrid-273.15,Tmingrid-273.15)[1]/gridnum
    #EDD[n,i]<-GDDEDD(Tmaxgrid-273.15,Tmingrid-273.15)[2]/gridnum
    print(n)
  }
}

for (n in 1:countynum){
  #write.table(GDD[n, ],paste("Data/MACAv2-METDATA/GDD_GFDL-ESM2M_rcp45_2081_2085/GDD_",countytoget[n],sep = ""), sep=" ")
  #write.table(EDD[n, ],paste("Data/MACAv2-METDATA/EDD_GFDL-ESM2M_rcp45_2081_2085/EDD_",countytoget[n],sep = ""), sep=" ")
  write.table(tmax[n, ],paste("Data/MACAv2-METDATA/Tmax_GFDL-ESM2M_rcp45_2081_2085/Tmax_",countytoget[n],sep = ""), sep=" ")
  write.table(tmin[n, ],paste("Data/MACAv2-METDATA/Tmin_GFDL-ESM2M_rcp45_2081_2085/Tmin_",countytoget[n],sep = ""), sep=" ")
}

  
  GDD<-matrix(NA,nrow=countynum,ncol=5) #want to get countywise GDD and EDD for each year
  EDD<-matrix(NA,nrow=countynum,ncol=5)
  Precip<-matrix(NA,nrow=countynum,ncol=5)
  tmax<-matrix(NA,nrow=countynum,ncol=5) #max min temperature
  tmin<-matrix(NA,nrow=countynum,ncol=5)
  for (n in 1:countynum){
    GDD[n, ]<-read.table(paste("Data/MACAv2-METDATA/GDD_GFDL-ESM2M_rcp45_2081_2085/GDD_",countytoget[n],sep = ""))$x
    EDD[n, ]<-read.table(paste("Data/MACAv2-METDATA/EDD_GFDL-ESM2M_rcp45_2081_2085/EDD_",countytoget[n],sep = ""))$x
    Precip[n, ]<-read.table(paste("Data/MACAv2-METDATA/Pr_GFDL-ESM2M_rcp45_2081_2085/Pr_",countytoget[n],sep = ""))$x
    tmax[n, ]<-read.table(paste("Data/MACAv2-METDATA/Tmax_GFDL-ESM2M_rcp45_2081_2085/Tmax_",countytoget[n],sep = ""))$x
    tmin[n, ]<-read.table(paste("Data/MACAv2-METDATA/Tmin_GFDL-ESM2M_rcp45_2081_2085/Tmin_",countytoget[n],sep = ""))$x
  }
  
  Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_pr_GFDL-ESM2M_r1i1p1_rcp45_2081_2085_CONUS_daily.nc")
  metpr<-nc_open(Pr_file)
  for (i in 1:5){
    k=60
    if (i>=4){ #when exist Feb29
      k=61
    }
    Pr<-ncvar_get(metpr,varid = "precipitation",start = c(596,1,k+365*(i-1)),count = c(dim_lon-595,dim_lat,184)) 
    #data size is too large, start from 596(lon), which is exactly 100 deg W
    for (n in 1:countynum){ 
      gridstoget<-which(countyname==countytoget[n],arr.ind = T) # （lat,lon）
      gridnum<-dim(gridstoget)[1]
      Prgrid<-matrix(NA,nrow=gridnum,ncol=184)
      for (j in 1:gridnum){
        Prgrid[j, ]<-Pr[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ] #dimension=(lon,lat,time), reverse order of lat because the first row is the northest lat
      }
      Prgrid<-Prgrid[complete.cases(Prgrid), ]
      Precip[n,i]<-sum(colMeans(Prgrid,na.rm=TRUE))
      print(n)
    }
  }
  
  for (n in 1:countynum){
    write.table(Precip[n, ],paste("Data/MACAv2-METDATA/Pr_GFDL-ESM2M_rcp45_2081_2085/Pr_",countytoget[n],sep = ""), sep=" ")
  }  
  
  countynames<-countytoget
  countynames[25]<-"alabama,dekalb"
  countynames[59]<-"alabama,st. clair"
  countynames[135]<-"arkansas,st. francis"
  countynames[188]<-"georgia,dekalb"
  countynames[323]<-"illinois,dekalb"
  countynames[326]<-"illinois,dupage"
  countynames[353]<-"illinois,lasalle"
  countynames[391]<-"illinois,st. clair"
  countynames[421]<-"indiana,dekalb"
  countynames[450]<-"indiana,laporte"
  countynames[480]<-"indiana,st. joseph"
  countynames[569]<-"iowa,o'brien"
  countynames[761]<-"louisiana,st. bernard"
  countynames[762]<-"louisiana,st. charles"
  countynames[763]<-"louisiana,st. helena"
  countynames[764]<-"louisiana,st. james"
  countynames[765]<-"louisiana,st. john the baptist"
  countynames[766]<-"louisiana,st. landry"
  countynames[767]<-"louisiana,st. martin"
  countynames[768]<-"louisiana,st. mary"
  countynames[769]<-"louisiana,st. tammany"
  countynames[798]<-"maryland,prince george's"
  countynames[799]<-"maryland,queen anne's"
  countynames[801]<-"maryland,st. mary's"
  countynames[882]<-"michigan,st. clair"
  countynames[883]<-"michigan,st. joseph"
  countynames[960]<-"minnesota,st. louis"
  countynames[992]<-"mississippi,desoto"
  countynames[1089]<-"missouri,dekalb"
  countynames[1155]<-"missouri,st. charles"
  countynames[1156]<-"missouri,st. clair"
  countynames[1157]<-"missouri,st. francois"
  countynames[1158]<-"missouri,st. louis"
  countynames[1159]<-"missouri,st. louis city"
  countynames[1160]<-"missouri,ste. genevieve"
  countynames[1243]<-"new york,st. lawrence"
  countynames[1576]<-"tennessee,dekalb"
  countynames[1693]<-"virginia,hampton city"
  countynames[1715]<-"virginia,newport news city"
  countynames[1716]<-"virginia,norfolk city"
  countynames[1741]<-"virginia,suffolk city"
  countynames[1745]<-"virginia,virginia beach city"
  countynames[1866]<-"wisconsin,st. croix"
  #convert county names to county ANSI code, a five-digit code, used to match the data in yielddata
  ANSI<-rep(NA,countynum)
  for (i in 1:countynum){
    state_county<-strsplit(countynames[i],",")
    State<-state_county[[1]][1]
    County<-state_county[[1]][2]
    if (State=="louisiana"){
      County<-paste(County,"parish",sep=" ")
    }
    ANSI[i]<-fips(State, County)
  }

load("Data_prism")
GDD_mean<-rep(NA,countynum)
EDD_mean<-rep(NA,countynum)
Pr_mean<-rep(NA,countynum)
Tmax_mean<-rep(NA,countynum)
Tmin_mean<-rep(NA,countynum)
for (i in 1:countynum){ #1981-2012 average annual GDD,EDD
  index<-which(Data_prism$fips==levels(Data_prism$fips)[i])
  GDD_mean[i]<-mean(Data_prism$GDD[index],na.rm=TRUE)
  EDD_mean[i]<-mean(Data_prism$EDD[index],na.rm=TRUE)
  Pr_mean[i]<-mean(Data_prism$Pr[index],na.rm=TRUE)
  Tmax_mean[i]<-mean(Data_prism$Tmax[index],na.rm=TRUE)
  Tmin_mean[i]<-mean(Data_prism$Tmin[index],na.rm=TRUE)
}

Data<-data.frame(fips=ANSI,GDD=rowMeans(GDD)-GDD_mean,EDD=rowMeans(EDD)-EDD_mean,Pr=rowMeans(Precip)-Pr_mean,
                 Tmax=rowMeans(tmax)-Tmax_mean,Tmin=rowMeans(tmin)-Tmin_mean)
a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data, values = "GDD") + labs(title = "GDD anomaly of each county between 2081-2085 and 1981-2012")+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-350,350),midpoint=0,name="GDD (deg C)")+theme(plot.title = element_text(size=14))
plot(a)
b=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data, values = "EDD") + labs(title = "EDD anomaly of each county between 2081-2085 and 1981-2012")+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-350,350),midpoint=0,name="EDD (deg C)")+theme(plot.title = element_text(size=14))
plot(b)
c=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data, values = "Pr") + labs(title = "Pr anomaly of each county between 2081-2085 and 1981-2012")+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-200,200),midpoint=0,name="Pr (mm)")+theme(plot.title = element_text(size=14))
plot(c)
d=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data, values = "Tmax") + labs(title = "Tmax anomaly of each county between 2081-2085 and 1981-2012")+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-5,5),midpoint=0,name="Tmax (deg C)")+theme(plot.title = element_text(size=14))
plot(d)
e=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data, values = "Tmin") + labs(title = "Tmin anomaly of each county between 2081-2085 and 1981-2012")+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-5,5),midpoint=0,name="Tmin (deg C)")+theme(plot.title = element_text(size=14))
plot(e)

Data$temp<-(Data$Tmax+Data$Tmin)/2
plot((Data$Tmax+Data$Tmin)/2,Data$Pr,pch=20,cex=0.5,xlab="Temperature change (deg C)",ylab="Precipitation change (mm)",main="",cex.axis=1.3,cex.lab=1.3)
abline(lm(Pr~temp, data = Data),lty=2,col="red")

plot(density((Tmax_mean+Tmin_mean)/2,na.rm=TRUE),xlim=c(7,30),xlab="Growing season mean temperature (deg C)",ylab="PDF",main="" ,col="blue",lwd=2,cex.axis=1.3,cex.lab=1.3)
lines(density(rowMeans(tmax+tmin)/2,na.rm=TRUE),lwd=2,col="red")
legend("topleft",col=c("blue","red"),lty = c(1,1),lwd=c(2,2),legend = c("1981-2012","2081-2085"),cex=1.3, bty = "n")
abline(v=mean(Tmax_mean+Tmin_mean,na.rm=TRUE)/2,lty=2,col="blue")
abline(v=mean(rowMeans(tmax+tmin),na.rm=TRUE)/2,lty=2,col="red")

plot(density(GDD_mean,na.rm=TRUE),xlim=c(700,3000),xlab="GDD (deg C)",ylab="PDF",main="" ,col="blue",lwd=2,cex.axis=1.3,cex.lab=1.3)
lines(density(rowMeans(GDD),na.rm=TRUE),lwd=2,col="red")
legend("topleft",col=c("blue","red"),lty = c(1,1),lwd=c(2,2),legend = c("1981-2012","2081-2085"),cex=1.3, bty = "n")
abline(v=mean(GDD_mean,na.rm=TRUE),lty=2,col="blue")
abline(v=mean(rowMeans(GDD),na.rm=TRUE),lty=2,col="red")

plot(density(EDD_mean,na.rm=TRUE),xlim=c(0,550),xlab="EDD (deg C)",ylab="PDF",main="" ,col="blue",lwd=2,cex.axis=1.3,cex.lab=1.3)
lines(density(rowMeans(EDD),na.rm=TRUE),lwd=2,col="red")
legend("topright",col=c("blue","red"),lty = c(1,1),lwd=c(2,2),legend = c("1981-2012","2081-2085"),cex=1.3, bty = "n")
abline(v=mean(EDD_mean,na.rm=TRUE),lty=2,col="blue")
abline(v=mean(rowMeans(EDD),na.rm=TRUE),lty=2,col="red")

plot(density(Pr_mean,na.rm=TRUE),xlim=c(300,1000),ylim=c(0,0.008),xlab="Growing season precipitation (mm)",ylab="PDF",main="" ,col="blue",lwd=2,cex.axis=1.3,cex.lab=1.3)
lines(density(rowMeans(Precip),na.rm=TRUE),lwd=2,col="red")
legend("topleft",col=c("blue","red"),lty = c(1,1),lwd=c(2,2),legend = c("1981-2012","2081-2085"),cex=1.3, bty = "n")
abline(v=mean(Pr_mean,na.rm=TRUE),lty=2,col="blue")
abline(v=mean(rowMeans(Precip),na.rm=TRUE),lty=2,col="red")