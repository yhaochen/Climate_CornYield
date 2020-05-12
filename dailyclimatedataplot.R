#plot dayly climate data
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
source("latlong2county.R")
source("GDDEDD.R")

#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmax.day.PRISM.AN81d.2p5min.1981.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon,each=dim_lat), y = rep(rev(grid_lat),dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)




#Statetoget=PA #2212-2278, 42
stateind<-c(756,854)

Tmax_day<-matrix(NA,nrow=stateind[2]+1-stateind[1],ncol=365*32)
Tmin_day<-matrix(NA,nrow=stateind[2]+1-stateind[1],ncol=365*32)
Precipitation_day<-matrix(NA,nrow=stateind[2]+1-stateind[1],ncol=365*32)

for (yeartoget in 1981:2012){
  Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmax.day.PRISM.AN81d.2p5min.",yeartoget,".nc",sep="")
  Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmin.day.PRISM.AN81d.2p5min.",yeartoget,".nc",sep="")
  Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/ppt.day.PRISM.AN81d.2p5min.",yeartoget,".nc",sep="")
  
  tmax<-nc_open(Tmax_file)
  tmin<-nc_open(Tmin_file)
  pr<-nc_open(Pr_file)
  
  Tmax<-ncvar_get(tmax,varid = "tmax",start = c(600,1,1),count = c(dim_lon-599,dim_lat,365))
  Tmin<-ncvar_get(tmin,varid = "tmin",start = c(600,1,1),count = c(dim_lon-599,dim_lat,365))
  Pr<-ncvar_get(pr,varid = "ppt",start = c(600,1,1),count = c(dim_lon-599,dim_lat,365))
  for (n in 1:(stateind[2]+1-stateind[1])){ 
    #determine which grids to get based on which counties we need
    gridstoget<-which(countyname==countyNames[n-1+stateind[1]],arr.ind = T)
    gridnum<-dim(gridstoget)[1]
    Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=365)
    Tmingrid<-matrix(NA,nrow=gridnum,ncol=365)
    Prgrid<-matrix(NA,nrow=gridnum,ncol=365)
    for (j in 1:gridnum){
      Tmaxgrid[j, ]<-Tmax[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ]
      Tmingrid[j, ]<-Tmin[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ]
      Prgrid[j, ]<-Pr[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ] #dimension=(lon,lat), and lat order is reversed
    }
    index<-c(((yeartoget-1980)*365-364):((yeartoget-1980)*365))
    Tmax_day[n, index]<-colMeans(Tmaxgrid)
    Tmin_day[n, index]<-colMeans(Tmingrid)
    Precipitation_day[n, index]<-colMeans(Prgrid)
    print(n)
  }
  print(yeartoget)
}
countymeanTmax<-matrix(NA,nrow=stateind[2]+1-stateind[1],ncol=365)
countymeanTmin<-matrix(NA,nrow=stateind[2]+1-stateind[1],ncol=365)
countymeanPr<-matrix(NA,nrow=stateind[2]+1-stateind[1],ncol=365)

for (i in 1:(stateind[2]+1-stateind[1])){
  for (j in 1:365){
    index<-seq(j,j+365*31,by=365)
    countymeanTmax[i,j]<-mean(Tmax_day[i, index],na.rm=TRUE)
    countymeanTmin[i,j]<-mean(Tmin_day[i, index],na.rm=TRUE)
    countymeanPr[i,j]<-mean(Precipitation_day[i, index],na.rm=TRUE)
  }
}
countyTmax_anomaly<-Tmax_day-matrix(rep(countymeanTmax,32),nrow=stateind[2]+1-stateind[1],ncol=365*32)
countyTmin_anomaly<-Tmin_day-matrix(rep(countymeanTmin,32),nrow=stateind[2]+1-stateind[1],ncol=365*32)
countyPr_anomaly<-Precipitation_day-matrix(rep(countymeanPr,32),nrow=stateind[2]+1-stateind[1],ncol=365*32)
stateTmax_anomaly<-colMeans(countyTmax_anomaly)
stateTmin_anomaly<-colMeans(countyTmin_anomaly)
statePr_anomaly<-colMeans(countyPr_anomaly)

index<-c((365*(2002-1980)-364):(365*(2002-1980)))
plot(c(1:365),stateTmax_anomaly[index],type="l",xlab="days of the year",ylab="Tmax anomaly (deg C)",cex.lab=1.4,cex.axis=1.4)
plot(c(1:365),stateTmin_anomaly[index],type="l",xlab="days of the year",ylab="Tmin anomaly (deg C)",cex.lab=1.4,cex.axis=1.4)
plot(c(1:365),statePr_anomaly[index],type="l",xlab="days of the year",ylab="Pr anomaly (mm)",cex.lab=1.4,cex.axis=1.4)
