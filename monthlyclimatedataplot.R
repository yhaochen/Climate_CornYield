#plot monthly climate data
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
source("latlong2county.R")
source("GDDEDD.R")

#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmax.mon.PRISM.AN81m.2p5min.1981.nc")
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

Tmax_month<-matrix(NA,nrow=countynum,ncol=12*32)
Tmin_month<-matrix(NA,nrow=countynum,ncol=12*32)
Precipitation_month<-matrix(NA,nrow=countynum,ncol=12*32)

for (yeartoget in 1981:2012){
Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmax.mon.PRISM.AN81m.2p5min.",yeartoget,".nc",sep="")
Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmin.mon.PRISM.AN81m.2p5min.",yeartoget,".nc",sep="")
Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/ppt.mon.PRISM.AN81m.2p5min.",yeartoget,".nc",sep="")
  
tmax<-nc_open(Tmax_file)
tmin<-nc_open(Tmin_file)
pr<-nc_open(Pr_file)
  
Tmax<-ncvar_get(tmax,varid = "tmax",start = c(600,1,1),count = c(dim_lon-599,dim_lat,12))
Tmin<-ncvar_get(tmin,varid = "tmin",start = c(600,1,1),count = c(dim_lon-599,dim_lat,12))
Pr<-ncvar_get(pr,varid = "ppt",start = c(600,1,1),count = c(dim_lon-599,dim_lat,12))
for (n in 1:countynum){ 
  #determine which grids to get based on which counties we need
  gridstoget<-which(countyname==countytoget[n],arr.ind = T)
  gridnum<-dim(gridstoget)[1]
  Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=12)
  Tmingrid<-matrix(NA,nrow=gridnum,ncol=12)
  Prgrid<-matrix(NA,nrow=gridnum,ncol=12)
  for (j in 1:gridnum){
    Tmaxgrid[j, ]<-Tmax[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ]
    Tmingrid[j, ]<-Tmin[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ]
    Prgrid[j, ]<-Pr[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ] #dimension=(lon,lat), and lat order is reversed
  }
  index<-c(((yeartoget-1980)*12-11):((yeartoget-1980)*12))
  Tmax_month[n, index]<-colMeans(Tmaxgrid)
  Tmin_month[n, index]<-colMeans(Tmingrid)
  Precipitation_month[n, index]<-colMeans(Prgrid)
  print(n)
}
print(yeartoget)
}

countymeanTmax<-matrix(NA,nrow=countynum,ncol=12)
countymeanTmin<-matrix(NA,nrow=countynum,ncol=12)
countymeanPr<-matrix(NA,nrow=countynum,ncol=12)
for (i in 1:countynum){
  for (j in 1:12){
    index<-seq(j,j+12*31,by=12)
    countymeanTmax[i,j]<-mean(Tmax_month[i, index],na.rm=TRUE)
    countymeanTmin[i,j]<-mean(Tmin_month[i, index],na.rm=TRUE)
    countymeanPr[i,j]<-mean(Precipitation_month[i, index],na.rm=TRUE)
  }
}
countyTmax_anomaly<-Tmax_month-matrix(rep(countymeanTmax,32),nrow=countynum,ncol=12*32)
countyTmin_anomaly<-Tmin_month-matrix(rep(countymeanTmin,32),nrow=countynum,ncol=12*32)
countyPr_anomaly<-Precipitation_month-matrix(rep(countymeanPr,32),nrow=countynum,ncol=12*32)

load("Data_prism")
fips<-levels(Data_prism$fips)
yeartoget<-2002
dir.create(paste("fig_weather/",yeartoget,sep=""))

for (i in 1:12){
  index<-(yeartoget-1981)*12+i
  tmaxfilename<-paste("fig_weather/",yeartoget,"/tmax_",yeartoget,"_",i,".jpeg",sep="")
  jpeg(file = tmaxfilename,width = 800,height=800)
  a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
               data=data.frame(tmax=countyTmax_anomaly[ ,index],fips=fips), values = "tmax") + labs(title = paste("Tmax anomaly of each county in ",yeartoget," ",monthtoget[i],sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-15,15),name="deg C")+theme(plot.title = element_text(size=14))
  plot(a)
  dev.off()
  tminfilename<-paste("fig_weather/",yeartoget,"/tmin_",yeartoget,"_",i,".jpeg",sep="")
  jpeg(file = tminfilename,width = 800,height=800)
  b=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
               data=data.frame(tmin=countyTmin_anomaly[ ,index],fips=fips), values = "tmin") + labs(title = paste("Tmin anomaly of each county in ",yeartoget," ",monthtoget[i],sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-15,15),name="deg C")+theme(plot.title = element_text(size=14))
  plot(b)
  dev.off()
  prfilename<-paste("fig_weather/",yeartoget,"/pr_",yeartoget,"_",i,".jpeg",sep="")
  jpeg(file = prfilename,width = 800,height=800)
  c=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
               data=data.frame(pr=countyPr_anomaly[ ,index],fips=fips), values = "pr") + labs(title = paste("Pr anomaly of each county in ",yeartoget," ",monthtoget[i],sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-120,120),name="mm")+theme(plot.title = element_text(size=14))
  plot(c)
  dev.off()
}


