#remove all data/variables and plots
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
source("latlong2county.R")
source("GDDEDD.R")

# 1979 to 2016 MET climate data



#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/public/METDATA/raw/tmmx_1980.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon,each=dim_lat), y = rep(grid_lat,dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

#INPUT 1: which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)
GDD<-matrix(NA,nrow=countynum,ncol=38) #want to get countywise GDD and EDD for each year
EDD<-matrix(NA,nrow=countynum,ncol=38)
tmax<-matrix(NA,nrow=countynum,ncol=38)
tmin<-matrix(NA,nrow=countynum,ncol=38)
for (i in 1:38){
  #read data from 1979 to 2016
  Tmax_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/tmmx_",i+1978,".nc",sep="")
  Tmin_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/tmmn_",i+1978,".nc",sep="")
  mettmax<-nc_open(Tmax_file)
  mettmin<-nc_open(Tmin_file)
  k=60
  if (i%%4==2){ #when exist Feb29
    k=61
  }
  #determine which grids to get based on which counties we need, and calculate county average Tmax,Tmin
  Tmax<-ncvar_get(mettmax,varid = "air_temperature",start = c(1,501,k),count = c(dim_lat,dim_lon-500,184)) #annual mean for each grid
  Tmin<-ncvar_get(mettmin,varid = "air_temperature",start = c(1,501,k),count = c(dim_lat,dim_lon-500,184))
  for (n in 1:countynum){
    gridstoget<-which(countyname==countytoget[n],arr.ind = T)
    gridnum<-dim(gridstoget)[1]
    Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=184)
    Tmingrid<-matrix(NA,nrow=gridnum,ncol=184)
    for (j in 1:gridnum){
      #Tmaxgrid is the daily data of Tmax in given county given year
      Tmaxgrid[j, ]<-Tmax[gridstoget[j,1],gridstoget[j,2]-500, ]-273.15 #dimension=(lat,lon)
      Tmingrid[j, ]<-Tmin[gridstoget[j,1],gridstoget[j,2]-500, ]-273.15
      
    }
    Tmaxgrid<-Tmaxgrid[complete.cases(Tmaxgrid), ]
    Tmingrid<-Tmingrid[complete.cases(Tmingrid), ]
    gridnum<-dim(Tmaxgrid)[1]
    tmax[n,i]<-mean(Tmaxgrid,NA.rm=TRUE)
    tmin[n,i]<-mean(Tmingrid,NA.rm=TRUE)
    GDD[n,i]<-GDDEDD(Tmaxgrid,Tmingrid)[1]/gridnum
    EDD[n,i]<-GDDEDD(Tmaxgrid,Tmingrid)[2]/gridnum
    print(n)
  }
  print(i)
}
for (n in 1:countynum){
  write.table(GDD[n, ],paste("Data/MET/GDD/GDD_MET_",countytoget[n],sep = ""), sep=" ")
  write.table(EDD[n, ],paste("Data/MET/EDD/EDD_MET_",countytoget[n],sep = ""), sep=" ")
  write.table(tmax[n, ],paste("Data/MET/Tmax/Tmax_MET_",countytoget[n],sep = ""), sep=" ")
  write.table(tmin[n, ],paste("Data/MET/Tmin/Tmin_MET_",countytoget[n],sep = ""), sep=" ")
}
write.table(c(1979:2016),"Data/MET/Year_MET",sep=" ")





