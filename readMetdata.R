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
countytoget<-countyNames[c(960,1281:1367,1742:1745)]
countynum<-length(countytoget)
GDD<-matrix(NA,nrow=countynum,ncol=38) #want to get countywise GDD and EDD for each year
EDD<-matrix(NA,nrow=countynum,ncol=38)
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
    Tmaxgrid<-Tmax[gridstoget[ ,1],gridstoget[ ,2]-500, ]-273.15
    Tmaxgrid<-na.omit(apply(Tmaxgrid,3,diag))
    Tmingrid<-Tmin[gridstoget[ ,1],gridstoget[ ,2]-500, ]-273.15
    Tmingrid<-na.omit(apply(Tmingrid,3,diag))
    gridnum<-dim(Tmaxgrid)[1]
    GDD[n,i]<-GDDEDD(Tmaxgrid,Tmingrid)[1]/gridnum
    EDD[n,i]<-GDDEDD(Tmaxgrid,Tmingrid)[2]/gridnum
    print(n)
  }
  print(i)
}
for (n in 1:countynum){
  write.table(GDD[n, ],paste("Data/MET/GDD/GDD_MET_",countytoget[n],sep = ""), sep=" ")
  write.table(EDD[n, ],paste("Data/MET/EDD/EDD_MET_",countytoget[n],sep = ""), sep=" ")
}
write.table(c(1979:2016),"Data/MET/Year_MET",sep=" ")





