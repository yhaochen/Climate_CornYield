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

# 1981 to 2013 PRISM climate data

#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmax.day.PRISM.AN81d.2p5min.1981.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon,each=dim_lat), y = rep(rev(grid_lat),dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

#INPUT 1: which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)
GDD<-matrix(NA,nrow=countynum,ncol=33) #want to get countywise GDD and EDD for each year
EDD<-matrix(NA,nrow=countynum,ncol=33)
tmax<-matrix(NA,nrow=countynum,ncol=33) #max min temperature
tmin<-matrix(NA,nrow=countynum,ncol=33)
SF<-read.table("LastSpringFrost") #SF[i.j]: ith county jth year
FF<-read.table("FirstFallFrost")
countymeanSF<-floor(rowMeans(SF))
countymeanFF<-floor(rowMeans(FF))



for (i in 1:32){
  #read data from 1981 to 2013
  Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmax.day.PRISM.AN81d.2p5min.",i+1980,".nc",sep="")
  Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/tmin.day.PRISM.AN81d.2p5min.",i+1980,".nc",sep="")
  mettmax<-nc_open(Tmax_file)
  mettmin<-nc_open(Tmin_file)
  k=60
  if (i%%4==0){ #when exist Feb29
    k=61
  }
  Tmax<-ncvar_get(mettmax,varid = "tmax",start = c(600,1,k),count = c(dim_lon-599,dim_lat,184)) #data of Tmax in each year
  Tmin<-ncvar_get(mettmin,varid = "tmin",start = c(600,1,k),count = c(dim_lon-599,dim_lat,184))
  #data size is too large, start from 600(lon), which is exactly 100 deg W
  for (n in 1:countynum){ 
    #determine which grids to get based on which counties we need, and calculate county average Tmax,Tmin,GDD,EDD,Spring/fall frost
    gridstoget<-which(countyname==countytoget[n],arr.ind = T)
    gridnum<-dim(gridstoget)[1]
    GS<-c((countymeanSF[n]-k+1):(countymeanFF[n]-k+1))
    Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=length(GS))
    Tmingrid<-matrix(NA,nrow=gridnum,ncol=length(GS))
    for (j in 1:gridnum){
    #Tmaxgrid is the daily data of Tmax in given county given year
    Tmaxgrid[j, ]<-Tmax[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1],GS] #dimension=(lon,lat), and lat order is reversed
    Tmingrid[j, ]<-Tmin[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1],GS]
    }
    Tmaxgrid<-Tmaxgrid[complete.cases(Tmaxgrid), ]
    Tmingrid<-Tmingrid[complete.cases(Tmingrid), ]
    gridnum<-dim(Tmaxgrid)[1]
    #Tdailygrid<-(Tmaxgrid+Tmingrid)/2
    #frostthres[n,i]<-quantile(density(Tdailygrid),probs=0.1)
    
    tmax[n,i]<-mean(Tmaxgrid,NA.rm=TRUE)
    
    #tmaxdaily<-colMeans(Tmaxgrid)
    
    tmin[n,i]<-mean(Tmingrid,NA.rm=TRUE)
    #tmindaily<-colMeans(Tmingrid)
    #frostdays<-which(tmindaily<0)+k-1
    #frostdaysTem<-tmindaily[frostdays-k+1]
    #frost<-matrix(data=c(frostdays,frostdaysTem),nrow=length(frostdays),ncol=2)
    #write.table(frost,paste("Data/PRISM/Frost/Frost_PRISM_",countytoget[n],i,sep = ""), sep=" ")

    GDD[n,i]<-GDDEDD(Tmaxgrid,Tmingrid)[1]/gridnum
    EDD[n,i]<-GDDEDD(Tmaxgrid,Tmingrid)[2]/gridnum
    print(n)
  }
  print(i)
}

for (n in 1:countynum){
  write.table(GDD[n, ],paste("Data/PRISM/GDD/GDD_PRISM_",countytoget[n],sep = ""), sep=" ")
  write.table(EDD[n, ],paste("Data/PRISM/EDD/EDD_PRISM_",countytoget[n],sep = ""), sep=" ")
  write.table(tmax[n, ],paste("Data/PRISM/Tmax/EDD_PRISM_",countytoget[n],sep = ""), sep=" ")
  write.table(tmin[n, ],paste("Data/PRISM/Tmin/EDD_PRISM_",countytoget[n],sep = ""), sep=" ")
}
write.table(c(1981:2013),"Data/PRISM/Year_PRISM",sep=" ")





