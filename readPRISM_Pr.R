#remove all data/variables and plots
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
source("latlong2county.R")

# 1981 to 2012 PRISM precip data

#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/ppt.day.PRISM.AN81d.2p5min.1981.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon,each=dim_lat), y = rep(rev(grid_lat),dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

# which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]
countynum<-length(countytoget)
Precipitation_month<-matrix(NA,nrow=countynum,ncol=32*12) # 32 years monthly county precip
Precipitation_GS<-matrix(NA,nrow=countynum,ncol=32) #Cumulative Pr in the growing season

SF<-read.table("LastSpringFrost") #SF[i.j]: ith county jth year
FF<-read.table("FirstFallFrost")
countymeanSF<-floor(rowMeans(SF))
countymeanFF<-floor(rowMeans(FF))
for (i in 1:32){
  #read data from 1981 to 2012
  precip_file<-paste("/gpfs/group/kzk10/default/private/data_archive/old_PRISM/ppt.day.PRISM.AN81d.2p5min.",i+1980,".nc",sep="")
  prismpr<-nc_open(precip_file)
  if (i%%4==0){ #when exist Feb29
    k=61
    days=c(0,31,29,31,30,31,30,31,31,30,31,30,31)
  }else {
    k=60
    days=c(0,31,28,31,30,31,30,31,31,30,31,30,31)
  }
  Pr<-ncvar_get(prismpr,varid = "ppt",start = c(600,1,1),count = c(dim_lon-599,dim_lat,k+305)) 
  
  #size too large, start from 600(lon), which is exactly 100 deg W
  for (n in 1:countynum){ 
    #determine which grids to get based on which counties we need
    gridstoget<-which(countyname==countytoget[n],arr.ind = T)
    gridnum<-dim(gridstoget)[1]
    GS<-c(k:(k+184))
    Prgrid<-matrix(NA,nrow=gridnum,ncol=k+305)
    for (j in 1:gridnum){
      #Prgrid is the daily data of Pr in given county given year
      Prgrid[j, ]<-Pr[gridstoget[j,2]-599,dim_lat+1-gridstoget[j,1], ] #dimension=(lon,lat), and lat order is reversed
    }
    Prgrid<-Prgrid[complete.cases(Prgrid), ]
    Prgrid<-colMeans(Prgrid) #county mean data in each day
    Precipitation_GS[n,i]<-sum(Prgrid[GS])
    for (m in 1:12){
      daystoget<-days[m]
      monthindex<-(i-1)*12+m
      Precipitation_month[n,monthindex]<-sum(Prgrid[c((sum(days[1:m])+1):(sum(days[1:(m+1)])))])
    }
    print(n)
  }
  print(i)
}
write.table(Precipitation_GS,paste("Pr"))
for (n in 1:countynum){
  write.table(Precipitation_GS[n, ],paste("Data/PRISM/Pr/Pr_PRISM_",countytoget[n],sep = ""), sep=" ")
}
for (n in 1:countynum){
  write.table(Precipitation_month[n, ],paste("Data/PRISM/Pr_monthly/Pr_monthly_PRISM_",countytoget[n],sep = ""), sep=" ")
}
write.table(c(1981:2012),"Data/PRISM/Year_PRISM",sep=" ")