#read projection daily macamet data: 2006-2099, save Tmax/min, Pr, RHmax/min
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
library(ggplot2)
library(usmap)
library(BMS)
library(housingData)
library(binaryLogic)
library(foreach)
source("latlong2county.R")
source("GDDEDD.R")

modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM",
              "MIROC-ESM-CHEM", "inmcm4", "HadGEM2-CC365", "CanESM2", "bcc-csm1-1-m")


file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_MIROC5_r1i1p1_rcp85_2006_2010_CONUS_daily.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon-360,each=dim_lat), y = rep(rev(grid_lat),dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

#which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)
years<-c(2006:2099)

cl<-makeCluster(10)
registerDoParallel(cl)

foreach (q = 1:18) %dopar%{
tmax<-matrix(NA,nrow=countynum,ncol=71*365+23*366)
tmin<-matrix(NA,nrow=countynum,ncol=71*365+23*366)
pr<-matrix(NA,nrow=countynum,ncol=71*365+23*366)
rhmax<-matrix(NA,nrow=countynum,ncol=71*365+23*366)
rhmin<-matrix(NA,nrow=countynum,ncol=71*365+23*366)


m1<-1 #used to record which columns to write in each loop
m2<-1
n1=1

dir.create(paste("/storage/work/h/hxy46/Countywise/SourceData/MACAv2-METDATA_proj_par/",modelnames[q],"_proj",sep=""),recursive = TRUE)

#read original data
for (i in 1:94){
  if ((i%%5==1)&(years[i]<2096)){
    Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_",modelnames[q],"_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmin_",modelnames[q],"_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_pr_",modelnames[q],"_r1i1p1_rcp85_"
                   ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    RHmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmax_",modelnames[q],"_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    RHmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmin_",modelnames[q],"_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    
    mettmax<-nc_open(Tmax_file)
    mettmin<-nc_open(Tmin_file)
    metpr<-nc_open(Pr_file)
    metrhmax<-nc_open(RHmax_file)
    metrhmin<-nc_open(RHmin_file)
    n1=1
  }
  
  if ((i%%5==1)&(years[i]>2095)){
    Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_",modelnames[q],"_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmin_",modelnames[q],"_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_pr_",modelnames[q],"_r1i1p1_rcp85_"
                   ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    RHmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmax_",modelnames[q],"_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    RHmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmin_",modelnames[q],"_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    
    mettmax<-nc_open(Tmax_file)
    mettmin<-nc_open(Tmin_file)
    metpr<-nc_open(Pr_file)
    metrhmax<-nc_open(RHmax_file)
    metrhmin<-nc_open(RHmin_file)
    n1=1
  }
  
  k=365
  if (i%%4==3){ #when exist Feb29
    k=366
  }
  
  m2<-m1+k-1
  #subset the grids to the east of 100W
  Tmax<-ncvar_get(mettmax,varid = "air_temperature",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  Tmin<-ncvar_get(mettmin,varid = "air_temperature",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  Pr<-ncvar_get(metpr,varid = "precipitation",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  RHmax<-ncvar_get(metrhmax,varid = "relative_humidity",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  RHmin<-ncvar_get(metrhmin,varid = "relative_humidity",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  
  n1<-n1+k
  #for each county find the grids within this county and take average
  for (n in 1:countynum){
    gridstoget<-which(countyname==countytoget[n],arr.ind = T)
    gridnum<-dim(gridstoget)[1]
    Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=k)
    Tmingrid<-matrix(NA,nrow=gridnum,ncol=k)
    Prgrid<-matrix(NA,nrow=gridnum,ncol=k)
    RHmaxgrid<-matrix(NA,nrow=gridnum,ncol=k)
    RHmingrid<-matrix(NA,nrow=gridnum,ncol=k)
    
    for (j in 1:gridnum){
      Tmaxgrid[j, ]<-Tmax[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ]
      Tmingrid[j, ]<-Tmin[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ]
      Prgrid[j, ]<-Pr[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ]
      RHmaxgrid[j, ]<-RHmax[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ]
      RHmingrid[j, ]<-RHmin[gridstoget[j,2]-595,dim_lat+1-gridstoget[j,1], ]
    }
    
    Tmaxgrid<-Tmaxgrid[complete.cases(Tmaxgrid), ]
    Tmingrid<-Tmingrid[complete.cases(Tmingrid), ]
    Prgrid<-Prgrid[complete.cases(Prgrid), ]
    RHmaxgrid<-RHmaxgrid[complete.cases(RHmaxgrid), ]
    RHmingrid<-RHmingrid[complete.cases(RHmingrid), ]
    
    tmax[n,c(m1:m2)]<-colMeans(Tmaxgrid)
    tmin[n,c(m1:m2)]<-colMeans(Tmingrid)
    pr[n,c(m1:m2)]<-colMeans(Prgrid)
    rhmax[n,c(m1:m2)]<-colMeans(RHmaxgrid)
    rhmin[n,c(m1:m2)]<-colMeans(RHmingrid)
    #print(n)
  }
  m1<-m2+1
}
save(pr,file=paste("SourceData/MACAv2-METDATA_proj_par/",modelnames[q],"_proj/pr",sep=""))
save(tmax,file=paste("SourceData/MACAv2-METDATA_proj_par/",modelnames[q],"_proj/tmax",sep=""))
save(tmin,file=paste("SourceData/MACAv2-METDATA_proj_par/",modelnames[q],"_proj/tmin",sep=""))
save(rhmax,file=paste("SourceData/MACAv2-METDATA_proj_par/",modelnames[q],"_proj/rhmax",sep=""))
save(rhmin,file=paste("SourceData/MACAv2-METDATA_proj_par/",modelnames[q],"_proj/rhmin",sep=""))
}