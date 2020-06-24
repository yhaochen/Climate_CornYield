#Read each climate projection
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(usmap)
library(ncdf4)
library(precintcon)
library(BMS)
library(housingData)
library(binaryLogic)
source("latlong2county.R")
source("GDDEDD.R")
#Determine which county each grid is in s
file<-nc_open("/gpfs/group/kzk10/default/public/METDATA/raw/tmmx_1980.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon,each=dim_lat), y = rep(grid_lat,dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

#which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

ANSI<-load("ANSI")
CountyANSI<-load("CountyANSI")
StateANSI<-load("StateANSI")
ma <- function(arr, n){  #moving avg function
  res = rep(NA,length(arr)-n+1)
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n+1):i])
  }
  res
}

su <- function(arr, n){  #cumulative sum
  res = rep(NA,length(arr)-n+1)
  for(i in n:length(arr)){
    res[i] = sum(arr[n:i])
  }
  res
}


#13 models
modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM")

#Metdata observation
load("Metdata_temp/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
Data$year=Data$year+1978
Data$year<-factor(Data$year)

#with each climate projection, save the hind/proj yield data of 64 structures
for (q in 1:length(modelnames)){
  foldername<-paste("Data/MACAv2-METDATA/",modelnames[q],"_proj/",sep="")
  load(paste(foldername,"tmax",sep=""))
  load(paste(foldername,"tmin",sep=""))
  load(paste(foldername,"pr",sep=""))
  load(paste(foldername,"rhmax",sep=""))
  load(paste(foldername,"rhmin",sep=""))
  tmax=tmax-273.15
  tmin=tmin-273.15
  totdays<-sum(!is.na(tmax[1, ]))
  #calculate VPD,GDD,EDD
  RHmean<-(rhmax+rhmin)/2
  Tmean<-(tmax+tmin)/2
  e_s<-6.112*exp((17.269*Tmean)/(Tmean-237.3)) #saturated vapor pressure (hPa)
  VPD<-e_s*(1-RHmean/100)
  GDD<-matrix(NA,nrow=countynum,ncol=totdays)
  EDD<-matrix(NA,nrow=countynum,ncol=totdays)
  for (i in 1: countynum){
    for (j in 1: totdays){
      GDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[1]
      EDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[2]
    }
  }
  #find growing season range for each county each year
  Tthres<-10
  GS_start<-matrix(NA,nrow=countynum,ncol=94)
  GS_end<-matrix(NA,nrow=countynum,ncol=94)
  m1=1
  m2=1
  for (i in 1:countynum){
    m1=1
    m2=1
    for (j in 1:94){
      k=365
      if (j%%4==3){ 
        k=366
      }
      m2<-m1+k-1
      Tmaxtoget<-tmax[i,c(m1:m2)]
      Tmintoget<-tmin[i,c(m1:m2)]
      Tmeantoget<-0.5*(Tmaxtoget+Tmintoget)
      TMA<-ma(Tmeantoget,21)
      GS_start[i,j]<-which(TMA >= Tthres)[1]
      GS_end[i,j]<-GS_start[i,j]+184
      m1=m2+1
    }
  }
  
  #create the climmate projection dataframe
  Tmax_GS<-rep(NA,countynum*94)
  Tmin_GS<-rep(NA,countynum*94)
  GDD_GS<-rep(NA,countynum*94)
  EDD_GS<-rep(NA,countynum*94)
  VPD_GS<-rep(NA,countynum*94)
  Pr_GS<-rep(NA,countynum*94)
  GS_length<-rep(NA,countynum*94)
  
  for (i in 1:countynum){
    m1=1
    m2=1
    for (j in 1:94){
      k=365
      if (j%%4==3){ #when exist Feb29
        k=366
      }
      m2<-m1+k-1
      GS_length[(i-1)*94+j]<-GS_end[i,j]-GS_start[i,j]
      index<-c(m1:m2)
      if (!is.na(GS_start[i,j])){
        GSindex<-index[GS_start[i,j]:GS_end[i,j]]
        Tmax_GS[(i-1)*94+j]<-mean(tmax[i,GSindex])
        Tmin_GS[(i-1)*94+j]<-mean(tmin[i,GSindex])
        GDD_GS[(i-1)*94+j]<-sum(GDD[i,GSindex])
        EDD_GS[(i-1)*94+j]<-sum(EDD[i,GSindex])
        VPD_GS[(i-1)*94+j]<-sum(VPD[i,GSindex])
        Pr_GS[(i-1)*94+j]<-sum(pr[i,GSindex])
      }
      m1=m2+1
    }
  }
  Data_proj<-data.frame(StateANSI=rep(StateANSI,each=94),countyANSI=rep(CountyANSI,each=94),fips=rep(ANSI,each=94),
                        year=rep(c(1:94),countynum),Tmax_GS=Tmax_GS,Tmin_GS=Tmin_GS,GDD_GS=GDD_GS,GS_length=GS_length,
                        EDD_GS=EDD_GS,VPD_GS=VPD_GS,Pr_GS=Pr_GS)
  save(Data_proj,file=paste("Metdata_temp/macaprojdataframe/Data_",modelnames[q],sep=""))
}

