#calculate the growing season and three growing phases from daily temperature data

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
library(housingData)
library(binaryLogic)
source("latlong2county.R")
source("GDDEDD.R")

#first read daily Tmax, Tmin, RHmax, RHmin, Pr

#yield data: 1981-2012 32yrs
#Metdata obs: 1979-2016
#MacaMet proj: 2066-2099
#
#
#
#
#
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
tmax<-matrix(NA,nrow=countynum,ncol=28*365+10*366)
tmin<-matrix(NA,nrow=countynum,ncol=28*365+10*366)
pr<-matrix(NA,nrow=countynum,ncol=28*365+10*366)
rhmax<-matrix(NA,nrow=countynum,ncol=28*365+10*366)
rhmin<-matrix(NA,nrow=countynum,ncol=28*365+10*366)

m1<-1 #used to record which columns to write in each loop
m2<-1

for (i in 1:38){
  #read data from 1979 to 2016
  Tmax_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/tmmx_",i+1978,".nc",sep="")
  Tmin_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/tmmn_",i+1978,".nc",sep="")
  Pr_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/pr_",i+1978,".nc",sep="")
  RHmax_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/rmax_",i+1978,".nc",sep="")
  RHmin_file<-paste("/gpfs/group/kzk10/default/public/METDATA/raw/rmin_",i+1978,".nc",sep="")
  
  mettmax<-nc_open(Tmax_file)
  mettmin<-nc_open(Tmin_file)
  metpr<-nc_open(Pr_file)
  metrhmax<-nc_open(RHmax_file)
  metrhmin<-nc_open(RHmin_file)
  
  k=365
  if (i%%4==2){ #when exist Feb29
    k=366
  }
  
  m2<-m1+k-1
  
  #determine which grids to get based on which counties we need, and calculate county average
  Tmax<-ncvar_get(mettmax,varid = "air_temperature",start = c(1,595,1),count = c(dim_lat,dim_lon-594,k)) 
  Tmin<-ncvar_get(mettmin,varid = "air_temperature",start = c(1,595,1),count = c(dim_lat,dim_lon-594,k))
  Pr<-ncvar_get(metpr,varid = "precipitation_amount",start = c(1,595,1),count = c(dim_lat,dim_lon-594,k))
  RHmax<-ncvar_get(metrhmax,varid = "relative_humidity",start = c(1,595,1),count = c(dim_lat,dim_lon-594,k))
  RHmin<-ncvar_get(metrhmin,varid = "relative_humidity",start = c(1,595,1),count = c(dim_lat,dim_lon-594,k))
  
  for (n in 1:countynum){
    gridstoget<-which(countyname==countytoget[n],arr.ind = T)
    gridnum<-dim(gridstoget)[1]
    Tmaxgrid<-matrix(NA,nrow=gridnum,ncol=k)
    Tmingrid<-matrix(NA,nrow=gridnum,ncol=k)
    Prgrid<-matrix(NA,nrow=gridnum,ncol=k)
    RHmaxgrid<-matrix(NA,nrow=gridnum,ncol=k)
    RHmingrid<-matrix(NA,nrow=gridnum,ncol=k)
    
    for (j in 1:gridnum){
      #Tmaxgrid is the daily data of Tmax in given county given year
      Tmaxgrid[j, ]<-Tmax[gridstoget[j,1],gridstoget[j,2]-594, ]-273.15 #dimension=(lat,lon)
      Tmingrid[j, ]<-Tmin[gridstoget[j,1],gridstoget[j,2]-594, ]-273.15
      Prgrid[j, ]<-Pr[gridstoget[j,1],gridstoget[j,2]-594, ]
      RHmaxgrid[j, ]<-RHmax[gridstoget[j,1],gridstoget[j,2]-594, ]
      RHmingrid[j, ]<-RHmin[gridstoget[j,1],gridstoget[j,2]-594, ]
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
    
    print(n)
  }
  m1<-m2+1
  print(i)
}

save(pr,file="Metdata_temp/Metpr_temp")
save(tmax,file="Metdata_temp/Mettmax_temp")
save(tmin,file="Metdata_temp/Mettmin_temp")
save(rhmax,file="Metdata_temp/Metrhmax_temp")
save(rhmin,file="Metdata_temp/Metrhmin_temp")



for (n in 1:countynum){
  write.table(tmax[n, ],paste("Data/MET/dailyTmax/Tmax_",countytoget[n],sep = ""), sep=" ")
  write.table(tmin[n, ],paste("Data/MET/dailyTmin/Tmin_",countytoget[n],sep = ""), sep=" ")
  write.table(pr[n, ],paste("Data/MET/dailyPr/Pr_",countytoget[n],sep = ""), sep=" ")
  write.table(rhmax[n, ],paste("Data/MET/dailyRHmax/RHmax_MET_",countytoget[n],sep = ""), sep=" ")
  write.table(rhmin[n, ],paste("Data/MET/dailyRHmin/RHmin_MET_",countytoget[n],sep = ""), sep=" ")
}


#second calculate VPD, GDD, EDD

load("Metdata_temp/Metpr")
load("Metdata_temp/Mettmax")
load("Metdata_temp/Mettmin")
load("Metdata_temp/Metrhmax")
load("Metdata_temp/Metrhmin")

RHmean<-(rhmax+rhmin)/2
Tmean<-(tmax+tmin)/2
e_s<-6.112*exp((17.269*Tmean)/(Tmean-237.3)) #saturated vapor pressure (hPa)
VPD<-e_s*(1-RHmean/100)

GDD<-matrix(NA,nrow=countynum,ncol=28*365+10*366)
EDD<-matrix(NA,nrow=countynum,ncol=28*365+10*366)
for (i in 1: countynum){
  for (j in 1: (28*365+10*366)){
    GDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[1]
    EDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[2]
  }
}
for (n in 1:countynum){
  write.table(VPD[n, ],paste("Data/MET/dailyVPD/VPD_",countytoget[n],sep = ""), sep=" ")
  write.table(GDD[n, ],paste("Data/MET/dailyGDD/GDD_",countytoget[n],sep = ""), sep=" ")
  write.table(EDD[n, ],paste("Data/MET/dailyEDD/EDD_",countytoget[n],sep = ""), sep=" ")
}
save(VPD,file="Metdata_temp/MetVPD")
save(GDD,file="Metdata_temp/MetGDD")
save(EDD,file="Metdata_temp/MetEDD")
load("Metdata_temp/MetVPD")
load("Metdata_temp/MetGDD")
load("Metdata_temp/MetEDD")
#third find growing season and three growing phases
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

Tthres<-10
GS_start<-matrix(NA,nrow=countynum,ncol=38)
GS_end<-matrix(NA,nrow=countynum,ncol=38)
m1=1
m2=1
for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:38){
    k=365
    if (j%%4==2){ #when exist Feb29
      k=366
    }
    m2<-m1+k-1
    Tmaxtoget<-tmax[i,c(m1:m2)]
    Tmintoget<-tmin[i,c(m1:m2)]
    Tmeantoget<-0.5*(Tmaxtoget+Tmintoget)
    GDDMA<-ma(Tmeantoget,21)
    GS_start[i,j]<-which(GDDMA >= Tthres)[1]
    GS_end[i,j]<-GS_start[i,j]+184
    m1=m2+1
  }
}
save(GS_start,file="Metdata_temp/MetGS_start")
save(GS_end,file="Metdata_temp/MetGS_end")


#fourth combine to a dataframe and regression
#dataframe should be countynum*yearnum
load("Metdata_temp/MetGS_start")
load("Metdata_temp/MetGS_end")
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#list of all US counties according to alphabetical order
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)

#states to include (eastern states)
statenames<-tolower(state.name[c(1,4,8,10,13:15,17,18,20,22:25,30,32,33,35,38,40,42,46,48,49)])
statenum<-length(statenames)
countynames<-rep(NA,length(countyNames))
countynum<-0
for (i in 1:length(countyNames)){
  state_county<-strsplit(countyNames[i],",")
  if (state_county[[1]][1] %in% statenames){
    countynum<-countynum+1 
    countynames[countynum]<-countyNames[i]
  }
}
countynames<-countynames[1:countynum]
originalcountynames<-countynames
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
CountyANSI<-as.numeric(substrRight(ANSI,3))
StateANSI<-as.numeric(substr(ANSI,1,2))


Tmax_GS<-rep(NA,countynum*38)
Tmin_GS<-rep(NA,countynum*38)
GDD_GS<-rep(NA,countynum*38)
EDD_GS<-rep(NA,countynum*38)
VPD_GS<-rep(NA,countynum*38)
Pr_GS<-rep(NA,countynum*38)
GS_length<-rep(NA,countynum*38)

for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:38){
    k=365
    if (j%%4==2){ #when exist Feb29
      k=366
    }
    m2<-m1+k-1
    GS_length[(i-1)*38+j]<-GS_end[i,j]-GS_start[i,j]
    index<-c(m1:m2)
    if (!is.na(GS_start[i,j])){
      GSindex<-index[GS_start[i,j]:GS_end[i,j]]
      Tmax_GS[(i-1)*38+j]<-mean(tmax[i,GSindex])
      Tmin_GS[(i-1)*38+j]<-mean(tmin[i,GSindex])
      GDD_GS[(i-1)*38+j]<-sum(GDD[i,GSindex])
      EDD_GS[(i-1)*38+j]<-sum(EDD[i,GSindex])
      VPD_GS[(i-1)*38+j]<-sum(VPD[i,GSindex])
      Pr_GS[(i-1)*38+j]<-sum(pr[i,GSindex])
    }
    m1=m2+1
  }
}
S<-read.table("harvest_area.csv",header=TRUE,sep=",")
area<-rep(NA,38*countynum)
for (i in 1:countynum){
  stateindex<-which(S$State.ANSI==StateANSI[i])
  countyindex<-which(S$County.ANSI==CountyANSI[i])
  yearindex<-intersect(stateindex,countyindex)
  year<-S$Year[yearindex]
  areatoget<-rep(NA,38)
  for (j in 1:38){
    ind<-which(year==1978+j)
    if (!identical(ind,integer(0))){
      areatoget[j]<-S$Value[yearindex[ind]]
    }
  }
  area[((i-1)*38+1):(i*38)]<-areatoget
}
W<-read.table("yielddata.csv",header=TRUE,sep=",")
yield_anomaly<-rep(NA,38*countynum)
yield<-rep(NA,38*countynum)
for (i in 1:countynum){
  stateindex<-which(W$State.ANSI==StateANSI[i])
  countyindex<-which(W$County.ANSI==CountyANSI[i])
  yearindex<-intersect(stateindex,countyindex)
  year<-W$Year[yearindex]
  yieldtoget<-rep(NA,38)
  for (j in 1:38){
    ind<-which(year==1978+j)
    if (!identical(ind,integer(0))){
      yieldtoget[j]<-W$Value[yearindex[ind]]
    }
  }
  if (!all(is.na(yieldtoget))){
    tempyear<-c(1979:2016)
    timetrend<-lm(yield~year+year^2,data=data.frame(yield=yieldtoget[which(!is.na(yieldtoget))]
                                                    ,year=tempyear[which(!is.na(yieldtoget))]))
    predictedmeanyield<-rep(NA,38)
    predictedmeanyield[which(!is.na(yieldtoget))]<-predict(timetrend,data=tempyear)
    yield_anomaly[((i-1)*38+1):(i*38)]<-yieldtoget-predictedmeanyield #anomaly for each state
  } else{
    yield_anomaly[((i-1)*38+1):(i*38)]<-yieldtoget
  }
  yield[((i-1)*38+1):(i*38)]<-yieldtoget
}

Data<-data.frame(StateANSI=rep(StateANSI,each=38),countyANSI=rep(CountyANSI,each=38),fips=rep(ANSI,each=38),
                 year=rep(c(1:38),countynum),Tmax_GS=Tmax_GS,Tmin_GS=Tmin_GS, GDD_GS=GDD_GS,
                 EDD_GS=EDD_GS,VPD_GS=VPD_GS,
                 Pr_GS=Pr_GS,yield=yield,GS_length=GS_length,
                 yield_anomaly=yield_anomaly,area=area)
Data<-Data[complete.cases(Data), ]
save(Data,file = "Metdata_temp/Data_Metobs")


#fifth read projection daily data: 2006-2099

#Determine which county each grid is in 
file<-nc_open("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_GFDL-ESM2M_r1i1p1_historical_1985_1989_CONUS_daily.nc")
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
years<-c(2006:2099)
load("Data/MACAv2-METDATA/Tmax_projmean")
load("Data/MACAv2-METDATA/Tmin_projmean")
load("Data/MACAv2-METDATA/Pr_projmean")
load("Data/MACAv2-METDATA/Rhmax_projmean")
load("Data/MACAv2-METDATA/Rhmin_projmean")
rhmax=Rhmax_mean[ ,c(1:33968)]
rhmin=Rhmin_mean[ ,c(1:33968)]
tmax=Tmax_mean[ ,c(1:33968)]
tmin=Tmin_mean[ ,c(1:33968)]
pr=Pr_mean[ ,c(1:33968)]
rm(Rhmax_mean,Rhmin_mean,Tmax_mean,Tmin_mean,Pr_mean)


RHmean<-(rhmax+rhmin)/2
Tmean<-(tmax+tmin)/2
e_s<-6.112*exp((17.269*Tmean)/(Tmean-237.3)) #saturated vapor pressure (hPa)
VPD<-e_s*(1-RHmean/100)

GDD<-matrix(NA,nrow=countynum,ncol=70*365+23*366)
EDD<-matrix(NA,nrow=countynum,ncol=70*365+23*366)
for (i in 1: countynum){
  for (j in 1: (70*365+23*366)){
    GDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[1]
    EDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[2]
  }
}
for (n in 1:countynum){
  write.table(VPD[n, ],paste("Data/MACAv2-METDATA/dailyVPD/VPD_",countytoget[n],sep = ""), sep=" ")
  write.table(GDD[n, ],paste("Data/MACAv2-METDATA/dailyGDD/GDD_",countytoget[n],sep = ""), sep=" ")
  write.table(EDD[n, ],paste("Data/MACAv2-METDATA/dailyEDD/EDD_",countytoget[n],sep = ""), sep=" ")
}
save(VPD,file="Data/MACAv2-METDATA/VPD_projmean")
save(GDD,file="Data/MACAv2-METDATA/GDD_projmean")
save(EDD,file="Data/MACAv2-METDATA/EDD_projmean")
load("Data/MACAv2-METDATA/VPD_projmean")
load("Data/MACAv2-METDATA/GDD_projmean")
load("Data/MACAv2-METDATA/EDD_projmean")
#GDDrange<-c(600,1200,2050)
Tthres<-10
GS_start<-matrix(NA,nrow=countynum,ncol=93)
GS_end<-matrix(NA,nrow=countynum,ncol=93)
m1=1
m2=1
for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:93){
    k=365
    if (j%%4==3){ #when exist Feb29
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
save(GS_start,file="Metdata_temp/MacaGS_start")
save(GS_end,file="Metdata_temp/MacaGS_end")
load("Metdata_temp/MacaGS_start")
load("Metdata_temp/MacaGS_end")


Tmax_GS<-rep(NA,countynum*93)
Tmin_GS<-rep(NA,countynum*93)
GDD_GS<-rep(NA,countynum*93)
EDD_GS<-rep(NA,countynum*93)
VPD_GS<-rep(NA,countynum*93)
Pr_GS<-rep(NA,countynum*93)
GS_length<-rep(NA,countynum*93)

for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:93){
    k=365
    if (j%%4==3){ #when exist Feb29
      k=366
    }
    m2<-m1+k-1
    GS_length[(i-1)*93+j]<-GS_end[i,j]-GS_start[i,j]
    index<-c(m1:m2)
    if (!is.na(GS_start[i,j])){
      GSindex<-index[GS_start[i,j]:GS_end[i,j]]
      Tmax_GS[(i-1)*93+j]<-mean(tmax[i,GSindex])
      Tmin_GS[(i-1)*93+j]<-mean(tmin[i,GSindex])
      GDD_GS[(i-1)*93+j]<-sum(GDD[i,GSindex])
      EDD_GS[(i-1)*93+j]<-sum(EDD[i,GSindex])
      VPD_GS[(i-1)*93+j]<-sum(VPD[i,GSindex])
      Pr_GS[(i-1)*93+j]<-sum(pr[i,GSindex])
    }
    m1=m2+1
  }
}
Data_proj<-data.frame(StateANSI=rep(StateANSI,each=93),countyANSI=rep(CountyANSI,each=93),fips=rep(ANSI,each=93),
                      year=rep(c(1:93),countynum),Tmax_GS=Tmax_GS,Tmin_GS=Tmin_GS,GDD_GS=GDD_GS,GS_length=GS_length,
                      EDD_GS=EDD_GS,VPD_GS=VPD_GS,Pr_GS=Pr_GS)
save(Data_proj,file = "Metdata_temp/Data_Macaproj")


#sixth regression and plot
load("Metdata_temp/Data_Metobs")
load("Metdata_temp/Data_Macaproj")
Data$StateANSI<-factor(Data$StateANSI)
Data_proj$StateANSI<-factor(Data_proj$StateANSI)
Data$year=Data$year+1978
Data_proj$year=Data_proj$year+2005
Data$year<-factor(Data$year)
Data_proj$year<-factor(Data_proj$year)

Tstructure<-c(5,6,7,8)
Pstructure<-c(9,10)
Tnum<-2^length(Tstructure)
Pnum<-2^length(Pstructure)
totlen=dim(Data)[1]
Tnames<-rep("/",Tnum)
for(i in 1:(Tnum-1)){
  Tindx<-as.binary(i,n=length(Tstructure))
  Tnames[i+1]<-paste(colnames(Data)[Tstructure[Tindx]],collapse = "+")
}
Pnames<-rep("/",Pnum)
for(i in 1:(Pnum-1)){
  Pindx<-as.binary(i,n=length(Pstructure))
  Pnames[i+1]<-paste(colnames(Data)[Pstructure[Pindx]],collapse = "+")
}

hind_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=32)
proj_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=93)
R_sqr<-rep(NA,Tnum*Pnum)
for (i in 1:Tnum){
  for (j in 1:Pnum){
    if (i==1 & j==1){
      model_fix<-lm(yield_anomaly~1,data=Data, weights=Data$area)
    } else if (i==1 & j>1){
      model_fix<-lm(as.formula(paste("yield_anomaly ~ ",Pnames[j], sep="") ),data=Data, weights=Data$area)
    } else if (j==1 & i>1){
      model_fix<-lm(as.formula(paste("yield_anomaly ~ ",Tnames[i], sep="") ),data=Data, weights=Data$area)
    } else {
      model_fix<-lm(as.formula(paste("yield_anomaly ~ ",paste(Tnames[i], Pnames[j],sep="+"), sep="") ),data=Data, weights=Data$area)
    }
    R_sqr[(i-1)*Pnum+j]<-summary(model_fix)$r.squared 
    proj<-predict(model_fix,Data_proj)
    for (k in 1:93){
      indx<-which(Data_proj$year==levels(Data_proj$year)[k])
      proj_fit[(i-1)*Pnum+j,k]<-mean(proj[indx],na.rm=TRUE)
    }
    hind<-predict(model_fix,Data)
    for (k in 1:32){
      indx<-which(Data$year==levels(Data$year)[k])
      hind_fit[(i-1)*Pnum+j,k]<-mean(hind[indx],na.rm=TRUE)
    }
  }
}

save(hind_fit, file="Metdata_temp/hind_fit")
save(proj_fit, file="Metdata_temp/proj_fit")
save(R_sqr,file="Metdata_temp/R_squared")
load("Metdata_temp/hind_fit")
load("Metdata_temp/proj_fit")
meanyield<-rep(NA,32)
meanyield_up<-rep(NA,32)
meanyield_low<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield[i]<-mean(Data$yield_anomaly[indx],na.rm=T)
  meanyield_up[i]<-meanyield[i]+abs(mean(Data$yield[indx],na.rm=T)*0.2)
  meanyield_low[i]<-meanyield[i]-abs(mean(Data$yield[indx],na.rm=T)*0.2)
}

index<-rep(NA,64)
k=1
for (i in 1:64){
  if (all(hind_fit[i, ] > meanyield_low) && all((hind_fit[i, ] < meanyield_up))){
    index[k]<-i
    k=k+1
  }
}
index<-index[!is.na(index)]

#quadratic yield trend in each state 
statenum<-length(levels(Data$StateANSI))
Year<-length(levels(Data$year))
Yield_stateyear<-rep(NA,statenum*Year) # the quadratic yield of each state each year
Yield_stateyear_proj<-rep(NA,statenum*93)
for (i in 1:statenum){
  index<-which(Data$StateANSI == levels(Data$StateANSI)[i])
  Yield_state<-Data$yield[index]
  Year_state<-Data$year[index]
  newdata<-data.frame(Yield=log(Yield_state),year=as.integer(Year_state))
  model<-lm(Yield~poly(year,2),data = newdata) #a quadratic time trend fit for each state
  Yield_predict<-predict(model,newdata) 
  for (j in 1:Year){
    index<-which(Year_state == levels(Data$year)[j])
    Yield_stateyear[(i-1)*Year+j]=mean(Yield_predict[index]) #each county each year has a "trend" value
  }
  
  index<-which(Data_proj$StateANSI == levels(Data_proj$StateANSI)[i])
  Year_state<-Data_proj$year[index]
  projyear<-data.frame(year=as.integer(Year_state))
  Yield_proj<-predict(model,projyear)
  for (j in 1:93){
    index<-which(Year_state == levels(Data_proj$year)[j])
    Yield_stateyear_proj[(i-1)*93+j]=mean(Yield_proj[index])
  }
}
trendindex<-(as.integer(Data$StateANSI)-1)*Year+as.integer(Data$year)  #find the index of detrending yield based on state and year
trendindex_proj<-(as.integer(Data_proj$StateANSI)-1)*93+as.integer(Data_proj$year)
Data$logYield<-log(Data$yield)-Yield_stateyear[trendindex]   #this anomaly is: yield - the quadratic trend

model<-lm(logYield~GDD_GS+EDD_GS+poly(Pr_GS,2),data=Data, weights=Data$area) #year and fips are fixed effects
hind<-predict(model,Data)
SRhind<-exp(hind+Yield_stateyear[trendindex])-exp(Yield_stateyear[trendindex]) #the anomaly hind
proj<-predict(model,Data_proj)
SRproj<-exp(proj+Yield_stateyear_proj[trendindex_proj])-exp(Yield_stateyear_proj[trendindex_proj])
SRyield<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  SRyield[i]<-mean(SRhind[indx],na.rm = TRUE)
}
SRyieldproj<-rep(NA,93)
for (i in 1:93){
  indx<-which(Data_proj$year==levels(Data_proj$year)[i])
  SRyieldproj[i]<-mean(SRproj[indx],na.rm = TRUE)
}
#Sam's model
model<-lm(yield_anomaly ~ GDD_GS+EDD_GS+ poly(Pr_GS,2)  + Tmax_GS + Tmin_GS + GDD_GS:EDD_GS + 
                 GDD_GS:poly(Pr_GS,2) + poly(Pr_GS,2):EDD_GS + GDD_GS:Tmax_GS +
                 GDD_GS:Tmin_GS +  EDD_GS:Tmax_GS + EDD_GS:Tmin_GS + 
                 poly(Pr_GS,2):Tmax_GS + poly(Pr_GS,2):Tmin_GS + Tmax_GS:Tmin_GS ,data= Data)
hind<-predict(model,Data)
proj<-predict(model,Data_proj)
Samhind<-rep(NA,32)
for (i in 1:32){
     indx<-which(Data$year==levels(Data$year)[i])
     Samhind[i]<-mean(hind[indx],na.rm = TRUE)
   }
Samproj<-rep(NA,93)
for (i in 1:93){
     indx<-which(Data_proj$year==levels(Data_proj$year)[i])
     Samproj[i]<-mean(proj[indx],na.rm = TRUE)
   }


par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2098),ylim = c(-45,45),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)
polygon(c(1981:2012,2012:1981),c(meanyield_up,rev(meanyield_low)),col="aliceblue")
for (i in 1:(Tnum*Pnum)){
  lines(c(1981:2012),hind_fit[i, ],col = "forestgreen",type = 'l',lwd=0.3)
  lines(c(2006:2098),proj_fit[i, ],col = "darkorchid",type = 'l',lwd=0.3)
}
points(c(1981:2012),meanyield,col="black",pch=20)
lines(c(1981:2012),SRyield,col="red",lwd=2)
lines(c(2006:2098),SRyieldproj,col="blue",lwd=2)
legend("topright",pch = c(20,NA,NA,NA,NA),lwd=c(NA,2,2,2,2),lty=c(NA,1,1,1,1),
       col=c("black",'forestgreen',"red","darkorchid","blue"),
       legend = c("Yield observation","Hindcasts of selected good models",
                  "Hindcasts of Schlenker and Roberts model","Projections of selected good models",
                  "Projections of Schlenker and Roberts"),bty="n",cex=1.4)


#All
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2098),ylim = c(-45,45),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)
polygon(c(1981:2012,2012:1981),c(meanyield_up,rev(meanyield_low)),col="aliceblue")
for (i in 1:length(index)){
     lines(c(1981:2012),hind_fit[index[i], ],col = "forestgreen",type = 'l',lwd=0.3)
     lines(c(2006:2098),proj_fit[index[i], ],col = "darkorchid",type = 'l',lwd=0.3)
   }
points(c(1981:2012),meanyield,col="black",pch=20)
lines(c(1981:2012),SRyield,col="red",lwd=2)
lines(c(1981:2012),Samhind,col="gold1",lwd=2)
lines(c(2006:2098),SRyieldproj,col="blue",lwd=2)
lines(c(2006:2098),Samproj,col="yellow1",lwd=2)
legend("topright",pch = c(20,NA,NA,NA,NA,NA,NA),lwd=c(NA,2,2,2,2,2,2),lty=c(NA,1,1,1,1,1,1),
                 col=c("black",'forestgreen',"red","gold1","darkorchid","blue","yellow1"),
                 legend = c("Yield observation","Hindcasts of selected good models",
                            "Hindcasts of Schlenker and Roberts model","Hindcasts of Sam's model","Projections of selected good models",
                            "Projections of Schlenker and Roberts","Projections of Sam's model"),bty="n",cex=1.4)

#Hind part
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2012),ylim = c(-25,25),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)
for (i in 1:length(index)){
     lines(c(1981:2012),hind_fit[index[i], ],col = "forestgreen",type = 'l',lwd=0.3)
   }
points(c(1981:2012),meanyield,col="black",pch=20)
lines(c(1981:2012),SRyield,col="red",lwd=2)
lines(c(1981:2012),Samhind,col="gold1",lwd=2)
legend("topleft",pch = c(20,NA,NA,NA),lwd=c(NA,2,2,2),lty=c(NA,1,1,1),col=c("black",'forestgreen',"red","gold1"),
                 legend = c("Yield observation","Hindcasts of selected good models",
                             "Hindcasts of Schlenker and Roberts model","Hindcasts of Sam's model"),bty="n",cex=1.4)

#proj part
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(2006,2098),ylim = c(-40,10),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)
for (i in 1:length(index)){
     lines(c(2006:2098),proj_fit[index[i], ],col = "darkorchid",type = 'l',lwd=0.3)
   }
lines(c(2006:2098),SRyieldproj,col="blue",lwd=2)
lines(c(2006:2098),Samproj,col="yellow1",lwd=2)
legend("topright",pch = c(NA,NA,NA),lwd=c(2,2,2),lty=c(1,1,1),col=c("darkorchid","blue","yellow1"),
                 legend = c("Projections of selected good models",
                            "Projections of Schlenker and Roberts","Projections of Sam's model"),bty="n",cex=1.4)

