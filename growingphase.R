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

load("Metdata_temp/Metpr_temp")
load("Metdata_temp/Mettmax_temp")
load("Metdata_temp/Mettmin_temp")
load("Metdata_temp/Metrhmax_temp")
load("Metdata_temp/Metrhmin_temp")

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
save(VPD,file="Metdata_temp/MetVPD_temp")
save(GDD,file="Metdata_temp/MetGDD_temp")
save(EDD,file="Metdata_temp/MetEDD_temp")

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

GDDrange<-c(600,1200,2050)
Tthres<-10
GS_start<-matrix(NA,nrow=countynum,ncol=38)
phase1<-matrix(NA,nrow=countynum,ncol=38)
phase2<-matrix(NA,nrow=countynum,ncol=38)
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
    GDDtoget<-GDD[i,c(m1:m2)]
    GDDMA<-ma(GDDtoget,21)
    GS_start[i,j]<-which(GDDMA >= Tthres)[1]
    if (!is.na(GS_start[i,j])){
      GDD_cumu<-su(GDDtoget,GS_start[i,j])
      phase1[i,j]<-which(GDD_cumu >= GDDrange[1])[1]
      phase2[i,j]<-which(GDD_cumu >= GDDrange[2])[1]
      GS_end[i,j]<-which(GDD_cumu >= GDDrange[3])[1]
      if (is.na(GS_end[i,j])){
        GS_end[i,j]<-k
      }
      if (is.na(phase2[i,j])){
        GS_start[i,j]<-NA
        phase1[i,j]<-NA
        phase2[i,j]<-NA
        GS_end[i,j]<-NA
      }
    }
    m1=m2+1
  }
}
save(GS_start,file="Metdata_temp/MetGS_start")
save(phase1,file="Metdata_temp/Metphase1")
save(phase2,file="Metdata_temp/Metphase2")
save(GS_end,file="Metdata_temp/MetGS_end")


#fourth combine to a dataframe and regression
#dataframe should be countynum*yearnum
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

Tmax_1<-rep(NA,countynum*38)
Tmax_2<-rep(NA,countynum*38)
Tmax_3<-rep(NA,countynum*38)
Tmax_GS<-rep(NA,countynum*38)
Tmin_1<-rep(NA,countynum*38)
Tmin_2<-rep(NA,countynum*38)
Tmin_3<-rep(NA,countynum*38)
Tmin_GS<-rep(NA,countynum*38)
EDD_1<-rep(NA,countynum*38)
EDD_2<-rep(NA,countynum*38)
EDD_3<-rep(NA,countynum*38)
EDD_GS<-rep(NA,countynum*38)
VPD_1<-rep(NA,countynum*38)
VPD_2<-rep(NA,countynum*38)
VPD_3<-rep(NA,countynum*38)
VPD_GS<-rep(NA,countynum*38)
Pr_1<-rep(NA,countynum*38)
Pr_2<-rep(NA,countynum*38)
Pr_3<-rep(NA,countynum*38)
Pr_GS<-rep(NA,countynum*38)

for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:38){
    k=365
    if (j%%4==2){ #when exist Feb29
      k=366
    }
    m2<-m1+k-1
    
    index<-c(m1:m2)
    if (!is.na(GS_start[i,j])){
      GSindex<-index[GS_start[i,j]:GS_end[i,j]]
      phase1index<-index[GS_start[i,j]:phase1[i,j]]
      phase2index<-index[(phase1[i,j]+1):phase2[i,j]]
      phase3index<-index[(phase2[i,j]+1):GS_end[i,j]]
      
      Tmax_GS[(i-1)*38+j]<-mean(tmax[i,GSindex])
      Tmax_1[(i-1)*38+j]<-mean(tmax[i,phase1index])
      Tmax_2[(i-1)*38+j]<-mean(tmax[i,phase2index])
      Tmax_3[(i-1)*38+j]<-mean(tmax[i,phase3index])
      
      Tmin_GS[(i-1)*38+j]<-mean(tmin[i,GSindex])
      Tmin_1[(i-1)*38+j]<-mean(tmin[i,phase1index])
      Tmin_2[(i-1)*38+j]<-mean(tmin[i,phase2index])
      Tmin_3[(i-1)*38+j]<-mean(tmin[i,phase3index])
      
      EDD_GS[(i-1)*38+j]<-sum(EDD[i,GSindex])
      EDD_1[(i-1)*38+j]<-sum(EDD[i,phase1index])
      EDD_2[(i-1)*38+j]<-sum(EDD[i,phase2index])
      EDD_3[(i-1)*38+j]<-sum(EDD[i,phase3index])
      
      VPD_GS[(i-1)*38+j]<-sum(VPD[i,GSindex])
      VPD_1[(i-1)*38+j]<-sum(VPD[i,phase1index])
      VPD_2[(i-1)*38+j]<-sum(VPD[i,phase2index])
      VPD_3[(i-1)*38+j]<-sum(VPD[i,phase3index])
      
      Pr_GS[(i-1)*38+j]<-sum(pr[i,GSindex])
      Pr_1[(i-1)*38+j]<-sum(pr[i,phase1index])
      Pr_2[(i-1)*38+j]<-sum(pr[i,phase2index])
      Pr_3[(i-1)*38+j]<-sum(pr[i,phase3index])
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
                 year=rep(c(1:38),countynum),Tmax_GS=Tmax_GS,Tmax_1=Tmax_1,Tmax_2=Tmax_2,Tmax_3=Tmax_3,Tmin_GS=Tmin_GS,
                 Tmin_1=Tmin_1,Tmin_2=Tmin_2,Tmin_3=Tmin_3,EDD_GS=EDD_GS,EDD_1=EDD_1,EDD_2=EDD_2,EDD_3=EDD_3,VPD_GS=VPD_GS,
                 VPD_1=VPD_1,VPD_2=VPD_2,VPD_3=VPD_3,Pr_GS=Pr_GS,Pr_1=Pr_1,Pr_2=Pr_2,Pr_3=Pr_3,yield=yield,
                 yield_anomaly=yield_anomaly,area=area)
Data<-Data[complete.cases(Data), ]
save(Data,file = "Metdata_temp/Data_Metobs")


#fifth read projection daily data: 2066-2099

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
years<-c(2066:2099)
tmax<-matrix(NA,nrow=countynum,ncol=26*365+8*366)
tmin<-matrix(NA,nrow=countynum,ncol=26*365+8*366)
pr<-matrix(NA,nrow=countynum,ncol=26*365+8*366)
rhmax<-matrix(NA,nrow=countynum,ncol=26*365+8*366)
rhmin<-matrix(NA,nrow=countynum,ncol=26*365+8*366)

m1<-1 #used to record which columns to write in each loop
m2<-1
n1<-1

for (i in 1:34){
  if ((i%%5==1)&(years[i]<2096)){
    Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_GFDL-ESM2M_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmin_GFDL-ESM2M_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_pr_GFDL-ESM2M_r1i1p1_rcp85_"
                   ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    RHmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmax_GFDL-ESM2M_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    RHmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmin_GFDL-ESM2M_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+4,"_CONUS_daily.nc",sep="")
    
    mettmax<-nc_open(Tmax_file)
    mettmin<-nc_open(Tmin_file)
    metpr<-nc_open(Pr_file)
    metrhmax<-nc_open(RHmax_file)
    metrhmin<-nc_open(RHmin_file)
    n1=1
  }
  
  if ((i%%5==1)&(years[i]>2095)){
    Tmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmax_GFDL-ESM2M_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    Tmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_tasmin_GFDL-ESM2M_r1i1p1_rcp85_"
                     ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    Pr_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_pr_GFDL-ESM2M_r1i1p1_rcp85_"
                   ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    RHmax_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmax_GFDL-ESM2M_r1i1p1_rcp85_"
                      ,years[i],"_",years[i]+3,"_CONUS_daily.nc",sep="")
    RHmin_file<-paste("/gpfs/group/kzk10/default/private/data_archive/MACAv2-METDATA/raw/macav2metdata_rhsmin_GFDL-ESM2M_r1i1p1_rcp85_"
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
  Tmax<-ncvar_get(mettmax,varid = "air_temperature",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  Tmin<-ncvar_get(mettmin,varid = "air_temperature",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  Pr<-ncvar_get(metpr,varid = "precipitation",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  RHmax<-ncvar_get(metrhmax,varid = "relative_humidity",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  RHmin<-ncvar_get(metrhmin,varid = "relative_humidity",start = c(596,1,n1),count = c(dim_lon-595,dim_lat,k))
  n1<-n1+k
  
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
    
    print(n)
  }
  m1<-m2+1
  print(i)
}

load("Metdata_temp/Macametpr_temp")
load("Metdata_temp/Macamettmax_temp")
load("Metdata_temp/Macamettmin_temp")
load("Metdata_temp/Macametrhmax_temp")
load("Metdata_temp/Macametrhmin_temp")


for (n in 1:countynum){
  write.table(tmax[n, ],paste("Data/MACAv2-METDATA/dailyTmax/Tmax_",countytoget[n],sep = ""), sep=" ")
  write.table(tmin[n, ],paste("Data/MACAv2-METDATA/dailyTmin/Tmin_",countytoget[n],sep = ""), sep=" ")
  #write.table(pr[n, ],paste("Data/MACAv2-METDATA/dailyPr/Pr_",countytoget[n],sep = ""), sep=" ")
  #write.table(rhmax[n, ],paste("Data/MACAv2-METDATA/dailyRHmax/RHmax_MET_",countytoget[n],sep = ""), sep=" ")
  #write.table(rhmin[n, ],paste("Data/MACAv2-METDATA/dailyRHmin/RHmin_MET_",countytoget[n],sep = ""), sep=" ")
}


RHmean<-(rhmax+rhmin)/2
Tmean<-(tmax+tmin)/2
e_s<-6.112*exp((17.269*Tmean)/(Tmean-237.3)) #saturated vapor pressure (hPa)
VPD<-e_s*(1-RHmean/100)

GDD<-matrix(NA,nrow=countynum,ncol=26*365+8*366)
EDD<-matrix(NA,nrow=countynum,ncol=26*365+8*366)
for (i in 1: countynum){
  for (j in 1: (26*365+8*366)){
    GDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[1]
    EDD[i,j]<-GDDEDD(tmax[i,j],tmin[i,j])[2]
  }
}
for (n in 1:countynum){
  write.table(VPD[n, ],paste("Data/MACAv2-METDATA/dailyVPD/VPD_",countytoget[n],sep = ""), sep=" ")
  write.table(GDD[n, ],paste("Data/MACAv2-METDATA/dailyGDD/GDD_",countytoget[n],sep = ""), sep=" ")
  write.table(EDD[n, ],paste("Data/MACAv2-METDATA/dailyEDD/EDD_",countytoget[n],sep = ""), sep=" ")
}
save(VPD,file="Metdata_temp/MacaVPD_temp")
save(GDD,file="Metdata_temp/MacaGDD_temp")
save(EDD,file="Metdata_temp/MacaEDD_temp")


GDDrange<-c(600,1200,2050)
Tthres<-10
GS_start<-matrix(NA,nrow=countynum,ncol=34)
phase1<-matrix(NA,nrow=countynum,ncol=34)
phase2<-matrix(NA,nrow=countynum,ncol=34)
GS_end<-matrix(NA,nrow=countynum,ncol=34)
m1=1
m2=1
for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:34){
    k=365
    if (j%%4==3){ #when exist Feb29
      k=366
    }
    m2<-m1+k-1
    GDDtoget<-GDD[i,c(m1:m2)]
    GDDMA<-ma(GDDtoget,21)
    GS_start[i,j]<-which(GDDMA >= Tthres)[1]
    if (!is.na(GS_start[i,j])){
      GDD_cumu<-su(GDDtoget,GS_start[i,j])
      phase1[i,j]<-which(GDD_cumu >= GDDrange[1])[1]
      phase2[i,j]<-which(GDD_cumu >= GDDrange[2])[1]
      GS_end[i,j]<-which(GDD_cumu >= GDDrange[3])[1]
      if (is.na(GS_end[i,j])){
        GS_end[i,j]<-k
      }
      if (is.na(phase2[i,j])){
        GS_start[i,j]<-NA
        phase1[i,j]<-NA
        phase2[i,j]<-NA
        GS_end[i,j]<-NA
      }
    }
    m1=m2+1
  }
}
save(GS_start,file="Metdata_temp/MacaGS_start")
save(phase1,file="Metdata_temp/Macaphase1")
save(phase2,file="Metdata_temp/Macaphase2")
save(GS_end,file="Metdata_temp/MacaGS_end")


Tmax_1<-rep(NA,countynum*34)
Tmax_2<-rep(NA,countynum*34)
Tmax_3<-rep(NA,countynum*34)
Tmax_GS<-rep(NA,countynum*34)
Tmin_1<-rep(NA,countynum*34)
Tmin_2<-rep(NA,countynum*34)
Tmin_3<-rep(NA,countynum*34)
Tmin_GS<-rep(NA,countynum*34)
EDD_1<-rep(NA,countynum*34)
EDD_2<-rep(NA,countynum*34)
EDD_3<-rep(NA,countynum*34)
EDD_GS<-rep(NA,countynum*34)
VPD_1<-rep(NA,countynum*34)
VPD_2<-rep(NA,countynum*34)
VPD_3<-rep(NA,countynum*34)
VPD_GS<-rep(NA,countynum*34)
Pr_1<-rep(NA,countynum*34)
Pr_2<-rep(NA,countynum*34)
Pr_3<-rep(NA,countynum*34)
Pr_GS<-rep(NA,countynum*34)

for (i in 1:countynum){
  m1=1
  m2=1
  for (j in 1:34){
    k=365
    if (j%%4==3){ #when exist Feb29
      k=366
    }
    m2<-m1+k-1
    
    index<-c(m1:m2)
    if (!is.na(GS_start[i,j])){
      GSindex<-index[GS_start[i,j]:GS_end[i,j]]
      phase1index<-index[GS_start[i,j]:phase1[i,j]]
      phase2index<-index[(phase1[i,j]+1):phase2[i,j]]
      phase3index<-index[(phase2[i,j]+1):GS_end[i,j]]
      
      Tmax_GS[(i-1)*34+j]<-mean(tmax[i,GSindex])
      Tmax_1[(i-1)*34+j]<-mean(tmax[i,phase1index])
      Tmax_2[(i-1)*34+j]<-mean(tmax[i,phase2index])
      Tmax_3[(i-1)*34+j]<-mean(tmax[i,phase3index])
      
      Tmin_GS[(i-1)*34+j]<-mean(tmin[i,GSindex])
      Tmin_1[(i-1)*34+j]<-mean(tmin[i,phase1index])
      Tmin_2[(i-1)*34+j]<-mean(tmin[i,phase2index])
      Tmin_3[(i-1)*34+j]<-mean(tmin[i,phase3index])
      
      EDD_GS[(i-1)*34+j]<-sum(EDD[i,GSindex])
      EDD_1[(i-1)*34+j]<-sum(EDD[i,phase1index])
      EDD_2[(i-1)*34+j]<-sum(EDD[i,phase2index])
      EDD_3[(i-1)*34+j]<-sum(EDD[i,phase3index])
      
      VPD_GS[(i-1)*34+j]<-sum(VPD[i,GSindex])
      VPD_1[(i-1)*34+j]<-sum(VPD[i,phase1index])
      VPD_2[(i-1)*34+j]<-sum(VPD[i,phase2index])
      VPD_3[(i-1)*34+j]<-sum(VPD[i,phase3index])
      
      Pr_GS[(i-1)*34+j]<-sum(pr[i,GSindex])
      Pr_1[(i-1)*34+j]<-sum(pr[i,phase1index])
      Pr_2[(i-1)*34+j]<-sum(pr[i,phase2index])
      Pr_3[(i-1)*34+j]<-sum(pr[i,phase3index])
    }
    m1=m2+1
  }
}

Data_proj<-data.frame(StateANSI=rep(StateANSI,each=34),countyANSI=rep(CountyANSI,each=34),fips=rep(ANSI,each=34),
                      year=rep(c(1:34),countynum),Tmax_GS=Tmax_GS,Tmax_1=Tmax_1,Tmax_2=Tmax_2,Tmax_3=Tmax_3,Tmin_GS=Tmin_GS,
                      Tmin_1=Tmin_1,Tmin_2=Tmin_2,Tmin_3=Tmin_3,EDD_GS=EDD_GS,EDD_1=EDD_1,EDD_2=EDD_2,EDD_3=EDD_3,VPD_GS=VPD_GS,
                      VPD_1=VPD_1,VPD_2=VPD_2,VPD_3=VPD_3,Pr_GS=Pr_GS,Pr_1=Pr_1,Pr_2=Pr_2,Pr_3=Pr_3)
save(Data_proj,file = "Metdata_temp/Data_Macaproj")


#sixth regression and plot
load("Metdata_temp/Data_Metobs")
load("Metdata_temp/Data_Macaproj")
Data$StateANSI<-factor(Data$StateANSI)
Data$year<-factor(Data$year)
Data_proj$StateANSI<-factor(Data_proj$StateANSI)
Data_proj$year<-factor(Data_proj$year)
Tstructure<-c(6,7,8,10,11,12,14,15,16)
Pstructure<-c(18,19,20,22,23,24)
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
proj_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=34)
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
    
    for (k in 1:34){
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

meanyield<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield[i]<-mean(Data$yield_anomaly[indx],na.rm=T)
}

index<-which(R_sqr>=0.32)
maxindex<-which(R_sqr==max(R_sqr))
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2012),ylim = c(-45,45),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)

for (i in 1:length(index)){
  lines(c(1981:2012),hind_fit[index[i], ],col = "darkseagreen1",type = 'l',lwd=0.2)
}
lines(c(1981:2012),hind_fit[maxindex, ],col = "firebrick",type = 'l',lwd=2)
points(c(1981:2012),meanyield,col="black",pch=20)
legend("topleft",pch = c(20,NA,NA),lwd=c(NA,2,2),col=c("black",'darkseagreen1',"firebrick"),
       legend = c(expression(paste("hindcasts of models with R"^"2", ">= 0.32")),expression(paste("hindcast of the model with maximum R"^2,"(0.343)"))),bty="n",cex=1.4)



par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2012),ylim = c(-45,45),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)

for (i in 1:(Pnum*Tnum)){
  lines(c(1981:2012),hind_fit[i, ],col = "darkseagreen1",type = 'l',lwd=0.2)
}
lines(c(1981:2012),hind_fit[maxindex, ],col = "firebrick",type = 'l',lwd=2)
points(c(1981:2012),meanyield,col="black",pch=20)
legend("topleft",pch = c(20,NA,NA),lwd=c(NA,2,2),col=c("black",'darkseagreen1',"firebrick"),
       legend = c("yield observation","hindcasts of all proposed model structures",expression(paste("hindcast of the model with maximum R"^2,"(0.343)"))),bty="n",cex=1.4)



par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(2066,2099),ylim = c(-45,45),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)

for (i in 1:(Pnum*Tnum)){
  lines(c(2066:2099),proj_fit[i, ],col = "darkseagreen1",type = 'l',lwd=0.2)
}
lines(c(2066:2099),proj_fit[maxindex, ],col = "firebrick",type = 'l',lwd=2)
legend("topleft",lwd=c(2,2),col=c('darkseagreen1',"firebrick"),
       legend = c("projections of all proposed model structures",expression(paste("projection of the model with maximum R"^2,"(0.343)"))),bty="n",cex=1.4)




par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(2066,2099),ylim = c(-45,45),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)

for (i in 1:length(index)){
  lines(c(2066:2099),proj_fit[index[i], ],col = "darkseagreen1",type = 'l',lwd=0.2)
}
lines(c(2066:2099),proj_fit[maxindex, ],col = "firebrick",type = 'l',lwd=2)
legend("topleft",lwd=c(2,2),col=c('darkseagreen1',"firebrick"),
       legend = c(expression(paste("projections of models with R"^"2", ">= 0.32")),
                  expression(paste("projection of the model with maximum R"^2,"(0.343)"))),bty="n",cex=1.4)

