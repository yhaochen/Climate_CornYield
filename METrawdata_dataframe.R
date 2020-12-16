#read the raw Metdata dataframe

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
source("latlong2county.R")
source("GDDEDD.R")

#first read daily Tmax, Tmin, RHmax, RHmin, Pr

#Metdata obs: 1979-2016
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
dir.create("/storage/work/h/hxy46/Countywise/Metdata/Metdataframe",recursive = TRUE)
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

save(pr,file="Metdata/Metdataframe/Metpr")
save(tmax,file="Metdata/Metdataframe/Mettmax")
save(tmin,file="Metdata/Metdataframe/Mettmin")
save(rhmax,file="Metdata/Metdataframe/Metrhmax")
save(rhmin,file="Metdata/Metdataframe/Metrhmin")




#second calculate VPD, GDD, EDD
# 
# load("Metdata/Metdataframe/Metpr")
# load("Metdata/Metdataframe/Mettmax")
# load("Metdata/Metdataframe/Mettmin")
# load("Metdata/Metdataframe/Metrhmax")
# load("Metdata/Metdataframe/Metrhmin")

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

save(VPD,file="Metdata/Metdataframe/MetVPD")
save(GDD,file="Metdata/Metdataframe/MetGDD")
save(EDD,file="Metdata/Metdataframe/MetEDD")
load("Metdata/Metdataframe/MetVPD")
load("Metdata/Metdataframe/MetGDD")
load("Metdata/Metdataframe/MetEDD")

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
save(GS_start,file="Metdata/Metdataframe/MetGS_start")
save(GS_end,file="Metdata/Metdataframe/MetGS_end")


#fourth combine to a dataframe
#dataframe should be countynum*yearnum
load("Metdata/Metdataframe/MetGS_start")
load("Metdata/Metdataframe/MetGS_end")
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
save(ANSI,file="ANSI")
save(CountyANSI,file="CountyANSI")
save(StateANSI,file="StateANSI")

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
  yield[((i-1)*38+1):(i*38)]<-yieldtoget
}

Data<-data.frame(StateANSI=rep(StateANSI,each=38),countyANSI=rep(CountyANSI,each=38),fips=rep(ANSI,each=38),
                 year=rep(c(1:38),countynum),Tmax_GS=Tmax_GS,Tmin_GS=Tmin_GS, GDD_GS=GDD_GS,
                 EDD_GS=EDD_GS,VPD_GS=VPD_GS,Pr_GS=Pr_GS,yield=yield,GS_length=GS_length,area=area)
Data<-Data[complete.cases(Data), ] #Yield data are 1981-2012
#calculate yield anomaly based on fixed effects
Data$StateANSI<-factor(Data$StateANSI)
Data$year=Data$year+1978
Data$year<-factor(Data$year)
model<-lm(yield~Tmax_GS+Tmin_GS+GDD_GS+EDD_GS+VPD_GS+Pr_GS+fips+year,data=Data)
Coef<-summary(model)$coefficients
Data$yield_anomaly<-rep(NA,dim(Data)[1])
for (i in 1:dim(Data)[1]){
  row_year<-which(row.names(Coef)==paste("year",Data$year[i],sep=""))
  if (length(row_year)==0){
    yeareffect<-0
  } else{
    yeareffect<-Coef[row_year,1]
  }
  Data$yield_anomaly[i]<-Data$yield[i]-yeareffect
}
Data$yield_anomaly<-Data$yield_anomaly-weighted.mean(Data$yield_anomaly,na.rm=T,w = Data$area)
save(Data,file = "Metdata/Metdataframe/Data_Metobs")


