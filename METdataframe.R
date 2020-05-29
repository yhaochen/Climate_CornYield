#remove all data/variables and plots
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
source("latlong2county.R")
source("GDDEDD.R")

# This file reads dataframe for Metdata observation 1979-2016 38 yrs

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
#find the counties based on given states
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

#read GDD EDD data
GDD<-rep(NA,38*countynum)
EDD<-rep(NA,38*countynum)

for (j in 1:countynum){
  gdd<-read.table(paste("Data/MET/GDD/GDD_MET_",originalcountynames[j],sep=""))
  edd<-read.table(paste("Data/MET/EDD/EDD_MET_",originalcountynames[j],sep=""))
  Year<-read.table(paste("Data/MET/Year_MET",sep=""))
  startpoint<-1980-Year[1,1]
  endpoint<-2017-Year[1,1]
  GDD[c(((j-1)*38+1):(j*38))]<-gdd[(startpoint:endpoint),1]
  EDD[c(((j-1)*38+1):(j*38))]<-edd[(startpoint:endpoint),1]
}



#read precip data and turn to SPI
SPI<-rep(NA,38*countynum) #annual mean SPI
PR<-rep(NA,38*countynum)

for (j in 1:countynum){
  monthly_pr<-read.table(paste("Data/MET/Pr_monthly/Pr_MET_",originalcountynames[j],sep=""))
  Year<-read.table(paste("Data/MET/Year_MET",sep=""))
  startpoint<-(1980-Year[1,1]-1)*12+1
  endpoint<-(2017-Year[1,1])*12
  #create the dataframe that can be used to calculate SPI, see in precintcon data(monthly)
  #annual mean precipitation (Mar-Aug) as another variable
  PrFrame<-data.frame(year=rep(c(1979:2016),each=12),month=rep(c(1:12),38),precipitation=monthly_pr[startpoint:endpoint,1])
  spi_monthly<-spi(as.monthly(PrFrame))
  Pr_monthly<-PrFrame[ ,3]
  for (k in 1:38){
    SPI[j*38-38+k]<-mean(spi_monthly[(k*12-11+2):(k*12-4),3],na.rm=TRUE)
    PR[j*38-38+k]<-sum(Pr_monthly[(k*12-11+2):(k*12-4)],na.rm=TRUE)
  }
}

#read area data
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

#read yield data
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

#save all data in a dataframe, by each county, then each year
Data_Met<-data.frame(GDD=GDD,EDD=EDD,Pr=PR,SPI=SPI,Area=area,Yield=yield,Yield_anomaly=yield_anomaly,
                         StateANSI=rep(StateANSI,each=38),countyANSI=rep(CountyANSI,each=38),fips=rep(ANSI,each=38),year=rep(c(1:38),countynum))

#lat lon of each county centroid
Data_Met$lat<-rep(NA,(38*countynum))
Data_Met$lon<-rep(NA,(38*countynum))
levelind<-as.integer(Data_Met$fips)
for (i in 1:(38*countynum)){
  c<-levels(Data_Met$fips)[levelind[i]] #which level (county)
  index<-which(geoCounty==c)
  Data_Met$lat[i]<-geoCounty[index,5]
  Data_Met$lon[i]<-geoCounty[index,4]
}

#Tmax, Tmin
Data_Met$Tmax<-rep(NA,(38*countynum))
Data_Met$Tmin<-rep(NA,(38*countynum))
for (i in 1:countynum){
  range<-c((i*38-37):(i*38))    
  tmax<-read.table(paste("Data/MET/Tmax/Tmax_MET_",originalcountynames[i],sep=""))
  Data_Met$Tmax[range]<-tmax[c(1:38),1]
  tmin<-read.table(paste("Data/MET/Tmin/Tmin_MET_",originalcountynames[i],sep=""))
  Data_Met$Tmin[range]<-tmin[c(1:38),1]
}



# July pr
Pr_Jul<-matrix(NA,nrow=countynum,ncol=38)
Julindx<-seq(7, 451, 12)
for (i in 1:countynum){
  Prdata<-read.table(paste("Data/MET/Pr_monthly/Pr_MET_",originalcountynames[i],sep = ""))
  Pr_Jul[i, ]<-Prdata[Julindx, ]
}
Data_Met$Pr_Jul<-NA
for (i in 1:length(Data_Met$GDD)){
  countyindx<-which(levels(Data_Met$fips)==Data_Met$fips[i])
  yearindx<-Data_Met$year[i]
  Data_Met$Pr_Jul[i]<-Pr_Jul[countyindx,yearindx]
}
save(Data_Met,file="Data_Met")
