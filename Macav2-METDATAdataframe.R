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

# This file reads dataframe for Maca-metdata projection 2066-2099 34 yrs

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
GDD<-rep(NA,34*countynum)
EDD<-rep(NA,34*countynum)

for (j in 1:countynum){
  gdd<-read.table(paste("Data/MACAv2-METDATA/GDD/GDD_",originalcountynames[j],sep=""))
  edd<-read.table(paste("Data/MACAv2-METDATA/EDD/EDD_",originalcountynames[j],sep=""))
  Year<-read.table(paste("Data/MACAv2-METDATA/Year_MACAv2-METDATA",sep=""))
  startpoint<-2067-Year[1,1]
  endpoint<-2100-Year[1,1]
  GDD[c(((j-1)*34+1):(j*34))]<-gdd[(startpoint:endpoint),1]
  EDD[c(((j-1)*34+1):(j*34))]<-edd[(startpoint:endpoint),1]
}



#read precip data and turn to SPI
SPI<-rep(NA,34*countynum) #annual mean SPI
PR<-rep(NA,34*countynum)

  for (j in 1:countynum){
    monthly_pr<-read.table(paste("Data/MACAv2-METDATA/Pr_monthly/Pr_monthly",originalcountynames[j],sep=""))
    Year<-read.table(paste("Data/MACAv2-METDATA/Year_MACAv2-METDATA",sep=""))
    startpoint<-(2067-Year[1,1]-1)*12+1
    endpoint<-(2100-Year[1,1])*12
    #create the dataframe that can be used to calculate SPI, see in precintcon data(monthly)
    #annual mean precipitation (Mar-Aug) as another variable
    PrFrame<-data.frame(year=rep(c(2066:2099),each=12),month=rep(c(1:12),34),precipitation=monthly_pr[startpoint:endpoint,1])
    spi_monthly<-spi(as.monthly(PrFrame))
    Pr_monthly<-PrFrame[ ,3]
    for (k in 1:34){
      SPI[j*34-34+k]<-mean(spi_monthly[(k*12-11+2):(k*12-4),3],na.rm=TRUE)
      PR[j*34-34+k]<-sum(Pr_monthly[(k*12-11+2):(k*12-4)],na.rm=TRUE)
    }
  }



#save all data in a dataframe, by each county, then each year
Data_Macamet<-data.frame(GDD=GDD,EDD=EDD,Pr=PR,SPI=SPI,
                       StateANSI=rep(StateANSI,each=34),countyANSI=rep(CountyANSI,each=34),fips=rep(ANSI,each=34),year=rep(c(1:34),countynum))

#lat lon of each county centroid
Data_Macamet$lat<-rep(NA,(34*countynum))
Data_Macamet$lon<-rep(NA,(34*countynum))
levelind<-as.integer(Data_Macamet$fips)
for (i in 1:(34*countynum)){
  c<-levels(Data_Macamet$fips)[levelind[i]] #which level (county)
  index<-which(geoCounty==c)
  Data_Macamet$lat[i]<-geoCounty[index,5]
  Data_Macamet$lon[i]<-geoCounty[index,4]
}

#Tmax, Tmin
Data_Macamet$Tmax<-rep(NA,(34*countynum))
Data_Macamet$Tmin<-rep(NA,(34*countynum))
for (i in 1:countynum){
  range<-c((i*34-33):(i*34))    
  tmax<-read.table(paste("Data/MACAv2-METDATA/Tmax/Tmax_",originalcountynames[i],sep=""))
  Data_Macamet$Tmax[range]<-tmax[c(1:34),1]
  tmin<-read.table(paste("Data/MACAv2-METDATA/Tmin/Tmin_",originalcountynames[i],sep=""))
  Data_Macamet$Tmin[range]<-tmin[c(1:34),1]
}



# July pr
Pr_Jul<-matrix(NA,nrow=countynum,ncol=34)
Julindx<-seq(7, 403, 12)
for (i in 1:countynum){
  Prdata<-read.table(paste("Data/MACAv2-METDATA/Pr_monthly/Pr_monthly",originalcountynames[i],sep = ""))
  Pr_Jul[i, ]<-Prdata[Julindx, ]
}
Data_Macamet$Pr_Jul<-NA
for (i in 1:length(Data_Macamet$GDD)){
  countyindx<-which(levels(Data_Macamet$fips)==Data_Macamet$fips[i])
  yearindx<-Data_Macamet$year[i]
  Data_Macamet$Pr_Jul[i]<-Pr_Jul[countyindx,yearindx]
}
save(Data_Macamet,file="Data_Macamet")
