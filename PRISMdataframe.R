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

# This file reads dataframe for PRISM

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#Read PRISM climate data
datanames<-c("PRISM")
datanum<-length(datanames)

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

#read area data
S<-read.table("harvest_area.csv",header=TRUE,sep=",")
area<-rep(NA,32*countynum)
for (i in 1:countynum){
  stateindex<-which(S$State.ANSI==StateANSI[i])
  countyindex<-which(S$County.ANSI==CountyANSI[i])
  yearindex<-intersect(stateindex,countyindex)
  year<-S$Year[yearindex]
  areatoget<-rep(NA,32)
  for (j in 1:32){
    ind<-which(year==1980+j)
    if (!identical(ind,integer(0))){
      areatoget[j]<-S$Value[yearindex[ind]]
    }
  }
  area[((i-1)*32+1):(i*32)]<-areatoget
}

#read yield data
W<-read.table("yielddata.csv",header=TRUE,sep=",")
yield_anomaly<-rep(NA,32*countynum)
yield<-rep(NA,32*countynum)
for (i in 1:countynum){
  stateindex<-which(W$State.ANSI==StateANSI[i])
  countyindex<-which(W$County.ANSI==CountyANSI[i])
  yearindex<-intersect(stateindex,countyindex)
  year<-W$Year[yearindex]
  yieldtoget<-rep(NA,32)
  for (j in 1:32){
    ind<-which(year==1980+j)
    if (!identical(ind,integer(0))){
      yieldtoget[j]<-W$Value[yearindex[ind]]
    }
  }
  if (!all(is.na(yieldtoget))){
    tempyear<-c(1981:2012)
    timetrend<-lm(yield~year+year^2,data=data.frame(yield=yieldtoget[which(!is.na(yieldtoget))]
                                                    ,year=tempyear[which(!is.na(yieldtoget))]))
    predictedmeanyield<-rep(NA,32)
    predictedmeanyield[which(!is.na(yieldtoget))]<-predict(timetrend,data=tempyear)
    yield_anomaly[((i-1)*32+1):(i*32)]<-yieldtoget-predictedmeanyield #anomaly for each state
  } else{
    yield_anomaly[((i-1)*32+1):(i*32)]<-yieldtoget
  }
  yield[((i-1)*32+1):(i*32)]<-yieldtoget
}



#read GDD EDD data
GDD<-matrix(NA,nrow = 32*countynum,ncol=datanum)
EDD<-matrix(NA,nrow = 32*countynum,ncol=datanum)
GDD_anomaly<-matrix(NA,nrow=32*countynum,ncol = datanum)
EDD_anomaly<-matrix(NA,nrow=32*countynum,ncol = datanum)
for (i in 1:datanum){
  for (j in 1:countynum){
    gdd<-read.table(paste("Data/",datanames[i],"/GDD/GDD_",datanames[i],"_",originalcountynames[j],sep=""))
    edd<-read.table(paste("Data/",datanames[i],"/EDD/EDD_",datanames[i],"_",originalcountynames[j],sep=""))
    Year<-read.table(paste("Data/",datanames[i],"/Year_",datanames[i],sep=""))
    startpoint<-1982-Year[1,1]
    endpoint<-2013-Year[1,1]
    GDD[c(((j-1)*32+1):(j*32)),i]<-gdd[(startpoint:endpoint),1]
    EDD[c(((j-1)*32+1):(j*32)),i]<-edd[(startpoint:endpoint),1]
    GDD_anomaly[c(((j-1)*32+1):(j*32)),i]<-gdd[(startpoint:endpoint),1]-mean(gdd[(startpoint:endpoint),1],na.rm=TRUE)#GDD and EDD also in anomaly space
    EDD_anomaly[c(((j-1)*32+1):(j*32)),i]<-edd[(startpoint:endpoint),1]-mean(edd[(startpoint:endpoint),1],na.rm=TRUE)#careful: anomaly space calculated linealry or quadratically
  }
}


#read precip data and turn to SPI
SPI<-matrix(NA,nrow=32*countynum,ncol=datanum) #annual mean SPI
PR<-matrix(NA,nrow=32*countynum,ncol=datanum)
for (i in 1:datanum){
  for (j in 1:countynum){
    monthly_pr<-read.table(paste("Data/",datanames[i],"/Pr_monthly/Pr_monthly_",datanames[i],"_",originalcountynames[j],sep=""))
    Year<-read.table(paste("Data/",datanames[i],"/Year_",datanames[i],sep=""))
    startpoint<-(1982-Year[1,1]-1)*12+1
    endpoint<-(2013-Year[1,1])*12
    #create the dataframe that can be used to calculate SPI, see in precintcon data(monthly)
    #annual mean precipitation (Mar-Aug) as another variable
    PrFrame<-data.frame(year=rep(c(1981:2012),each=12),month=rep(c(1:12),32),precipitation=monthly_pr[startpoint:endpoint,1])
    spi_monthly<-spi(as.monthly(PrFrame))
    Pr_monthly<-PrFrame[ ,3]
    for (k in 1:32){
      SPI[j*32-32+k,i]<-mean(spi_monthly[(k*12-11+2):(k*12-4),3],na.rm=TRUE)
      PR[j*32-32+k,i]<-sum(Pr_monthly[(k*12-11+2):(k*12-4)],na.rm=TRUE)
    }
  }
}


#save all data in a dataframe, by each county, then each year
Data<-data.frame(Yield=yield,Yield_anomaly=yield_anomaly,GDD=GDD[ ,1],GDD_anomaly=GDD_anomaly[ ,1],EDD=EDD[ ,1],EDD_anomaly=EDD_anomaly[ ,1],Pr=PR[ ,1],SPI=SPI[ ,1],Area=area,
                       StateANSI=rep(StateANSI,each=32),countyANSI=rep(CountyANSI,each=32),fips=rep(ANSI,each=32),year=rep(c(1:32),countynum))


#lat lon of each county centroid
Data_prism$lat<-rep(NA,(32*countynum))
Data_prism$lon<-rep(NA,(32*countynum))
levelind<-as.integer(Data_prism$fips)
for (i in 1:(32*countynum)){
  c<-levels(Data_prism$fips)[levelind[i]] #which level (county)
  index<-which(geoCounty==c)
  Data_prism$lat[i]<-geoCounty[index,5]
  Data_prism$lon[i]<-geoCounty[index,4]
}

#Tmax, Tmin
Data_prism$Tmax<-rep(NA,(32*countynum))
Data_prism$Tmin<-rep(NA,(32*countynum))
for (i in 1:countynum){
  range<-c((i*32-31):(i*32))    
  tmax<-read.table(paste("Data/PRISM/Tmax/Tmax_PRISM_",originalcountynames[i],sep=""))
  Data_prism$Tmax[range]<-tmax[c(1:32),1]
  tmin<-read.table(paste("Data/PRISM/Tmin/Tmin_PRISM_",originalcountynames[i],sep=""))
  Data_prism$Tmin[range]<-tmin[c(1:32),1]
}

#lastSF, firstFF
SF<-matrix(NA,nrow = countynum,ncol = Year) # the last spring frost
FF<-matrix(NA,nrow = countynum,ncol = Year) # the first fall frost
#1. get last/first sf,ff of each state
for (i in 1:countynum){
  for (j in 1:32){
    frostdata<-read.table(paste("Data/PRISM/Frost/Frost_PRISM_",originalcountynames[i],j,sep=""))
    SF[i,j]<-frostdata$V1[max(which(frostdata$V1<183))]-1
    FF[i,j]<-frostdata$V1[min(which(frostdata$V1>183))]-1
  }
}
SF[is.na(SF)]<-61
FF[is.na(FF)]<-334
#save SF, FF
#write.table(SF,"LastSpringFrost")
#write.table(FF,"FirstFallFrost")
GSlength<-FF-SF
#2. get mean date, anomaly, growing season
countymeanSF<-rowMeans(SF)
countymeanFF<-rowMeans(FF)
countymeanGSlength<-countymeanFF-countymeanSF
#SF_anomaly<- SF-matrix(rep(countymeanSF,32),nrow=countynum,ncol=Year)
#FF_anomaly<- FF-matrix(rep(countymeanFF,32),nrow=countynum,ncol=Year)
GDD_changingGS<-matrix(NA,nrow=countynum,ncol=32)
EDD_changingGS<-matrix(NA,nrow=countynum,ncol=32)
Pr_changingGS<-matrix(NA,nrow=countynum,ncol=32)
for (i in 1:countynum){
  GDDdata<-read.table(paste("Data/PRISM/GDD_changingGS/GDD_changingGS_PRISM_",originalcountynames[i],sep = ""))
  EDDdata<-read.table(paste("Data/PRISM/EDD_changingGS/EDD_changingGS_PRISM_",originalcountynames[i],sep = ""))
  Prdata<-read.table(paste("Data/PRISM/Pr_changingGS/Pr_PRISM_",originalcountynames[i],sep = ""))
  Pr_changingGS[i, ]<-Prdata[c(1:32), ]
  GDD_changingGS[i, ]<-GDDdata[c(1:32), ]
  EDD_changingGS[i, ]<-EDDdata[c(1:32), ]
}
Data_prism$GDD_changingGS<-NA
Data_prism$EDD_changingGS<-NA
Data_prism$Pr_changingGS<-NA
Data_prism$lastSF<-NA
Data_prism$firstFF<-NA
for (i in 1:length(Data_prism$GDD)){
  countyindx<-which(levels(Data_prism$fips)==Data_prism$fips[i])
  yearindx<-which(levels(Data_prism$year)==Data_prism$year[i])
  Data_prism$lastSF[i]<-SF[countyindx,yearindx]
  Data_prism$firstFF[i]<-FF[countyindx,yearindx]
  Data_prism$GDD_changingGS[i]<-GDD_changingGS[countyindx,yearindx]
  Data_prism$EDD_changingGS[i]<-EDD_changingGS[countyindx,yearindx]
  Data_prism$Pr_changingGS[i]<-Pr_changingGS[countyindx,yearindx]
}


# July pr
Pr_Jul<-matrix(NA,nrow=countynum,ncol=32)
Julindx<-seq(7, 379, 12)
for (i in 1:countynum){
  Prdata<-read.table(paste("Data/PRISM/Pr_monthly/Pr_monthly_PRISM_",originalcountynames[i],sep = ""))
  Pr_Jul[i, ]<-Prdata[Julindx, ]
}
Data_prism$Pr_Jul<-NA
for (i in 1:length(Data_prism$GDD)){
  countyindx<-which(levels(Data_prism$fips)==Data_prism$fips[i])
  yearindx<-which(levels(Data_prism$year)==Data_prism$year[i])
  Data_prism$Pr_Jul[i]<-Pr_Jul[countyindx,yearindx]
}
save(Data_prism,file="Data_prism")
