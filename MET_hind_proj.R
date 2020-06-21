
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
library(adaptMCMC)
library(geoR) #for variogram
library(housingData) #for returning county centroid
library(gstat)
library(fields)
library(binaryLogic)
load("Data_Met")
load("Data_Macamet")
Data_Met<-Data_Met[complete.cases(Data_Met), ]
Data_Met$StateANSI<-factor(Data_Met$StateANSI)
Data_Met$year<-factor(Data_Met$year)
Data_Macamet<-Data_Macamet[complete.cases(Data_Macamet), ]
Data_Macamet$StateANSI<-factor(Data_Macamet$StateANSI)
Data_Macamet$year<-factor(Data_Macamet$year)

missingcounties=c("01037", "05015", "05027", "05049", "05051", "05053", "05065", "05101", "05109", "05113", "05129", 
                  "05133", "05135", "05137", "05139", "13053", "13067", "21095", "21131", "21133", "21159", "21193", 
                  "21195", "22071", "22075", "22087", "22089", "22095", "22109", "24510", "26083", "27031", "27075", 
                  "29510", "34003", "34013", "34017", "34031", "34039", "36005", "36041", "36047", "36059", "36061", 
                  "36079", "36081", "36085", "36087", "36119", "42101", "51013", "51650", "51700", "51710", "54047", "54101", "55078")
for (i in 1:57){
  index<-which(Data_Macamet$fips==missingcounties[i])
  Data_Macamet<-Data_Macamet[-index, ]
}

statenum<-length(levels(Data_Met$StateANSI))
countynum<-length(levels(Data_Met$fips))
Year<-length(levels(Data_Met$year))
Yield_stateyear<-rep(NA,statenum*Year)
Year_proj<-length(levels(Data_Macamet$year))
Yield_stateyear_proj<-rep(NA,statenum*Year_proj)# the quadratic yield of each state each year
for (i in 1:statenum){
  index<-which(Data_Met$StateANSI == levels(Data_Met$StateANSI)[i])
  Yield_state<-Data_Met$Yield[index]
  Year_state<-Data_Met$year[index]
  newdata<-data.frame(Yield=log(Yield_state),year=as.integer(Year_state)+1978,year_sqr=(as.integer(Year_state)+1978)^2)
  model<-lm(Yield~year+year_sqr,data = newdata)
  Yield_predict<-predict(model,newdata) #then take average yield in each year
  for (j in 1:Year){
    index<-which(Year_state == levels(Data_Met$year)[j])
    Yield_stateyear[(i-1)*Year+j]=mean(Yield_predict[index])
  }
  
  #find the trend for projection years
  index_proj<-which(Data_Macamet$StateANSI == levels(Data_Macamet$StateANSI)[i])
  Yield_state_proj<-Data_Macamet$Yield[index_proj]
  Year_state_proj<-Data_Macamet$year[index_proj]
  newdata<-data.frame(year=as.integer(Year_state_proj)+2065,year_sqr=(as.integer(Year_state_proj)+2065)^2)
  Yield_predict_proj<-predict(model,newdata)
  for (k in 1:Year_proj){
    index<-which(Year_state_proj == levels(Data_Macamet$year)[k])
    Yield_stateyear_proj[(i-1)*Year_proj+k]=mean(Yield_predict_proj[index])
  }
}
trendindex<-(as.integer(Data_Met$StateANSI)-1)*Year+as.integer(Data_Met$year)  #find the index of detrending yield based on state and year
trendindex_proj<-(as.integer(Data_Macamet$StateANSI)-1)*Year_proj+as.integer(Data_Macamet$year)
Data_Met$logYield<-log(Data_Met$Yield)-Yield_stateyear[trendindex]   #this anomaly is: yield - the quadratic trend



# a list of model structures and CV errors
Tstructure<-c(1,2)
Pstructure<-c(3,4)
Tnum<-2^length(Tstructure)
Pnum<-2^length(Pstructure)
totlen=dim(Data_Met)[1]
AIC<-matrix(NA,nrow = Tnum,ncol=Pnum)
BIC<-matrix(NA,nrow = Tnum,ncol=Pnum)

Tnames<-rep("/",Tnum)
for(i in 1:(Tnum-1)){
  Tindx<-as.binary(i,n=length(Tstructure))
  Tnames[i+1]<-paste(colnames(Data_Met)[Tstructure[Tindx]],collapse = "+")
}
Pnames<-rep("/",Pnum)
for(i in 1:(Pnum-1)){
  Pindx<-as.binary(i,n=length(Pstructure))
  Pnames[i+1]<-paste(colnames(Data_Met)[Pstructure[Pindx]],collapse = "+")
}
hind_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=Year)
hind_low<-matrix(NA,nrow = Tnum*Pnum,ncol=Year)
hind_high<-matrix(NA,nrow = Tnum*Pnum,ncol=Year)
pred_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=Year_proj)
pred_low<-matrix(NA,nrow = Tnum*Pnum,ncol=Year_proj)
pred_high<-matrix(NA,nrow = Tnum*Pnum,ncol=Year_proj)
for (i in 1:Tnum){
  for (j in 1:Pnum){
    if (i==1 & j==1){
      model_fix<-lm(logYield~fips,data=Data_Met, weights=Data_Met$Area)
    } else if (i==1 & j>1){
      model_fix<-lm(as.formula(paste("logYield ~ fips+",Pnames[j], sep="") ),data=Data_Met, weights=Data_Met$Area)
    } else if (j==1 & i>1){
      model_fix<-lm(as.formula(paste("logYield ~ fips+",Tnames[i], sep="") ),data=Data_Met, weights=Data_Met$Area)
    } else {
      model_fix<-lm(as.formula(paste("logYield ~ fips+",paste(Tnames[i], Pnames[j],sep="+"), sep="") ),data=Data_Met, weights=Data_Met$Area)
    }
     pred<-predict(model_fix,Data_Macamet,interval = 'confidence')
     for (k in 1:Year_proj){
       indx<-which(Data_Macamet$year==levels(Data_Macamet$year)[k])
       pred_fit[(i-1)*Tnum+j,k]<-mean(exp(pred[indx,1]+Yield_stateyear_proj[trendindex_proj[indx]]))
       pred_low[(i-1)*Tnum+j,k]<-mean(exp(pred[indx,2]+Yield_stateyear_proj[trendindex_proj[indx]]))
       pred_high[(i-1)*Tnum+j,k]<-mean(exp(pred[indx,3]+Yield_stateyear_proj[trendindex_proj[indx]]))
     }
    
    hind<-predict(model_fix,Data_Met,interval = 'confidence')
    for (k in 1:Year){
      indx<-which(Data_Met$year==levels(Data_Met$year)[k])
      hind_fit[(i-1)*Tnum+j,k]<-mean(exp(hind[indx,1]+Yield_stateyear[trendindex[indx]]))
      hind_low[(i-1)*Tnum+j,k]<-mean(exp(hind[indx,2]+Yield_stateyear[trendindex[indx]]))
      hind_high[(i-1)*Tnum+j,k]<-mean(exp(hind[indx,3]+Yield_stateyear[trendindex[indx]]))
    }
  }
}

names<-rep(NA,16)
for (i in 1:4){
  for (j in 1:4){
    if (i==1 & j==1){
      names[(i-1)*4+j]="/"
    } else if (i==1 & j>1){
      names[(i-1)*4+j]=Pnames[j]
    } else if (j==1 & i>1){
      names[(i-1)*4+j]=Tnames[i]
    } else {
      names[(i-1)*4+j]=paste(Tnames[i], Pnames[j],sep="+")
    }
  }
}
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(2066,2099),ylim = c(580,7800),xlab="Year",ylab="Yield (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)
cl <- rainbow(16)
for (i in 1:(Tnum*Pnum)){
  lines(c(2066:2099),pred_fit[i, ],col = cl[i],type = 'l',lwd=2)
  lines(c(2066:2099),pred_high[i, ],col = cl[i],type = 'l',lty=2,lwd=1)
  lines(c(2066:2099),pred_low[i, ],col = cl[i],type = 'l',lty=2,lwd=1)
}
#legend("topleft",col=cl,lty = rep(1,16),legend = names)

meanyield<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data_Met$year==levels(Data_Met$year)[i])
  meanyield[i]<-mean(Data_Met$Yield[indx],na.rm=T)
}
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2012),ylim = c(60,150),xlab="Year",ylab="Yield (bush/acre)",type = "n",cex.axis=1.4,cex.lab=1.4)
cl <- rainbow(16)
for (i in 1:(Tnum*Pnum)){
  lines(c(1981:2012),hind_fit[i, ],col = cl[i],type = 'l',lwd=2)
  lines(c(1981:2012),hind_high[i, ],col = cl[i],type = 'l',lty=2,lwd=1)
  lines(c(1981:2012),hind_low[i, ],col = cl[i],type = 'l',lty=2,lwd=1)
}
points(c(1981:2012),meanyield,col="black",pch=20)
legend("topleft",pch = 20,legend = "yield observation",cex=2)


confint(model_fix, 'GDD', level=0.95)

