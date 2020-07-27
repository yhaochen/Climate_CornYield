#Sample the parametric uncertainty
rm(list = ls())
graphics.off()
library(binaryLogic)
library(ggplot2)
source("parameteruncertainty.R")

colMax <- function(data) {
  apply(data, 2, max)
}
colMin <- function(data) {
  apply(data, 2, min)
}
rowMax <- function(data) {
  apply(data, 1, max)
}
rowMin <- function(data) {
  apply(data, 1, min)
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
parasamplenum<-1000

Tstructure<-c(5,6,7,8)
Pstructure<-c(9,10)
Tnum<-2^length(Tstructure)
Pnum<-2^length(Pstructure)

hind_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=32)
proj_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=94*13)
hind_parasample<-array(NA,c(Tnum*Pnum-1,32,parasamplenum))
proj_parasample<-array(NA,c(Tnum*Pnum-1,94*13,parasamplenum))

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
R_sqr<-rep(NA,Tnum*Pnum)


window of hindcasts
meanyield_anomaly<-rep(NA,32)
windowyield_up<-rep(NA,32)
windowyield_low<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-mean(Data$yield_anomaly[indx],na.rm=T)
}

model<-lm(yield_anomaly~Tmax_GS+Tmin_GS+GDD_GS+EDD_GS+VPD_GS+Pr_GS,data=Data,weights = Data$area)
hind<-predict(model,Data)
parasamplenum=1000
besthind_parasample<-matrix(NA,nrow=32,ncol=parasamplenum)


variablenames<-variable.names(model)
variablenum<-length(variablenames)
  #for hindcast of each parametric sample
  col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
  for (m in 1:(variablenum-1)){
    col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
  }
  covered_percent<-rep(NA,40)
for (i in 1:40){
  sigma<-i*0.05
  MCpara<-parasample(model,parasamplenum,sigma)
  included<-rep(NA,32)
  for (k in 1:32){
    indx<-which(Data$year==levels(Data$year)[k])
    for (n in 1:parasamplenum){
      besthind_parasample[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
    }
    included[k]<-(meanyield_anomaly[k]<max(besthind_parasample[k, ])) && (meanyield_anomaly[k]>min(besthind_parasample[k, ]))
  }
  covered_percent[i]<-sum(included)/32
}


sigma<-which.min(abs(covered_percent-0.95))*0.05
MCpara<-parasample(model,10000,sigma)
for (k in 1:32){
  indx<-which(Data$year==levels(Data$year)[k])
  for (n in 1:parasamplenum){
    besthind_parasample[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
  }
}
windowyield_up<-rowMax(besthind_parasample)
windowyield_low<-rowMin(besthind_parasample)
 save(windowyield_up,file="Metdata_temp/windowyield_up")
 save(windowyield_low,file="Metdata_temp/windowyield_low")
# 
# par(mar=c(4,5.1,1.6,2.1))
# plot(0,0,xlim = c(1981,2099),ylim = c(-50,40),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
# polygon(c(1981:2012,2012:1981),c(windowyield_up,rev(windowyield_low)),col="darkseagreen1",border=NA) #window
# points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=2)

load("Metdata_temp/windowyield_up")
load("Metdata_temp/windowyield_low")

#Uncertainty analysis part
for (i in 1:Tnum){
  for (j in 1:Pnum){
    if (i==1 & j==1){
      model<-lm(yield_anomaly~1,data=Data, weights=Data$area)
    } else if (i==1 & j>1){
      model<-lm(as.formula(paste("yield_anomaly ~ ",Pnames[j], sep="") ),data=Data, weights=Data$area)
    } else if (j==1 & i>1){
      model<-lm(as.formula(paste("yield_anomaly ~ ",Tnames[i], sep="") ),data=Data, weights=Data$area)
    } else {
      model<-lm(as.formula(paste("yield_anomaly ~ ",paste(Tnames[i], Pnames[j],sep="+"), sep="") ),data=Data, weights=Data$area)
    }
    R_sqr[(i-1)*Pnum+j]<-summary(model)$adj.r.squared
    MCpara<-parasample(model,parasamplenum,3)
    variablenames<-variable.names(model)
    variablenum<-length(variablenames)
    if (variablenum!=1) {
      #for hindcast of each parametric sample
      col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
      for (m in 1:(variablenum-1)){
        col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
      }
    }
    hind<-predict(model,Data)
    for (k in 1:32){
      indx<-which(Data$year==levels(Data$year)[k])
      hind_fit[(i-1)*Pnum+j,k]<-mean(hind[indx],na.rm=TRUE)
      if (variablenum!=1) {
        for (n in 1:parasamplenum){
          hind_parasample[(i-1)*Pnum+j-1,k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
        }
      }
    }
    
    
    #projection of each sample
    if (variablenum!=1) {
      col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
    }
    for (q in 1:length(modelnames)){
      load(paste("Metdata_temp/macaprojdataframe/Data_",modelnames[q],sep=""))
      Data_proj$StateANSI<-factor(Data_proj$StateANSI)
      Data_proj$year=Data_proj$year+2005
      Data_proj$year<-factor(Data_proj$year)
      if (variablenum!=1) {
        for (m in 1:(variablenum-1)){
          col_data_proj[m]<-which(colnames(Data_proj)==variablenames[m+1])
        }
      }
      proj<-predict(model,Data_proj)
      for (k in 1:94){
        indx<-which(Data_proj$year==levels(Data_proj$year)[k])
        proj_fit[(i-1)*Pnum+j,(q-1)*94+k]<-mean(proj[indx],na.rm=TRUE)
        if (variablenum!=1) {
          for (n in 1:parasamplenum){
            proj_parasample[(i-1)*Pnum+j-1,(q-1)*94+k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data_proj[indx,col_data_proj]))
          }
        }
      }
    }
    
    #keep the models if hindcast pass the window
    if (variablenum!=1) {
      for (n in 1:parasamplenum){
        for (k in 1:32){
          if ((hind_parasample[(i-1)*Pnum+j-1,k,n] > windowyield_up[k]) || (hind_parasample[(i-1)*Pnum+j-1,k,n] < windowyield_low[k])){
            proj_parasample[(i-1)*Pnum+j-1, ,n] <- NA
            k=32
          }
        }
      }
    }
  }
}

save(hind_fit,file="Metdata_temp/hind_bestfit")
save(proj_fit,file="Metdata_temp/proj_bestfit")
save(hind_parasample,file="Metdata_temp/hind_parasample")
save(proj_parasample,file="Metdata_temp/proj_parasample")



