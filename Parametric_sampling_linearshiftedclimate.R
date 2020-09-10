#Sample the parametric uncertainty for two linear shifted climate scenario
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



#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
Data$year=Data$year+1978
Data$year<-factor(Data$year)
meanyield_anomaly<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,w = Data$area[indx])
}
#with each climate projection, save the hind/proj yield data of 64 structures
parasamplenum<-1000
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
#for the additional linear shifted climate proj (2020-2049; 2070-2099)
proj_linearshifted_fit<-matrix(NA,nrow = Tnum*Pnum,ncol=30)
proj_linearshifted_parasample<-array(NA,c(Tnum*Pnum,30,parasamplenum))
load("Metdata/window_sigma")

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
   
   variablenames<-variable.names(model)
   variablenum<-length(variablenames)

   #projection of each sample
   if (variablenum!=1) {
     MCpara<-parasample(model,parasamplenum,window_sigma[(i-1)*Pnum+j])
     col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
     load("Metdata/macaprojdataframe/Data_linearshifted_2070_2099")
     Data_proj_linearshifted_2070_2099$StateANSI<-factor(Data_proj_linearshifted_2070_2099$StateANSI)
     Data_proj_linearshifted_2070_2099$year=Data_proj_linearshifted_2070_2099$year+2069
     Data_proj_linearshifted_2070_2099$year<-factor(Data_proj_linearshifted_2070_2099$year)
       for (m in 1:(variablenum-1)){
         col_data_proj[m]<-which(colnames(Data_proj_linearshifted_2070_2099)==variablenames[m+1])
       }
     proj<-predict(model,Data_proj_linearshifted_2070_2099)
     for (k in 1:30){
       indx<-which(Data_proj_linearshifted_2070_2099$year==levels(Data_proj_linearshifted_2070_2099$year)[k])
       proj_linearshifted_fit[(i-1)*Pnum+j,k]<-mean(proj[indx],na.rm=TRUE)
         for (n in 1:parasamplenum){
           proj_linearshifted_parasample[(i-1)*Pnum+j,k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data_proj_linearshifted_2070_2099[indx,col_data_proj]))
         }
     }
   } else {
      proj_linearshifted_fit[1, ]<-mean(meanyield_anomaly)
      for (n in 1:parasamplenum){
        hind_parasample[i, ,n]<-mean(meanyield_anomaly)-windowsigma[1]+n/parasamplenum*2*windowsigma[1]
      }
    }
  }
}
 save(proj_linearshifted_fit,file="Metdata/proj_linearshifted_bestfit_2070_2099")
 save(proj_linearshifted_parasample,file="Metdata/proj_linearshifted_parasample_2070_2099")