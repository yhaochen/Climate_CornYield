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


#18 models
modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM",
              "MIROC-ESM-CHEM", "inmcm4", "HadGEM2-CC365", "CanESM2", "bcc-csm1-1-m")

#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
Data$year=Data$year+1978
Data$year<-factor(Data$year)

#with each climate projection, save the hind/proj yield data of 64 structures
parasamplenum<-10000
hindyear<-32
projyear<-94
modelnum<-length(modelnames)
Tstructure<-c(5,6,7,8)
Pstructure<-c(9,10)
Tnum<-2^length(Tstructure)
Pnum<-2^length(Pstructure)
strunum<-Tnum*Pnum

hind_fit<-matrix(NA,nrow = strunum,ncol=hindyear)
proj_fit<-matrix(NA,nrow = strunum,ncol=projyear*modelnum)
hind_parasample<-array(NA,c(strunum,hindyear,parasamplenum))
proj_parasample<-array(NA,c(strunum,projyear*modelnum,parasamplenum))

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

meanyield_anomaly<-rep(NA,hindyear)
for (i in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,w = Data$area[indx])
}
window_sigma<-rep(NA,strunum)
#Uncertainty sampling part
set.seed(666)
for (i in 1:Tnum){ #T
  for (j in 1:Pnum){ #P
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
    if (variablenum!=1){
      #for hindcast of each parametric sample
      col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
      for (m in 1:(variablenum-1)){
        col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
      }
    
    hind<-predict(model,Data)
    
    covered_percent<-rep(NA,300)
    hind_parasample1<-rep(NA,parasamplenum)
    for (q in 1:300){
      sigma<-q*0.1
      MCpara<-parasample(model,parasamplenum,sigma)
      included<-rep(NA,hindyear)
      for (k in 1:hindyear){
        indx<-which(Data$year==levels(Data$year)[k])
        for (n in 1:parasamplenum){
          hind_parasample1[n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
        }
        included[k]<-(meanyield_anomaly[k]<max(hind_parasample1)) && (meanyield_anomaly[k]>min(hind_parasample1))
      }
      covered_percent[q]<-sum(included)/hindyear
    }
    window_sigma[(i-1)*Pnum+j]<-which.min(abs(covered_percent-0.95))*0.1
    MCpara<-parasample(model,parasamplenum,window_sigma[(i-1)*Pnum+j])
    
     for (k in 1:hindyear){
       indx<-which(Data$year==levels(Data$year)[k])
       hind_fit[(i-1)*Pnum+j,k]<-mean(hind[indx],na.rm=TRUE)
       if (variablenum!=1) {
         for (n in 1:parasamplenum){
           hind_parasample[(i-1)*Pnum+j,k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
         }
       }
     }
    
    #projection of each sample
    col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
    for (q in 1:modelnum){
      load(paste("Metdata/macaprojdataframe/Data_",modelnames[q],sep=""))
      Data_proj$StateANSI<-factor(Data_proj$StateANSI)
      Data_proj$year=Data_proj$year+2005
      Data_proj$year<-factor(Data_proj$year)
      for (m in 1:(variablenum-1)){
        col_data_proj[m]<-which(colnames(Data_proj)==variablenames[m+1])
      }
      proj<-predict(model,Data_proj)
      for (k in 1:projyear){
        indx<-which(Data_proj$year==levels(Data_proj$year)[k])
        proj_fit[(i-1)*Pnum+j,(q-1)*projyear+k]<-mean(proj[indx],na.rm=TRUE)
        for (n in 1:parasamplenum){
          proj_parasample[(i-1)*Pnum+j,(q-1)*projyear+k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data_proj[indx,col_data_proj]))
        }
      }
    }
  
    } else{ # the base model: no weather variables
      hind_fit[1, ]<-mean(meanyield_anomaly)
      proj_fit[1, ]<-mean(meanyield_anomaly)
      windowwidth<-sort(abs(meanyield_anomaly-mean(meanyield_anomaly)))[hindyear-1]
      window_sigma[1]<-windowwidth
      for (n in 1:parasamplenum){
        hind_parasample[i, ,n]<-mean(meanyield_anomaly)-windowwidth+n/parasamplenum*2*windowwidth
      }
    }
    print((i-1)*Pnum+j)
  }
}
save(window_sigma,file="Metdata/window_sigma")
save(hind_fit,file="Metdata/hind_bestfit")
save(proj_fit,file="Metdata/proj_bestfit")
save(hind_parasample,file="Metdata/hind_parasample")
save(proj_parasample,file="Metdata/proj_parasample")


