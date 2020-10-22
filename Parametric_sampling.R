
#Sample the parametric uncertainty
rm(list = ls())
graphics.off()
library(lhs)

frac_between <- function (vect1,vect2,vect3){ #vect2 is lower, vect3 is higher
  totlen<-length(vect1)
  out<-0
  for (i in 1:totlen){
    if ((vect1[i]<vect2[i]) | (vect1[i]>vect3[i])){
      out<-out+1
    }
  }
  return((totlen-out)/totlen)
}


#18 models
modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM",
              "MIROC-ESM-CHEM", "inmcm4", "HadGEM2-CC365", "CanESM2", "bcc-csm1-1-m")

#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)

#with each climate projection, save the hind/proj yield data of 64 structures
parasamplenum<-100000
hindyear<-32
projyear<-94
modelnum<-length(modelnames)

#matrix and array to save the hindcast/projection samples
#_fit is the best fit (for each structure and each climate), _parasample is the parameter samples
hind_fit<-rep(NA,hindyear)
proj_fit<-matrix(NA,nrow = projyear,ncol=modelnum)
hind_parasample<-matrix(NA,nrow=hindyear,ncol=parasamplenum)
proj_parasample<-array(NA,c(projyear,modelnum,parasamplenum))
totlen=dim(Data)[1]
meanyield_anomaly<-rep(NA,hindyear)
for (i in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,w = Data$area[indx])
} 



#Uncertainty sampling part
set.seed(666)
model<-lm(yield_anomaly~Tmax_GS+Tmin_GS+GDD_GS+EDD_GS+Pr_GS+VPD_GS,data=Data)
variablenames<-variable.names(model)
variablenum<-length(variablenames)
#for hindcast of each parametric sample
col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
for (m in 1:(variablenum-1)){
  col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
}
hind<-predict(model,Data)
#use model's best estimate to find a plausible window that covers 95% annual yield observation
Besthind<-data.frame(yield_anomaly=hind,area=Data$area)
hind_yield_anomaly<-rep(NA,hindyear)
for (k in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[k])
  hind_yield_anomaly[k]<-weighted.mean(Besthind$yield_anomaly[indx],na.rm=T,w = Besthind$area[indx])
} 
difference<-abs(hind_yield_anomaly-meanyield_anomaly)
step<-quantile(difference,0.95)
upperbound<-hind_yield_anomaly+step
lowerbound<-hind_yield_anomaly-step
    
bestestimate<-summary(model)$coefficient[ ,1]
#Latin hypercube sampling in a range for all parameters defined above
MCpara<-randomLHS(parasamplenum,variablenum) #range [0,1]
stepwidth<-0.25 #sample parameters by +- certain percent
for (m in 1:variablenum){
  MCpara[ ,m]<-MCpara[ ,m]*bestestimate[m]*stepwidth*2+bestestimate[m]*(1-stepwidth)
}

for (k in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[k])
  hind_fit[k]<-mean(hind[indx],na.rm=TRUE,w = Data$area[indx])
  for (n in 1:parasamplenum){
    hind_parasample[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
  }
}
  
#keep the samples with hindcast passing the plausible window
for (n in 1:parasamplenum){
  if (frac_between(hind_parasample[ ,n],lowerbound,upperbound)!=1){
    hind_parasample[ ,n]<-NA
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
    proj_fit[k,q]<-mean(proj[indx],na.rm=TRUE,w = Data$area[indx])
    for (n in 1:parasamplenum){
      if (!is.na(hind_parasample[1,n])){
        proj_parasample[k,q,n]<-MCpara[n,1]+
          MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data_proj[indx,col_data_proj]))
      }
    }
  }
}
validsample<-which(!is.na(hind_parasample[1, ]))
hind_parasample<-hind_parasample[ ,validsample]
proj_parasample<-proj_parasample[ , ,validsample]
save(hind_fit,file="Metdata/hind_bestfit")
save(proj_fit,file="Metdata/proj_bestfit")
save(hind_parasample,file="Metdata/hind_parasample")
save(proj_parasample,file="Metdata/proj_parasample")


