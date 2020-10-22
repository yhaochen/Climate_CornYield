#Sample the parametric uncertainty for two linear shifted climate scenario
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

#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)

#with each climate projection, save the hind/proj yield data of 64 structures
parasamplenum<-100000
hindyear<-32
projyear<-94

#for the additional linear shifted climate proj (2020-2049; 2070-2099)
proj_linearshifted_fit_2020_2049<-rep(NA,30)
proj_linearshifted_parasample_2020_2049<-matrix(NA,nrow=30,ncol=parasamplenum)
proj_linearshifted_fit_2070_2099<-rep(NA,30)
proj_linearshifted_parasample_2070_2099<-matrix(NA,nrow=30,ncol=parasamplenum)
hind_parasample<-matrix(NA,nrow=hindyear,ncol=parasamplenum)
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
     
#linearshifted proj
for (projmodel in 1:2){
  if (projmodel==1){
    col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
    load("Metdata/macaprojdataframe/Data_linearshifted_2020_2049")
    Data_proj_linearshifted_2020_2049$StateANSI<-factor(Data_proj_linearshifted_2020_2049$StateANSI)
    Data_proj_linearshifted_2020_2049$year=Data_proj_linearshifted_2020_2049$year+2019
    Data_proj_linearshifted_2020_2049$year<-factor(Data_proj_linearshifted_2020_2049$year)
    for (m in 1:(variablenum-1)){
      col_data_proj[m]<-which(colnames(Data_proj_linearshifted_2020_2049)==variablenames[m+1])
    }
    proj<-predict(model,Data_proj_linearshifted_2020_2049)
    for (k in 1:30){
      indx<-which(Data_proj_linearshifted_2020_2049$year==levels(Data_proj_linearshifted_2020_2049$year)[k])
      proj_linearshifted_fit_2020_2049[k]<-mean(proj[indx],na.rm=TRUE)
      for (n in 1:parasamplenum){
        if (!is.na(hind_parasample[1,n])){
          proj_linearshifted_parasample_2020_2049[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data_proj_linearshifted_2020_2049[indx,col_data_proj]))
        }
      }
    }
  }
  if (projmodel==2){
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
      proj_linearshifted_fit_2070_2099[k]<-mean(proj[indx],na.rm=TRUE)
      for (n in 1:parasamplenum){
        if (!is.na(hind_parasample[1,n])){
          proj_linearshifted_parasample_2070_2099[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data_proj_linearshifted_2070_2099[indx,col_data_proj]))
        } 
      }
    }
  } 
}
validsample<-which(!is.na(hind_parasample[1, ]))
proj_linearshifted_parasample_2020_2049<-proj_linearshifted_parasample_2020_2049[ ,validsample]
proj_linearshifted_parasample_2070_2099<-proj_linearshifted_parasample_2070_2099[ ,validsample]
save(proj_linearshifted_fit_2020_2049,file="Metdata/proj_linearshifted_bestfit_2020_2049")
save(proj_linearshifted_parasample_2020_2049,file="Metdata/proj_linearshifted_parasample_2020_2049")
save(proj_linearshifted_fit_2070_2099,file="Metdata/proj_linearshifted_bestfit_2070_2099")
save(proj_linearshifted_parasample_2070_2099,file="Metdata/proj_linearshifted_parasample_2070_2099")