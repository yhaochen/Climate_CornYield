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
load("Metdata/kept_sample_num")
Data$StateANSI<-factor(Data$StateANSI)

#directly load the valid sample index of the 1000000 samples
parasamplenum<-5000000
stepwidth<-3 #sample parameters by +- std
hindyear<-32
projyear<-94

#for the additional linear shifted climate proj (2020-2049; 2070-2099)
proj_linearshifted_fit_2020_2049<-rep(NA,30)
proj_linearshifted_parasample_2020_2049<-matrix(NA,nrow=30,ncol=length(validsample))
proj_linearshifted_fit_2070_2099<-rep(NA,30)
proj_linearshifted_parasample_2070_2099<-matrix(NA,nrow=30,ncol=length(validsample))

#Uncertainty sampling part
set.seed(666)
Data$GDD_sqr<-Data$GDD_GS^2
Data$EDD_sqr<-Data$EDD_GS^2
Data$Tmax_sqr<-Data$Tmax_GS^2
Data$Tmin_sqr<-Data$Tmin_GS^2
Data$Pr_sqr<-Data$Pr_GS^2
Data$VPD_sqr<-Data$VPD_GS^2
model<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmax_GS+Tmax_sqr+Tmin_GS+Tmin_sqr+
            Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data)
variablenames<-variable.names(model)
variablenum<-length(variablenames)
#for hindcast of each parametric sample
col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
for (m in 1:(variablenum-1)){
  col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
}

bestestimate<-summary(model)$coefficient[ ,1]
bestestimatestd<-summary(model)$coefficient[ ,2]
#Latin hypercube sampling in a range for all parameters defined above
MCpara<-randomLHS(parasamplenum,variablenum) #range [0,1]
for (m in 1:variablenum){
  MCpara[ ,m]<-MCpara[ ,m]*2-1 #map to [-1,1]
  MCpara[ ,m]<-MCpara[ ,m]*(bestestimatestd[m]*stepwidth) + bestestimate[m]
}


#linearshifted proj
for (projmodel in 1:2){
  if (projmodel==1){
    col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
    load("Metdata/macaprojdataframe/Data_linearshifted_2020_2049")
    Data_proj_linearshifted_2020_2049$StateANSI<-factor(Data_proj_linearshifted_2020_2049$StateANSI)
    Data_proj_linearshifted_2020_2049$year=Data_proj_linearshifted_2020_2049$year+2019
    Data_proj_linearshifted_2020_2049$year<-factor(Data_proj_linearshifted_2020_2049$year)
    Data_proj_linearshifted_2020_2049$GDD_sqr<-Data_proj_linearshifted_2020_2049$GDD_GS^2
    Data_proj_linearshifted_2020_2049$EDD_sqr<-Data_proj_linearshifted_2020_2049$EDD_GS^2
    Data_proj_linearshifted_2020_2049$Tmax_sqr<-Data_proj_linearshifted_2020_2049$Tmax_GS^2
    Data_proj_linearshifted_2020_2049$Tmin_sqr<-Data_proj_linearshifted_2020_2049$Tmin_GS^2
    Data_proj_linearshifted_2020_2049$Pr_sqr<-Data_proj_linearshifted_2020_2049$Pr_GS^2
    Data_proj_linearshifted_2020_2049$VPD_sqr<-Data_proj_linearshifted_2020_2049$VPD_GS^2
    for (m in 1:(variablenum-1)){
      col_data_proj[m]<-which(colnames(Data_proj_linearshifted_2020_2049)==variablenames[m+1])
    }
    proj<-predict(model,Data_proj_linearshifted_2020_2049)
    for (k in 1:30){
      indx<-which(Data_proj_linearshifted_2020_2049$year==levels(Data_proj_linearshifted_2020_2049$year)[k])
      proj_linearshifted_fit_2020_2049[k]<-mean(proj[indx],na.rm=TRUE)
      for (n in 1:length(validsample)){
        validindx<-validsample[n]
        proj_linearshifted_parasample_2020_2049[k,n]<-MCpara[validindx,1]+MCpara[validindx,2:variablenum]%*%colMeans(as.matrix(Data_proj_linearshifted_2020_2049[indx,col_data_proj]))
      }
    }
  }
  if (projmodel==2){
    col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
    load("Metdata/macaprojdataframe/Data_linearshifted_2070_2099")
    Data_proj_linearshifted_2070_2099$StateANSI<-factor(Data_proj_linearshifted_2070_2099$StateANSI)
    Data_proj_linearshifted_2070_2099$year=Data_proj_linearshifted_2070_2099$year+2069
    Data_proj_linearshifted_2070_2099$year<-factor(Data_proj_linearshifted_2070_2099$year)
    Data_proj_linearshifted_2070_2099$GDD_sqr<-Data_proj_linearshifted_2070_2099$GDD_GS^2
    Data_proj_linearshifted_2070_2099$EDD_sqr<-Data_proj_linearshifted_2070_2099$EDD_GS^2
    Data_proj_linearshifted_2070_2099$Tmax_sqr<-Data_proj_linearshifted_2070_2099$Tmax_GS^2
    Data_proj_linearshifted_2070_2099$Tmin_sqr<-Data_proj_linearshifted_2070_2099$Tmin_GS^2
    Data_proj_linearshifted_2070_2099$Pr_sqr<-Data_proj_linearshifted_2070_2099$Pr_GS^2
    Data_proj_linearshifted_2070_2099$VPD_sqr<-Data_proj_linearshifted_2070_2099$VPD_GS^2
    for (m in 1:(variablenum-1)){
      col_data_proj[m]<-which(colnames(Data_proj_linearshifted_2070_2099)==variablenames[m+1])
    }
    proj<-predict(model,Data_proj_linearshifted_2070_2099)
    for (k in 1:30){
      indx<-which(Data_proj_linearshifted_2070_2099$year==levels(Data_proj_linearshifted_2070_2099$year)[k])
      proj_linearshifted_fit_2070_2099[k]<-mean(proj[indx],na.rm=TRUE)
      for (n in 1:length(validsample)){
        validindx<-validsample[n]
        proj_linearshifted_parasample_2070_2099[k,n]<-MCpara[validindx,1]+MCpara[validindx,2:variablenum]%*%colMeans(as.matrix(Data_proj_linearshifted_2070_2099[indx,col_data_proj]))
      }
    }
  } 
}
save(proj_linearshifted_fit_2020_2049,file="Metdata/proj_linearshifted_bestfit_2020_2049")
save(proj_linearshifted_parasample_2020_2049,file="Metdata/proj_linearshifted_parasample_2020_2049")
save(proj_linearshifted_fit_2070_2099,file="Metdata/proj_linearshifted_bestfit_2070_2099")
save(proj_linearshifted_parasample_2070_2099,file="Metdata/proj_linearshifted_parasample_2070_2099")