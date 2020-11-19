
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

#with each climate projection, save the hind/proj yield data
parasamplenum<-5000000
stepwidth<-3 #sample parameters by +- std
hindyear<-32
projyear<-94
modelnum<-length(modelnames)

#matrix and array to save the hindcast/projection samples
#_fit is the best fit (for each structure and each climate), _parasample is the parameter samples
hind_fit<-rep(NA,hindyear)
proj_fit<-matrix(NA,nrow = projyear,ncol=modelnum)
hind_parasample<-matrix(NA,nrow=hindyear,ncol=parasamplenum)

totlen=dim(Data)[1]
meanyield_anomaly<-rep(NA,hindyear)
for (i in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,w = Data$area[indx])
} 



#Uncertainty sampling part
set.seed(666)
Data$GDD_sqr<-Data$GDD_GS^2
Data$EDD_sqr<-Data$EDD_GS^2
Data$Tmax_sqr<-Data$Tmax_GS^2
Data$Tmin_sqr<-Data$Tmin_GS^2
Data$Pr_sqr<-Data$Pr_GS^2
Data$VPD_sqr<-Data$VPD_GS^2
model<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+Tmin_GS+Tmin_sqr+
            Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data) #best model
fullmodel<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmax_GS+Tmax_sqr+Tmin_GS+Tmin_sqr+
                Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data)
variablenames<-variable.names(fullmodel)
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
    
bestestimate<-summary(fullmodel)$coefficient[ ,1]
bestestimatestd<-summary(fullmodel)$coefficient[ ,2]
#EDD terms best estimates from the best model are outside the 3-sigma range of the full model
extrawidth<-c(0,0,0,0.15,0.0005,0,0,0,0,0,0,0,0)
#Latin hypercube sampling in a range for all parameters defined above
MCpara<-randomLHS(parasamplenum,variablenum) #range [0,1]
for (m in 1:variablenum){
  MCpara[ ,m]<-MCpara[ ,m]*2-1 #map to [-1,1]
  MCpara[ ,m]<-MCpara[ ,m]*(bestestimatestd[m]*stepwidth+extrawidth[m]) + bestestimate[m]
}

for (k in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[k])
  hind_fit[k]<-mean(hind[indx],na.rm=TRUE,w = Data$area[indx])
  for (n in 1:parasamplenum){
    hind_parasample[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%
      apply(as.matrix(Data[indx,col_data_hind]),2,function(v) weighted.mean(x=v,w=Data$area[indx]))
  }
}
  
#keep the samples with hindcast passing the plausible window
for (n in 1:parasamplenum){
  if (frac_between(hind_parasample[ ,n],lowerbound,upperbound)!=1){
    hind_parasample[ ,n]<-NA
  }
}
validsample<-which(!is.na(hind_parasample[1, ]))

#projection of each sample
proj_parasample<-array(NA,c(projyear,modelnum,length(validsample)))
col_data_proj<-rep(NA,variablenum-1) #first variable is intercept
for (q in 1:modelnum){
  load(paste("Metdata/macaprojdataframe/Data_",modelnames[q],sep=""))
  Data_proj$StateANSI<-factor(Data_proj$StateANSI)
  Data_proj$year=Data_proj$year+2005
  Data_proj$year<-factor(Data_proj$year)
  Data_proj$GDD_sqr<-Data_proj$GDD_GS^2
  Data_proj$EDD_sqr<-Data_proj$EDD_GS^2
  Data_proj$Tmax_sqr<-Data_proj$Tmax_GS^2
  Data_proj$Tmin_sqr<-Data_proj$Tmin_GS^2
  Data_proj$Pr_sqr<-Data_proj$Pr_GS^2
  Data_proj$VPD_sqr<-Data_proj$VPD_GS^2
  for (m in 1:(variablenum-1)){
    col_data_proj[m]<-which(colnames(Data_proj)==variablenames[m+1])
  }
  proj<-predict(model,Data_proj)
  for (k in 1:projyear){
    indx<-which(Data_proj$year==levels(Data_proj$year)[k])
    proj_fit[k,q]<-mean(proj[indx],na.rm=TRUE)
    proj_parasample[k,q, ]<-MCpara[validsample,1]+MCpara[validsample,2:variablenum]%*%colMeans(as.matrix(Data_proj[indx,col_data_proj]))  
  }
}
hind_parasample<-hind_parasample[ ,validsample]


save(validsample,file="Metdata/kept_sample_num")
save(hind_fit,file="Metdata/hind_bestfit")
save(proj_fit,file="Metdata/proj_bestfit")
save(hind_parasample,file="Metdata/hind_parasample")
save(proj_parasample,file="Metdata/proj_parasample")


