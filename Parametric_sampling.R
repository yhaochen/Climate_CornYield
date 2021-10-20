#Sample the parametric uncertainty
rm(list = ls())
graphics.off()
library(lhs) #Latin Hypercube sampling package
library(parallel)
library(MASS)

#this function calculates how much of vect1 is between vect2 and vect3
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

#this function calculates the hindcasts based on sampled parameters and historical climate
hindcalc<-function(MC){
  result<-variables %*% MC
  if (frac_between(result,lowerbound,upperbound)!=1){
    return (rep(NA,hindyear))
  } else{
    return (result)
  }
}

#18 models
modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM",
              "MIROC-ESM-CHEM", "inmcm4", "HadGEM2-CC365", "CanESM2", "bcc-csm1-1-m")

#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
#sampling number in each run (too large would cause memory issue)
parasamplenum<-10000000
stepwidth<-c(10,10,10,50,50,10,10,10,10,10,10,10,10) #sample parameters range: defined by +- standard deviations
hindyear<-40
projyear<-94
modelnum<-length(modelnames)

dir.create("Metdata/samples/valid_samples",recursive = TRUE)
dir.create("Metdata/samples/proj_parasample",recursive = TRUE)
dir.create("Metdata/samples/hind_parasample",recursive = TRUE)

#_fit is the best fit (for each year each climate projection)
hind_fit<-rep(NA,hindyear)
proj_fit<-matrix(NA,nrow = projyear,ncol=modelnum)
totlen=dim(Data)[1]

#historical mean yield anomaly each year
meanyield_anomaly<-rep(NA,hindyear)
for (i in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,w = Data$area[indx])
} 
#best model and full model regression
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
#which column of Data corresponds to each parameter
col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
for (m in 1:(variablenum-1)){
  col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
}
hind<-predict(model,Data)
#use model's best estimate to find a plausible window that covers 95% annual yield observation
Besthind<-data.frame(yield_anomaly=hind,area=Data$area)
hind_fit<-rep(NA,hindyear)

#calculate the annual variables in each year to save calculation time, and save the best hindcast fit
variables<-data.frame(Intercept=rep(1,hindyear),GDD=rep(NA,hindyear),GDD_sqr=rep(NA,hindyear),EDD=rep(NA,hindyear),EDD_sqr=rep(NA,hindyear),
                      Tmax=rep(NA,hindyear),Tmax_sqr=rep(NA,hindyear),Tmin=rep(NA,hindyear),Tmin_sqr=rep(NA,hindyear),
                      Pr=rep(NA,hindyear),Pr_sqr=rep(NA,hindyear),VPD=rep(NA,hindyear),VPD_sqr=rep(NA,hindyear))
variables<-matrix(NA,nrow=40,ncol=13)
variables[ ,1]<-1
for (k in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[k])
  hind_fit[k]<-weighted.mean(Besthind$yield_anomaly[indx],na.rm=T,w = Besthind$area[indx])
  Datatoget<-Data[indx,col_data_hind[c(1,3,5,7,9,11)]]
  variables[k,c(2,4,6,8,10,12)]<-apply(Datatoget,2,function(v) weighted.mean(x=v,w=Data$area[indx]))
  variables[k,c(3,5,7,9,11,13)]<-variables[k,c(2,4,6,8,10,12)]^2
} 
save(hind_fit,file="Metdata/hind_bestfit")

difference<-abs(hind_fit-meanyield_anomaly)
step<-quantile(difference,0.95)
upperbound<-hind_fit+step
lowerbound<-hind_fit-step
bestestimate<-summary(fullmodel)$coefficient[ ,1]
bestestimatestd<-summary(fullmodel)$coefficient[ ,2]
save(step,file="Metdata/step")

#Uncertainty sampling part
#Set different seeds in each run
for (trial in 405:456){
  set.seed(trial)
  #Latin hypercube sampling in a range for all parameters defined above
  MCpara<-randomLHS(parasamplenum,variablenum) #range [0,1]
  
  for (m in 1:variablenum)  {
    MCpara[ ,m]<-MCpara[ ,m]*2-1 #map to [-1,1]
    MCpara[ ,m]<-MCpara[ ,m]*(bestestimatestd[m]*stepwidth[m]) + bestestimate[m]
  }
  MCparalist<-lapply(seq_len(nrow(MCpara)), function(i) MCpara[i, ]) #turn into list
  #the parallel computation part
  hindtoget<-mcmapply(hindcalc,MCparalist,mc.cores = 4)
  #find the hindcasts that pass the plausible window
  validsample<-which(!is.na(hindtoget[1, ]))
  hind_parasample<-hindtoget[ ,validsample]
  
  #projection of each accepted sample
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
      #in the first run calculate the projection fit and save it
      if (trial==1){
        proj_fit[k,q]<-mean(proj[indx],na.rm=TRUE)
      }
      proj_parasample[k,q, ]<-MCpara[validsample,1]+MCpara[validsample,2:variablenum]%*%colMeans(as.matrix(Data_proj[indx,col_data_proj]))  
    }
  }
  MCpara<-MCparalist[validsample]
  save(MCpara,file=paste("Metdata/samples/valid_samples/para",trial,sep=""))
  save(hind_parasample,file=paste("Metdata/samples/hind_parasample/hind",trial,sep=""))
  save(proj_parasample,file=paste("Metdata/samples/proj_parasample/proj",trial,sep=""))
  if (trial==1){
    save(proj_fit,file="Metdata/proj_bestfit")
  }
}

#combind all files into a large one
 num<-c(1:1000)
 library(abind)
 for (i in 1:456){
   paralistname<-paste("Metdata/samples/valid_samples/para",num[i],sep="")
   load(paralistname)
   if (i==1){
     Para<-MCpara
   } else{
     Para<-append(Para,MCpara)
   }
   
   hindname<-paste("Metdata/samples/hind_parasample/hind",num[i],sep="")
   load(hindname)
   if (i==1){
     Hind<-hind_parasample
   } else{
     Hind<-cbind(Hind,hind_parasample)
   }
   
   projname<-paste("Metdata/samples/proj_parasample/proj",num[i],sep="")
   load(projname)
   if (i==1){
     Proj<-proj_parasample
   } else{
     Proj<-abind(Proj,proj_parasample,along=3)
   }
 }

 
 for (i in 457:1000){
   paralistname<-paste("Metdata/samples/valid_samples/para",num[i],sep="")
   load(paralistname)
   if (i==1){
     Para<-MCpara
   } else{
     Para<-append(Para,MCpara)
   }
   hind_parasample<-variables %*% matrix(unlist(MCpara),nrow=13,byrow=FALSE)
   Hind<-cbind(Hind,hind_parasample)
   
   
   projname<-paste("Metdata/samples/proj_parasample/proj",num[i],sep="")
   load(projname)
   if (i==1){
     Proj<-proj_parasample
   } else{
     Proj<-abind(Proj,proj_parasample,along=3)
   }
 }
 save(Para,file="Metdata/valid_samples")
 hind_parasample<-Hind
 proj_parasample<-Proj
 save(hind_parasample,file="Metdata/hind_parasample")
 save(proj_parasample,file="Metdata/proj_parasample")
 
 
 
 for (i in 1:32){
   indx<-which(Data$year==levels(Data$year)[i+2])
   if (i==1){
     old_indx<-indx
   }
   old_indx<-append(old_indx,indx)
 }

 old_model<-model<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+Tmin_GS+Tmin_sqr+
                        Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data[old_indx, ]) 
 old_fit<-c(-0.8699407,3.8729372,-6.9691059,1.0510181,1.4837100,-4.2384899,
            -7.6086790,-15.1707768,3.6049754,5.6470485,-0.3883320,6.3727977,
            -0.3235119,3.0485193,1.9649107,3.5188601,5.3067566,3.5163875,
            1.3933228,6.0683284,5.9629284,-1.4350786,5.9281361,8.4526365,
            3.5307753,6.3300494,-2.6024171,6.9728213,11.7803955,-4.3948945,
            -20.7319639,-25.9410716)
 old_proj<-predict(old_model,Data)
 old_proj_annual<-rep(NA,8)
 for (i in 1:8){
   j<-i
   if (i>2){
     j<-i+32
   }
   indx<-which(Data$year==levels(Data$year)[j])
   old_proj_annual[i]<-weighted.mean(old_proj[indx],na.rm=T,w = Data$area[indx])
 }
 years<-c(1979:2018)
 par(mar=c(4,5.1,1.6,2.1))
 plot(years,hind_fit,type="l",xlab="Year",ylab="Yield anomaly (bu/acre)",cex.axis=1.8,cex.lab=1.8,lwd=2)
 points(years,hind_fit,pch=20)
 lines(c(1981:2012),old_fit,type = "l",col="blue")
 points(c(1981:2012),old_fit,pch=20,col="blue",lwd=2)
 points(c(1979:1980,2013:2018),old_proj_annual,col="red",pch=20)
 legend("bottomleft",col=c("black","blue","red"),lty=c(1,1,NA),lwd=c(2,2,NA),pch=c(NA,NA,20),bty="n",
        legend=c("Best fit with 40-year data","Best fit with 32-year data"),cex=1.8)
 