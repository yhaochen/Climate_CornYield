
#Take the mean climate projection (from 18 models), and combine with hindcast to form 1980-2099 data
rm(list = ls())
graphics.off()
source("plots.R")
source("GDDEDD.R")

modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM",
              "MIROC-ESM-CHEM", "inmcm4", "HadGEM2-CC365", "CanESM2", "bcc-csm1-1-m")

#calculate a mean climate proj/hind

#projection period: 2006-2099
Tmax_mean_proj<-matrix(0,nrow=1878,ncol=71*365+23*366)
Tmin_mean_proj<-matrix(0,nrow=1878,ncol=71*365+23*366)
Pr_mean_proj<-matrix(0,nrow=1878,ncol=71*365+23*366)
Rhmax_mean_proj<-matrix(0,nrow=1878,ncol=71*365+23*366)
Rhmin_mean_proj<-matrix(0,nrow=1878,ncol=71*365+23*366)

#read the projection data
for (i in 1:length(modelnames)){
   foldername<-paste("SourceData/MACAv2-METDATA_proj/",modelnames[i],"_proj/",sep="")
   load(paste(foldername,"tmax",sep=""))
   load(paste(foldername,"tmin",sep=""))
   load(paste(foldername,"pr",sep=""))
   load(paste(foldername,"rhmax",sep=""))
   load(paste(foldername,"rhmin",sep=""))
   Tmax_mean_proj<-Tmax_mean_proj+tmax
   Tmin_mean_proj<-Tmin_mean_proj+tmin
   Pr_mean_proj<-Pr_mean_proj+pr
   Rhmax_mean_proj<-Rhmax_mean_proj+rhmax
   Rhmin_mean_proj<-Rhmin_mean_proj+rhmin
   print(i)
   print(which(is.na(Pr_mean_proj[1, ]))[1])
}
#Calculate projection mean
Tmax_mean_proj<-Tmax_mean_proj/length(modelnames)-273.15
Tmin_mean_proj<-Tmin_mean_proj/length(modelnames)-273.15
Pr_mean_proj<-Pr_mean_proj/length(modelnames)
Rhmax_mean_proj<-Rhmax_mean_proj/length(modelnames)
Rhmin_mean_proj<-Rhmin_mean_proj/length(modelnames)
#Save projection mean data
save(Tmax_mean_proj,file="SourceData/MACAv2-METDATA_proj/mean_proj/Tmax")
save(Tmin_mean_proj,file="SourceData/MACAv2-METDATA_proj/mean_proj/Tmin")
save(Pr_mean_proj,file="SourceData/MACAv2-METDATA_proj/mean_proj/Pr")
save(Rhmax_mean_proj,file="SourceData/MACAv2-METDATA_proj/mean_proj/Rhmax")
save(Rhmin_mean_proj,file="SourceData/MACAv2-METDATA_proj/mean_proj/Rhmin")


#hindcast part (1980-2005)   
Tmax_mean_hind<-matrix(0,nrow=1878,ncol=19*365+7*366)
Tmin_mean_hind<-matrix(0,nrow=1878,ncol=19*365+7*366)
Pr_mean_hind<-matrix(0,nrow=1878,ncol=19*365+7*366)
Rhmax_mean_hind<-matrix(0,nrow=1878,ncol=19*365+7*366)
Rhmin_mean_hind<-matrix(0,nrow=1878,ncol=19*365+7*366)
#read hindcast data
for (i in 1:length(modelnames)){
  foldername<-paste("SourceData/MACAv2-METDATA_hind/",modelnames[i],"_hind/",sep="")
  load(paste(foldername,"tmax",sep=""))
  load(paste(foldername,"tmin",sep=""))
  load(paste(foldername,"pr",sep=""))
  load(paste(foldername,"rhmax",sep=""))
  load(paste(foldername,"rhmin",sep=""))
  Tmax_mean_hind<-Tmax_mean_hind+tmax
  Tmin_mean_hind<-Tmin_mean_hind+tmin
  Pr_mean_hind<-Pr_mean_hind+pr
  Rhmax_mean_hind<-Rhmax_mean_hind+rhmax
  Rhmin_mean_hind<-Rhmin_mean_hind+rhmin
  print(i)
  print(which(is.na(Pr_mean_hind[1, ]))[1])
}
#calculate hindcast mean
Tmax_mean_hind<-Tmax_mean_hind/length(modelnames)-273.15
Tmin_mean_hind<-Tmin_mean_hind/length(modelnames)-273.15
Pr_mean_hind<-Pr_mean_hind/length(modelnames)
Rhmax_mean_hind<-Rhmax_mean_hind/length(modelnames)
Rhmin_mean_hind<-Rhmin_mean_hind/length(modelnames)
#Save hindcast mean
save(Tmax_mean_hind,file="SourceData/MACAv2-METDATA_hind/mean_hind/Tmax")
save(Tmin_mean_hind,file="SourceData/MACAv2-METDATA_hind/mean_hind/Tmin")
save(Pr_mean_hind,file="SourceData/MACAv2-METDATA_hind/mean_hind/Pr")
save(Rhmax_mean_hind,file="SourceData/MACAv2-METDATA_hind/mean_hind/Rhmax")
save(Rhmin_mean_hind,file="SourceData/MACAv2-METDATA_hind/mean_hind/Rhmin")


#Relative humidity and temperature mean
RHmean_proj<-(Rhmax_mean_proj+Rhmin_mean_proj)/2
Tmean_proj<-(Tmax_mean_proj+Tmin_mean_proj)/2
RHmean_hind<-(Rhmax_mean_hind+Rhmin_mean_hind)/2
Tmean_hind<-(Tmax_mean_hind+Tmin_mean_hind)/2

#combine hindcast+projection data (mean of 18 models)
Tmean_all<-cbind(Tmean_hind,Tmean_proj)
RHmean_all<-cbind(RHmean_hind,RHmean_proj)
Pr_mean_all<-cbind(Pr_mean_hind,Pr_mean_proj)
save(Tmean_all,file="SourceData/MACAv2-METDATA_meandata_1980-2099/T")
save(RHmean_all,file="SourceData/MACAv2-METDATA_meandata_1980-2099/RH")
save(Pr_mean_all,file="SourceData/MACAv2-METDATA_meandata_1980-2099/Pr")
