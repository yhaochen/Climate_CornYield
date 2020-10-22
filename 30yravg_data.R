#This script calculates 30-year mean yield for two time windows (2020-2049 and 2070-2099) under different uncertainty samples
rm(list = ls())
graphics.off()

#read hindcast/projection data and calculate 30-year mean for two time windows (2020-2049 and 2070-2099)
load("Metdata/hind_bestfit")
load("Metdata/proj_bestfit")
load("Metdata/hind_parasample")
load("Metdata/proj_parasample")
totalyears<-dim(proj_parasample)[1]
climnum<-dim(proj_parasample)[2]+1
parasamplenum<-dim(proj_parasample)[3]
for (k in 1:2){
  if (k==1){
    selectedyears<-c(15:44) # 2020-2049
    load("Metdata/proj_linearshifted_bestfit_2020_2049")
    load("Metdata/proj_linearshifted_parasample_2020_2049")
    proj_linearshifted_parasample<-proj_linearshifted_parasample_2020_2049
    proj_linearshifted_fit<-proj_linearshifted_fit_2020_2049
  }
  if (k==2){
    selectedyears<-c(65:94) # 2070-2099
    load("Metdata/proj_linearshifted_bestfit_2070_2099")
    load("Metdata/proj_linearshifted_parasample_2070_2099")
    proj_linearshifted_parasample<-proj_linearshifted_parasample_2070_2099
    proj_linearshifted_fit<-proj_linearshifted_fit_2070_2099
  }

  # set up the 3_D matrix of year*clim*para samples
  proj_all<-array(NA,c(30,climnum,parasamplenum))
  proj_all[ ,c(1:18), ]<-proj_parasample[selectedyears, , ]
  proj_all[ ,19, ]<-proj_linearshifted_parasample
  proj_30yravg<-apply(proj_all,c(2,3),mean,na.rm=TRUE) #30-yr average for each combination

  
  proj_fit_all<-matrix(NA,nrow=30,ncol=climnum)
  proj_fit_all[ ,c(1:18)]<-proj_fit[selectedyears, ]
  proj_fit_all[ ,19]<-proj_linearshifted_fit
  proj_fit_30yravg<-colMeans(proj_fit_all,na.rm=TRUE)
  if (k==1){
    save(proj_fit_30yravg,file="Metdata/proj_fit_30yravg_2020_2049")
    save(proj_30yravg,file="Metdata/proj_30yravg_2020_2049")
  }
  if (k==2){
    save(proj_fit_30yravg,file="Metdata/proj_fit_30yravg_2070_2099")
    save(proj_30yravg,file="Metdata/proj_30yravg_2070_2099")
  }
}