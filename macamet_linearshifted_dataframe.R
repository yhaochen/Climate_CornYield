
rm(list = ls())
graphics.off()
source("plots.R")
source("GDDEDD.R")

ma <- function(arr, n){  #moving avg function
  res = rep(NA,length(arr)-n+1)
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n+1):i])
  }
  res
}

#load combined hindcast + mean climate data and use them to calculate temperature difference, pr and RH ratio
load("SourceData/MACAv2-METDATA_meandata_1980-2099/T")
load("SourceData/MACAv2-METDATA_meandata_1980-2099/RH")
load("SourceData/MACAv2-METDATA_meandata_1980-2099/Pr")
#load Met observation data and shift based on the difference between mean projection data and hindcast data
load("Metdata/Metdataframe/Mettmax")
load("Metdata/Metdataframe/Mettmin")
load("Metdata/Metdataframe/Metpr")
load("Metdata/Metdataframe/Metrhmax")
load("Metdata/Metdataframe/Metrhmin")
totdays<-30*365+8
countynum<-1878

#calculate monthly weather difference (1980-2009 VS 2020-2049 2070-2099)
days_month<-c(0,31,29,31,30,31,30,31,31,30,31,30,31)
Tdiff_2020<-rep(NA,12)
Prratio_2020<-rep(NA,12)
Rhratio_2020<-rep(NA,12)
Tdiff_2070<-rep(NA,12)
Prratio_2070<-rep(NA,12)
Rhratio_2070<-rep(NA,12)
indx_2020<-40*365+10
indx_2070<-90*365+23-1

#shifted weather data
shiftedMettmax_2020_2049<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMettmin_2020_2049<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMetpr_2020_2049<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMetrhmax_2020_2049<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMetrhmin_2020_2049<-matrix(NA,nrow = 1878,ncol = totdays)

shiftedMettmax_2070_2099<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMettmin_2070_2099<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMetpr_2070_2099<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMetrhmax_2070_2099<-matrix(NA,nrow = 1878,ncol = totdays)
shiftedMetrhmin_2070_2099<-matrix(NA,nrow = 1878,ncol = totdays)

for (i in 1:12){
  m1<-sum(days_month[c(1:i)])+1
  for (j in 1:30){ #each 30-yr period: 1980-2009, 2020-2049, 2070-2099
    if ( ((j%%4==1) & (i<=2)) || ((j%%4==0) & (i>2)) ) {  
      days_year<-366
    } else{
      days_year<-365
    }
    #find the days index of each month
    if (j==1){
      m2<-m1+days_month[i+1]-1
      daysindx_1980<-c(m1:m2)
      daysindx_2020<-c(m1:m2)+indx_2020
      daysindx_2070<-c(m1:m2)+indx_2070
    }
    if (j>1){
      daysindx_1980<-append(daysindx_1980,c(m1:m2),after=length(daysindx_1980))
      daysindx_2020<-append(daysindx_2020,c(m1:m2)+indx_2020,after=length(daysindx_2020))
      daysindx_2070<-append(daysindx_2070,c(m1:m2)+indx_2070,after=length(daysindx_2070))
    }
    m1<-m1+days_year
    m2<-m2+days_year
    if ((i==2)&(j%%4==1)){
      m2=m2-1
    }
    if ((i==2)&(j%%4==0)){
      m2=m2+1
    }
  }
  #calculate T difference and P, RH ratio
  Tdiff_2020[i]<-mean(Tmean_all[ ,daysindx_2020])-mean(Tmean_all[ ,daysindx_1980])
  Tdiff_2070[i]<-mean(Tmean_all[ ,daysindx_2070])-mean(Tmean_all[ ,daysindx_1980])
  Prratio_2020[i]<-sum(Pr_mean_all[ ,daysindx_2020])/sum(Pr_mean_all[ ,daysindx_1980])
  Prratio_2070[i]<-sum(Pr_mean_all[ ,daysindx_2070])/sum(Pr_mean_all[ ,daysindx_1980])
  Rhratio_2020[i]<-mean(RHmean_all[ ,daysindx_2020])/mean(RHmean_all[ ,daysindx_1980])
  Rhratio_2070[i]<-mean(RHmean_all[ ,daysindx_2070])/mean(RHmean_all[ ,daysindx_1980])
  #Shift data in two 30-yr period
  shiftedMettmax_2020_2049[ ,daysindx_1980]<-tmax[ ,(daysindx_1980+730)]+Tdiff_2020[i]
  shiftedMettmin_2020_2049[ ,daysindx_1980]<-tmin[ ,(daysindx_1980+730)]+Tdiff_2020[i]
  shiftedMetpr_2020_2049[ ,daysindx_1980]<-pr[ ,(daysindx_1980+730)]*Prratio_2020[i]
  shiftedMetrhmax_2020_2049[ ,daysindx_1980]<-rhmax[ ,(daysindx_1980+730)]*Rhratio_2020[i]
  shiftedMetrhmin_2020_2049[ ,daysindx_1980]<-rhmin[ ,(daysindx_1980+730)]*Rhratio_2020[i]
  
  shiftedMettmax_2070_2099[ ,daysindx_1980]<-tmax[ ,(daysindx_1980+730)]+Tdiff_2070[i]
  shiftedMettmin_2070_2099[ ,daysindx_1980]<-tmin[ ,(daysindx_1980+730)]+Tdiff_2070[i]
  shiftedMetpr_2070_2099[ ,daysindx_1980]<-pr[ ,(daysindx_1980+730)]*Prratio_2070[i]
  shiftedMetrhmax_2070_2099[ ,daysindx_1980]<-rhmax[ ,(daysindx_1980+730)]*Rhratio_2070[i]
  shiftedMetrhmin_2070_2099[ ,daysindx_1980]<-rhmin[ ,(daysindx_1980+730)]*Rhratio_2070[i]
}

#Save the shifted data
save(shiftedMettmax_2020_2049,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/tmax_2020_2049")
save(shiftedMettmin_2020_2049,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/tmin_2020_2049")
save(shiftedMetpr_2020_2049,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/pr_2020_2049")
save(shiftedMetrhmax_2020_2049,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmax_2020_2049")
save(shiftedMetrhmin_2020_2049,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmin_2020_2049")
save(shiftedMettmax_2070_2099,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/tmax_2070_2099")
save(shiftedMettmin_2070_2099,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/tmin_2070_2099")
save(shiftedMetpr_2070_2099,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/pr_2070_2099")
save(shiftedMetrhmax_2070_2099,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmax_2070_2099")
save(shiftedMetrhmin_2070_2099,file="SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmin_2070_2099")


#can also start from loading after saving 
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/tmax_2020_2049")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/tmin_2020_2049")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/pr_2020_2049")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmax_2020_2049")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmin_2020_2049")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/tmax_2070_2099")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/tmin_2070_2099")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/pr_2070_2099")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmax_2070_2099")
# load("SourceData/MACAv2-METDATA_proj/linearshift_proj/rhmin_2070_2099")

#Then calculate VPD,EDD,GDD
for (k in 1:2){
  if (k==1){
    yr<-"2020_2049"
  }
  if (k==2){
    yr<-"2070_2099"
  }
  RHmean<-(get(paste("shiftedMetrhmax_",yr,sep=""))+get(paste("shiftedMetrhmin_",yr,sep="")))/2
  Tmean<-(get(paste("shiftedMettmax_",yr,sep=""))+get(paste("shiftedMettmin_",yr,sep="")))/2
  e_s<-6.112*exp((17.269*Tmean)/(Tmean-237.3)) #saturated vapor pressure (hPa)
  VPD<-e_s*(1-RHmean/100)
  GDD<-matrix(NA,nrow=countynum,ncol=totdays)
  EDD<-matrix(NA,nrow=countynum,ncol=totdays)
  for (i in 1: countynum){
    for (j in 1: totdays){
      GDD[i,j]<-GDDEDD(get(paste("shiftedMettmax_",yr,sep=""))[i,j],get(paste("shiftedMettmin_",yr,sep=""))[i,j])[1]
      EDD[i,j]<-GDDEDD(get(paste("shiftedMettmax_",yr,sep=""))[i,j],get(paste("shiftedMettmin_",yr,sep=""))[i,j])[2]
    }
  }
  
  
  Tthres<-10
  GS_start<-matrix(NA,nrow=countynum,ncol=30)
  GS_end<-matrix(NA,nrow=countynum,ncol=30)
  m1=1
  m2=1
  for (i in 1:countynum){
    m1=1
    m2=1
    for (j in 1:30){
      k=365
      if (j%%4==0){ 
        k=366
      }
      m2<-m1+k-1
      Tmeantoget<-Tmean[i,c(m1:m2)]
      TMA<-ma(Tmeantoget,21)
      GS_start[i,j]<-which(TMA >= Tthres)[1]
      GS_end[i,j]<-GS_start[i,j]+184
      m1=m2+1
    }
  }
  
  #create the climate projection dataframe
  Tmax_GS<-rep(NA,countynum*30)
  Tmin_GS<-rep(NA,countynum*30)
  GDD_GS<-rep(NA,countynum*30)
  EDD_GS<-rep(NA,countynum*30)
  VPD_GS<-rep(NA,countynum*30)
  Pr_GS<-rep(NA,countynum*30)
  GS_length<-rep(NA,countynum*30)
  
  for (i in 1:countynum){
    m1=1
    m2=1
    for (j in 1:30){
      k=365
      if (j%%4==0){ #when exist Feb29
        k=366
      }
      m2<-m1+k-1
      GS_length[(i-1)*30+j]<-GS_end[i,j]-GS_start[i,j]
      index<-c(m1:m2)
      if (!is.na(GS_start[i,j])){
        GSindex<-index[GS_start[i,j]:GS_end[i,j]]
        Tmax_GS[(i-1)*30+j]<-mean(get(paste("shiftedMettmax_",yr,sep=""))[i,GSindex])
        Tmin_GS[(i-1)*30+j]<-mean(get(paste("shiftedMettmin_",yr,sep=""))[i,GSindex])
        GDD_GS[(i-1)*30+j]<-sum(GDD[i,GSindex])
        EDD_GS[(i-1)*30+j]<-sum(EDD[i,GSindex])
        VPD_GS[(i-1)*30+j]<-sum(VPD[i,GSindex])
        Pr_GS[(i-1)*30+j]<-sum(get(paste("shiftedMetpr_",yr,sep=""))[i,GSindex])
      }
      m1=m2+1
    }
  }
  
  ANSI<-load("ANSI")
  CountyANSI<-load("CountyANSI")
  StateANSI<-load("StateANSI")
  if (k==1){
    Data_proj_linearshifted_2020_2049<-data.frame(StateANSI=rep(StateANSI,each=30),countyANSI=rep(CountyANSI,each=30),fips=rep(ANSI,each=30),
                                                  year=rep(c(1:30),countynum),Tmax_GS=Tmax_GS,Tmin_GS=Tmin_GS,GDD_GS=GDD_GS,GS_length=GS_length,
                                                  EDD_GS=EDD_GS,VPD_GS=VPD_GS,Pr_GS=Pr_GS)
    save(Data_proj_linearshifted_2020_2049,file="Metdata/macaprojdataframe/Data_linearshifted_2020_2049")
  }
  if (k==2){
    Data_proj_linearshifted_2070_2099<-data.frame(StateANSI=rep(StateANSI,each=30),countyANSI=rep(CountyANSI,each=30),fips=rep(ANSI,each=30),
                                                  year=rep(c(1:30),countynum),Tmax_GS=Tmax_GS,Tmin_GS=Tmin_GS,GDD_GS=GDD_GS,GS_length=GS_length,
                                                  EDD_GS=EDD_GS,VPD_GS=VPD_GS,Pr_GS=Pr_GS)
    save(Data_proj_linearshifted_2070_2099,file="Metdata/macaprojdataframe/Data_linearshifted_2070_2099")
  }
}


