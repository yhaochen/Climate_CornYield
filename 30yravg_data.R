#This script calculates 30-year mean yield for two time windows (2020-2049 and 2070-2099) under different uncertainty samples

#ser function is similar to seq(from,by,length.out) except "from" can be a vector
ser <- function(data,by,times){ 
  totlength<-length(data)*times
  output<-rep(NA,totlength)
  for (i in 1:times){
    cols<-c( (length(data)*(i-1)+1) : (length(data)*i) )
    output[cols]<-data+by*(i-1)
  }
  output
}


#read hindcast/projection data and calculate 30-year mean for two time windows (2020-2049 and 2070-2099)
load("Metdata/hind_bestfit")
load("Metdata/proj_bestfit")
load("Metdata/hind_parasample")
load("Metdata/proj_parasample")
totalyears<-94
for (k in 1:2){
  if (k==1){
    selectedyears<-c(15:44) # 2020-2049
    load("Metdata/proj_linearshifted_bestfit_2020_2049")
    load("Metdata/proj_linearshifted_parasample_2020_2049")
    
  }
  if (k==2){
    selectedyears<-c(65:94) # 2070-2099
    load("Metdata/proj_linearshifted_bestfit_2070_2099")
    load("Metdata/proj_linearshifted_parasample_2070_2099")
  }
  parasamplenum<-dim(proj_parasample)[3]
  strunum<-dim(proj_parasample)[1]
  climnum<-dim(proj_parasample)[2]/totalyears+1
  selectedcols<-ser(selectedyears,totalyears,climnum-1)
  
  # set up the 3_D matrix of para*stru*clim 1000*63*19 *30 samples
  proj_all<-array(NA,c(strunum,climnum*30,parasamplenum))
  proj_all[ ,c(1:(climnum*30-30)), ]<-proj_parasample[ ,selectedcols, ]
  proj_all[ ,c((climnum*30-29):(climnum*30)), ]<-proj_linearshifted_parasample
  proj_30yravg<-array(NA,c(strunum,climnum,parasamplenum)) #30-yr average for each combination
  for (i in 1:(climnum)){
    colindx<-c(1:30)+(i-1)*30
    proj_30yravg[ ,i, ]<-apply(proj_all[ ,colindx, ],c(1,3),mean)
  }
  
  proj_fit_all<-matrix(NA,nrow=strunum+1,ncol=climnum*30)
  proj_fit_all[ ,c(1:(climnum*30-30))]<-proj_fit[ ,selectedcols]
  proj_fit_all[ ,c((climnum*30-29):(climnum*30))]<-proj_linearshifted_fit
  proj_fit_30yravg<-array(NA,c(strunum+1,climnum))
  for (i in 1:climnum){
    colindx<-c(1:30)+(i-1)*30
    proj_fit_30yravg[ ,i]<-rowMeans(proj_fit_all[ ,colindx])
  }
  if (k==1){
    save(proj_fit_30yravg,file="Metdata/proj_fit_30yravg_2020_2049")
    save(proj_30yravg,file="Metdata/proj_30yravg_2020_2049")
  }
  if (k==2){
    save(proj_fit_30yravg,file="Metdata/proj_fit_30yravg_2070_2099")
    save(proj_30yravg,file="Metdata/proj_30yravg_2070_2099")
  }
}