# Model selection by choosing the model with least cross validation error
rm(list = ls())
graphics.off()
load("Metdata/Metdataframe/Data_Metobs")

set.seed(1)
#model is the model structure in text, data is the dataframe used to calculate hindcast, obs is observation data
CV<-function(model,data,obs){ 
  totlen<-length(obs)
  totindx<-c(1:totlen)
  groups<-split(totindx,sample(c(1:10),size=totlen,replace = TRUE))
  CV_each<-rep(0,10)
  for (i in 1:10) {  #10 cross-validation tests
    test_indx<-groups[[i]]
    test_obs<-obs[test_indx]
    train_indx<-totindx[-test_indx]
    train_model<-lm(as.formula(model),weights = area,data=data[train_indx, ])
    options(warn=-1)
    pred<-predict(train_model,data[test_indx, ]) #which variables to test?
    CV_each[i]<-mean((pred-obs[test_indx])^2)
  }
  return(mean(CV_each))
}
#must include GDD,EDD
GDD<-c("GDD_GS","poly(GDD_GS,2)")
EDD<-c("EDD_GS","poly(EDD_GS,2)")
Tmax<-c(" ","Tmax_GS","poly(Tmax_GS,2)")
Tmin<-c(" ","Tmin_GS","poly(Tmin_GS,2)")
Pr<-c(" ","Pr_GS","poly(Pr_GS,2)")
VPD<-c(" ","VPD_GS","poly(VPD_GS,2)")


models<-expand.grid(GDD,EDD,Tmax,Tmin,Pr,VPD)
strunum<-dim(models)[1]
modelnames<-rep("/",strunum)
for (i in 1:strunum){
  modelnames[i]<-paste(models[i, ]$Var1,models[i, ]$Var2,models[i, ]$Var3,
                         models[i, ]$Var4,models[i, ]$Var5,models[i, ]$Var6,sep="+")
  modelnames[i]<-gsub("[+] ","",modelnames[i])
  modelnames[i]<-paste("yield~",modelnames[i],sep="")
}


cv<-rep(NA,strunum)
for (i in 1:strunum){
  cv[i]<-CV(modelnames[i],Data,Data$yield)
}

model<-lm(as.formula(paste(modelnames[which.min(cv)],"+fips+year")),weights = area,data=Data) #best model
Coef<-summary(model)$coefficients
Data$yield_anomaly<-rep(NA,dim(Data)[1])
for (i in 1:dim(Data)[1]){
  row_year<-which(row.names(Coef)==paste("year",Data$year[i],sep=""))
  if (length(row_year)==0){
    yeareffect<-0
  } else{
    yeareffect<-Coef[row_year,1]
  }
  Data$yield_anomaly[i]<-Data$yield[i]-yeareffect
}
Data$yield_anomaly<-Data$yield_anomaly-weighted.mean(Data$yield_anomaly,na.rm=T,w = Data$area)
save(Data,file = "Metdata/Metdataframe/Data_Metobs")
