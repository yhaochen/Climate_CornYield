
#remove all data/variables and plots
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(ncdf4)
library(ggplot2)
library(usmap)
library(precintcon)
library(BMS)
library(adaptMCMC)
library(geoR) #for variogram
library(housingData) #for returning county centroid
library(gstat)
library(fields)
load("Data_prism")
Data_prism<-Data_prism[complete.cases(Data_prism), ]
Data_prism$StateANSI<-factor(Data_prism$StateANSI)
Data_prism$year<-factor(Data_prism$year)
### reimplement S&R paper's model using GDD/EDD (piecewise linear model)
#quadratic yield trend in each state 
statenum<-length(levels(Data_prism$StateANSI))
Year<-length(levels(Data_prism$year))
Yield_stateyear<-rep(NA,statenum*Year) # the quadratic yield of each state each year
for (i in 1:statenum){
  index<-which(Data_prism$StateANSI == levels(Data_prism$StateANSI)[i])
  Yield_state<-Data_prism$Yield[index]
  Year_state<-Data_prism$year[index]
  
  newdata<-data.frame(Yield=log(Yield_state),year=as.integer(Year_state),year_sqr=as.integer(Year_state)^2)
  model<-lm(Yield~year+year_sqr,data = newdata)
  Yield_predict<-predict(model,newdata) #then take average yield in each year
  for (j in 1:Year){
    index<-which(Year_state == levels(Data_prism$year)[j])
    Yield_stateyear[(i-1)*Year+j]=mean(Yield_predict[index])
  }
}

trendindex<-(as.integer(Data_prism$StateANSI)-1)*Year+as.integer(Data_prism$year)  #find the index of detrending yield based on state and year
Data_prism$logYield<-log(Data_prism$Yield)-Yield_stateyear[trendindex]   #this anomaly is: yield - the quadratic trend



Datatouse<-Data_prism[complete.cases(Data_prism), ]
plot(Datatouse$FF_anomaly,Datatouse$resfix,pch=20,xlab="Number of spring frost events lasts for at least 3 days",ylab="S&R model residuals (bush/acre)",cex.lab=1.3,cex.axis=1.3)
abline(lm(resfix~FF_anomaly, data = Datatouse),col="red",lwd=3)

Datatouse$SF_anomaly<-round(Datatouse$SF_anomaly/5)*5
Datatouse$FF_anomaly<-round(Datatouse$FF_anomaly/5)*5
par(mar=c(5,5,3,2))
boxplot(resfix~FF_anomaly, data = Datatouse,
        main="",xlab="First fall frost date anomaly (days)", ylab="Yield residuals (bush/acre)",
        col="orange",border="brown", boxwex = 0.5,cex.lab=1.5,cex.axis=1.5
)



#model_fix<-lm(logYield~GDD_anomaly+EDD_anomaly+Pr+Pr^2+fips+year,data=Data_prism)
model_fix<-lm(logYield~GDD_changingGS_anomaly+EDD_changingGS_anomaly+lastSF+firstFF+Pr_changingGS+Pr_changingGS^2+fips+year,data=Data_prism,weights=Data_prism$Area)
#model_fix<-lm(logYield~GDD_changingGS_anomaly+SF_anomaly+FF_anomaly,data=Data_prism)
hind_fix<-predict(model_fix, Data_prism)
Data_prism$hind<-exp(hind_fix+Yield_stateyear[trendindex])
Data_prism$res<-Data_prism$Yield-Data_prism$hind
#CV error test for each model 

totlen=dim(Data_prism)[1]
cvind<-c(1:totlen)


groups<-split(cvind,sample(rep(1:10,each = floor(totlen/10))))
CV_each<-rep(0,10)
for (j in 1:10) {  #10 cross-validation tests
  test_indx<-groups[[j]]
  test_data<-Data_prism[test_indx, "Yield"]
  train_indx<-cvind[-test_indx]
  train_model<-model_fix
  pred<-predict(train_model,Data_prism[test_indx, ]) #which variables to test?
  trendindex_1<-trendindex[test_indx]
  CV_each[j]<-mean((exp(pred+Yield_stateyear[trendindex_1])-test_data)^2)  #used for log model
  #CV_each[j]<-mean((pred+mean(Data_prism$Yield)-test_data)^2)
}
CV<-mean(CV_each) 
CV<-sqrt(CV)/mean(Data_prism[ ,"Yield"],na.rm=TRUE)
CV

i=3
index<-which(Data_prism$year==i)
plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
           data=Data_prism[index,c(12,25)], values = "res") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-45,45),name="bu/acre")+theme(plot.title = element_text(size=14))

index<-which(Data_prism$year==i)
plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
           data=Data_prism[index,c(12,28)], values = "resfix") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-45,45),name="bu/acre")+theme(plot.title = element_text(size=14))

