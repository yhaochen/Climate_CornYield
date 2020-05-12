
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
library(binaryLogic)
load("Data_prism")
Data_prism<-Data_prism[complete.cases(Data_prism), ]
Data_prism$StateANSI<-factor(Data_prism$StateANSI)
Data_prism$year<-factor(Data_prism$year)
### reimplement S&R paper's model using GDD/EDD (piecewise linear model)
#quadratic yield trend in each state 

statenum<-length(levels(Data_prism$StateANSI))
countynum<-length(levels(Data_prism$fips))
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


# a list of model structures and CV errors
Tstructure<-c(19,18,22,21)
Pstructure<-c(24,8,23)
Tnum<-2^length(Tstructure)
Pnum<-2^length(Pstructure)
totlen=dim(Data_prism)[1]
cvind<-c(1:totlen)
CV<-matrix(NA,nrow = Tnum,ncol=Pnum)
AIC<-matrix(NA,nrow = Tnum,ncol=Pnum)
BIC<-matrix(NA,nrow = Tnum,ncol=Pnum)

Tnames<-rep("/",Tnum)
for(i in 1:(Tnum-1)){
  Tindx<-as.binary(i,n=length(Tstructure))
  Tnames[i+1]<-paste(colnames(Data_prism)[Tstructure[Tindx]],collapse = "+")
}
Pnames<-rep("/",Pnum)
for(i in 1:(Pnum-1)){
  Pindx<-as.binary(i,n=length(Pstructure))
  Pnames[i+1]<-paste(colnames(Data_prism)[Pstructure[Pindx]],collapse = "+")
}

for (i in 1:Tnum){
  for (j in 1:Pnum){
    if (i==1 & j==1){
      model_fix<-lm(logYield~1,data=Data_prism, weights=Data_prism$Area)
    } else if (i==1 & j>1){
      model_fix<-lm(as.formula(paste("logYield ~ ",Pnames[j], sep="") ),data=Data_prism, weights=Data_prism$Area)
    } else if (j==1 & i>1){
      model_fix<-lm(as.formula(paste("logYield ~ ",Tnames[i], sep="") ),data=Data_prism, weights=Data_prism$Area)
    } else {
      model_fix<-lm(as.formula(paste("logYield ~ ",paste(Tnames[i], Pnames[j],sep="+"), sep="") ),data=Data_prism, weights=Data_prism$Area)
    }
    groups<-split(cvind,sample(rep(1:10,each = floor(totlen/10))))
    CV_each<-rep(0,10)
    for (k in 1:10) {  #10 cross-validation tests
      test_indx<-groups[[k]]
      test_data<-Data_prism[test_indx, "logYield"]
      train_indx<-cvind[-test_indx]
      train_model<-model_fix
      pred<-predict(train_model,Data_prism[test_indx, ]) #which variables to test?
      trendindex_1<-trendindex[test_indx]
      CV_each[k]<-mean( (exp(pred+Yield_stateyear[trendindex_1])-exp(test_data+Yield_stateyear[trendindex_1]))^2 )
    }
    CV[i,j]<-sqrt(mean(CV_each))/mean(Data_prism[ ,"Yield"],na.rm=TRUE)
    AIC[i,j]<-AIC(model_fix)
    BIC[i,j]<-BIC(model_fix)
  }
}

Tstructure<-c("FF","SF","EDD","GDD")
Tnames<-rep("/",Tnum)
for(i in 1:(Tnum-1)){
  Tindx<-as.binary(i,n=length(Tstructure))
  Tnames[i+1]<-paste(Tstructure[Tindx],collapse = "+")
}
Pstructure<-c("Pr_Jul","SPI","Pr")
Pnames<-rep("/",Pnum)
for(i in 1:(Pnum-1)){
  Pindx<-as.binary(i,n=length(Pstructure))
  Pnames[i+1]<-paste(Pstructure[Pindx],collapse = "+")
}

par(mar = c(5, 20, 4.5, 2))
colorlow<-"darkseagreen1"
colorhigh<-"white"
color2D.matplot(AIC, show.values = TRUE,axes = FALSE,show.legend=TRUE,xlab = "",ylab = "",vcex = 2,vcol = "black",extremes = c(colorlow,colorhigh))
axis(3, at = seq_len(ncol(CV)) - 0.5,
     labels = Pnames, tick = FALSE, cex.axis = 2)
axis(2, at = seq_len(nrow(CV)) - 0.5,
     labels = rev(Tnames), tick = FALSE, las = 1, cex.axis = 2)






#each year's res plot
dir.create("fig_res")
load("S&Rimplementation/Prism_S&Rmodel")
for (i in 1:32){
  filename<-paste("fig_res/res",i,".jpeg",sep="")
  jpeg(file = filename,width = 800,height=800)
  index<-which(Data_prism$year==i)
  a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data_prism[index,c(12,16)], values = "resfix") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-45,45),name="bu/acre")+theme(plot.title = element_text(size=14))
  plot(a)
  dev.off()
}

#each year's GS length
SF<-read.table("LastSpringFrost") #SF[i.j]: ith county jth year
FF<-read.table("FirstFallFrost")
countymeanSF<-floor(rowMeans(SF))
countymeanFF<-floor(rowMeans(FF))
countymeanGS<-countymeanFF-countymeanSF
Data_prism$GS<-Data_prism$firstFF-Data_prism$lastSF
Data_prism$GS_anomaly<-NA
for (i in 1:47498){
  countyindx<-which(levels(Data_prism$fips)==Data_prism$fips[i])
  Data_prism$GS_anomaly[i]<-Data_prism$GS[i]-countymeanGS[countyindx]
}

dir.create("fig_GS")
for (i in 1:32){
  filename<-paste("fig_GS/GS",i,".jpeg",sep="")
  jpeg(file = filename,width = 800,height=800)
  index<-which(Data_prism$year==i)
  a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
               data=Data_prism[index,c(12,26)], values = "GS_anomaly") + labs(title = paste("Growing season length anomaly of each county in ",i+1980,sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-60,60),midpoint=0,name="days")+theme(plot.title = element_text(size=14))
  plot(a)
  dev.off()
}

# A time bar
year<-c(1981:2012)
for (i in 1:32){
  filename<-paste("fig_res/year",i,".jpeg",sep="")
  jpeg(file = filename,width = 900,height=500)
  plot(0,xlim=c(1980,2013),type="n",axes=FALSE,xlab="",ylab="")
  axis(1, at = year, labels = year)
  points(x=year[i],y=-1.05,pch=20,col="red",cex=2)
  dev.off()
}






par(mfrow=c(6,6))
for(i in 1:32){
  index<-which(Data_prism$year==i)
  plot(plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
               data=Data_prism[index,c(12,16)], values = "resfix") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
    scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-45,45),name="bu/acre")+theme(plot.title = element_text(size=14)))
}




#correlation between growing season and residuals
Data_prism$GS<-Data_prism$firstFF-Data_prism$lastSF
plot(Data_prism$GS,Data_prism$resfix,pch=20)
   #1. correlation through time at each county, shown by map
Cordataframe<-data.frame(fips=levels(Data_prism$fips),cor=rep(NA,countynum))
for (i in 1:countynum){
  indx<-which(Data_prism$fips==levels(Data_prism$fips)[i])
  Resdata<-Data_prism$resfix[indx]
  GSdata<-Data_prism$GS[indx]
  Cordataframe$cor[i]<-cor(Resdata,GSdata)
}
plot(plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
                data=Cordataframe, values = "cor") + labs(title = "")+ 
       scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-1,1),name="")+theme(plot.title = element_text(size=14)))

  #2. correlation of all county data in each year, shown by time series
Corr<-rep(NA,Year)
for (i in 1: Year){
  indx<-which(Data_prism$year==levels(Data_prism$year)[i])
  Resdata<-Data_prism$resfix[indx]
  GSdata<-Data_prism$GS[indx]
  Corr[i]<-cor(Resdata,GSdata)
}
par(mar=c(5,5,3,2))
plot(c(1981:2012),Corr,pch=20,xlab="Year",ylab="Correlation coefficient",cex.lab=1.4,cex.axis=1.4)
