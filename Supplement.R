# This script generates the supplement plots
rm(list = ls())
graphics.off()
library(usmap)
library(MASS)
library(fields)
library(ggplot2)
load("Metdata/valid_samples")
hindyear<-40
projyear<-94
parasamplenum<-length(Para)
Para<-matrix(unlist(Para),nrow = parasamplenum,ncol=13,byrow = TRUE)
load("Metdata/Metdataframe/Data_Metobs")

#Best model
model<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmax_GS+Tmin_GS+Tmin_sqr+
            +VPD_GS,weights = area,data=Data) 
bestestimate<-summary(model)$coefficient[ ,1]
bestestimatestd<-summary(model)$coefficient[ ,2]
fullmodel<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmax_GS+Tmax_sqr+Tmin_GS+Tmin_sqr+
                Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,weights = area,data=Data)
bestestimate_full<-summary(fullmodel)$coefficient[ ,1]
bestestimatestd_full<-summary(fullmodel)$coefficient[ ,2]

# (1,2). shifted/downscaled climate histogram comparison
load("Metdata/macaprojdataframe/Data_linearshifted_2020_2049")
load("Metdata/macaprojdataframe/Data_linearshifted_2070_2099")
load("Metdata/macaprojdataframe/Data_MIROC5")
EDD_proj_2070_2099<-Data_proj$EDD_GS[which((Data_proj$year>64)&Data_proj$year<95)]

png(paste("Plots/supplement/EDD_comparison.png"), width = 1000, height = 618)
par(mar=c(4.5,5.1,1.6,4.1))
plot(0,0,xlim=c(0,max(EDD_proj_2070_2099)),ylim = c(-500,4500),
     xlab = "Extreme Degree Day (Degree days)",ylab="Frequency",type = "n",yaxt="n",cex.axis=2,cex.lab=2)
axis(2,at=seq(0,4500,by=1000),cex.axis=2) #y axis labels
hist(Data_proj_linearshifted_2070_2099$EDD_GS,col=rgb(0,0,1,0.3),add=TRUE,breaks = seq(0,1040,20))
hist(EDD_proj_2070_2099,col=rgb(1,0,0,0.3),add=TRUE,breaks = seq(0,1040,20))
boxplot(Data_proj_linearshifted_2070_2099$EDD_GS, horizontal = TRUE,na.rm=TRUE, xaxt="n",col=rgb(0,0,1,0.3),
        frame=F, pch=20, main="", add=TRUE,at=-250,boxwex=200,cex=0.5)
boxplot(EDD_proj_2070_2099, horizontal = TRUE,na.rm=TRUE, xaxt="n",col=rgb(1,0,0,0.3),
        frame=F, pch=20, main="", add=TRUE,at=-450,boxwex=200,cex=0.5)
legend(450,4000,fill=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3)),legend=c("2070-2099 EDD projection of linear shifted climate",
        "2070-2099 EDD projection of MACAv2-METDATA"),bty="n",cex=1.5)
dev.off()

Tavg_proj_2070_2099<-(Data_proj$Tmax_GS[which((Data_proj$year>64)&Data_proj$year<95)]+
                        Data_proj$Tmin_GS[which((Data_proj$year>64)&Data_proj$year<95)])/2
png(paste("Plots/supplement/Tavg_comparison.png"), width = 1000, height = 618)
par(mar=c(4.5,5.1,1.6,4.1))
plot(0,0,xlim=c(min(Tavg_proj_2070_2099),max(Tavg_proj_2070_2099)),ylim = c(-1000,7000),
     xlab = "Average temperature (Degree Celcius)",ylab="Frequency",type = "n",yaxt="n",cex.axis=2,cex.lab=2)
axis(2,at=seq(0,7000,by=1000),cex.axis=2) #y axis labels
hist((Data_proj_linearshifted_2070_2099$Tmax_GS+Data_proj_linearshifted_2070_2099$Tmin_GS)/2,
     col=rgb(0,0,1,0.3),add=TRUE,breaks = seq(0,30,0.5))
hist(Tavg_proj_2070_2099,col=rgb(1,0,0,0.3),add=TRUE,breaks = seq(0,30,0.5))
boxplot((Data_proj_linearshifted_2070_2099$Tmax_GS+Data_proj_linearshifted_2070_2099$Tmin_GS)/2, 
        horizontal = TRUE,na.rm=TRUE, xaxt="n",col=rgb(0,0,1,0.3),frame=F, pch=20, main="", add=TRUE,at=-500,boxwex=400,cex=0.5)
boxplot(Tavg_proj_2070_2099, horizontal = TRUE,na.rm=TRUE, xaxt="n",col=rgb(1,0,0,0.3),
        frame=F, pch=20, main="", add=TRUE,at=-850,boxwex=400,cex=0.5)
legend(15,7300,fill=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3)),legend=c("2070-2099 temperature projection of linear shifted climate",
                                                               "2070-2099 temperature projection of MACAv2-METDATA"),bty="n",cex=1.5)
dev.off()

# (3).variable distribution plots
varnames<-c("Intercept~(bu/acre)","GDD~(bu/acre/de*gree~day)","GDD^2~(bu/acre/de*gree~day^2)",
            "EDD~(bu/acre/de*gree~day)","EDD^2~(bu/acre/de*greee~day^2)","Tmax~(bu/acre/de*gree)",
            "Tmax^2~(bu/acre/de*gree^2)","Tmin~(bu/acre/de*gree)","Tmin^2~(bu/acre/de*gree^2)",
            "Pr~(bu/acre/mm)","Pr^2~(bu/acre/mm^2)","VPD~(bu/acre/hPa)","VPD^2~(bu/acre/hPa^2)")
panel<-c("A","B","C","D","E","F","G","H","I","J","K","L","M")
png(paste("Plots/supplement/para_distribution_new.png"), width = 1200, height = 1200)
m <- rbind(c(1:4),c(5:8),c(9:12),c(13,14,14,14))
layout(m)
par(mar=c(6,6,2,2))
for (i in 1:13){ #Skip the 7, 10, 11, 13th parameter
  # vect1 is the full posterior parameter distribution
  # vect2 is the best model's parameter normal distribution
  vect1<-Para[ ,i]
  if ((i!=7) & (i!=10) & (i!=11) & (i!=13)){
    if (i<7){
      vect2<-dnorm(seq(min(vect1),max(vect1),len=512),mean=bestestimate[i],sd=bestestimatestd[i])
    }
    if ((i>7) & (i<10)){
      vect2<-dnorm(seq(min(vect1),max(vect1),len=512),mean=bestestimate[i-1],sd=bestestimatestd[i-1])
    }
    if (i==12){
      vect2<-dnorm(seq(min(vect1),max(vect1),len=512),mean=bestestimate[i-3],sd=bestestimatestd[i-3])
    }
    plot(0,0,xlim=c(min(vect1),max(vect1)),ylim=c(0,max(max(density(vect1)$y),max(vect2))),
         xlab=parse(text=varnames[i]),ylab="",type = "n",cex.axis=1.7,cex.lab=1.8)
    lines(density(vect1),xlim=c(min(vect1,0),max(vect1,0)),
          lwd=2,xlab=variable.names(model)[i],ylab="",col="red",main="")
    lines(seq(min(vect1),max(vect1),len=512),vect2,xlim=c(min(vect1,0),max(vect1,0)),
          lwd=2,xlab=variable.names(model)[i],ylab="",col="black",main="")
    text(min(vect1),max(max(density(vect1)$y),max(vect2)),panel[i],cex=2)
  } else{
    plot(density(vect1),xlim=c(min(vect1,0),max(vect1,0)),cex.axis=1.7,cex.lab=1.8,
         type="l",lwd=2,xlab=parse(text=varnames[i]),ylab=" ",col="red",main="")
    text(min(vect1),max(density(vect1)$y),panel[i],cex=2)
    abline(v=0,lty=2,lwd=2)
  }
}
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(0.55,0.8, lty=1,lwd=2,col="red",bty="n",cex=2.8,seg.len=2.4,
       legend=c("Precalibration samples"))
legend(0.55,0.3, lty=1,lwd=2,col="black",bty="n",cex=2.8,seg.len=2.4,
       legend=c("The best model's linear regression"))
dev.off()


# (4, 5). 2D paremeter distributions pair plot
hm_col_scale<-colorRampPalette(c("black","blue","green","orange","red"))(1000)
linearindx<-c(1,2,4,6,8,10,12)
vec<-Para[  ,linearindx]
variablenames<-c("Intercept","GDD","EDD","Tmax","Tmin","Pr","VPD")
units<-c("(bu/acre)","(bu/acre/degree day)","(bu/acre/degree day)",
         "(bu/acre/degree)","(bu/acre/degree)","(bu/acre/mm)","(bu/acre/hPa)")
distance_name<-c(0:5)*22 + c(-1,1,1,2,3,5)
distance<-c(0:5)*22+c(-1,-3,-3,-1,0,0)
height<-c(0,-13.5,-13.5,-12,-12,-12,-12)
png(paste("Plots/supplement/scatter_new.png"), width = 1000*1.5, height = 618*1.5)
par(mfrow=c(7,7),mar=c(3,3,1,1))
for (i in 1:7){
  for (j in 1:7){
    vec_density <- kde2d(vec[ ,j],vec[ ,i], n=2000)
    if (j<i){
      image(vec_density,col = hm_col_scale,xlab = "",ylab="",cex.axis=1.4)
      points(bestestimate_full[linearindx[j]],bestestimate_full[linearindx[i]],col="black",pch=20,cex=2)
    } else{
      plot(0,type='n',axes=FALSE,ann=FALSE)
    }
  }
}
#Labels at the top
for (j in 1:6){
  mtext(variablenames[j], outer=FALSE, side=4, at=18.5, cex=1.4, las=1, line=-147+distance_name[j])
  mtext(units[j], outer=FALSE, side=4, at=17.5, cex=1.2, las=1, line=-147+distance[j])
}
#Labels on the right
for (i in 2:7){
  mtext(variablenames[i], outer=FALSE, side=4, at=3.2*(7-i), cex=1.4, las=1, line=-8)
  mtext(units[i], outer=FALSE, side=4, at=3.2*(7-i)-1.1, cex=1.2, las=1, line=height[i])
}
dev.off()

png(paste("Plots/supplement/scatterGDDEDD_new.png"), width = 1000, height = 618)
par(mar=c(5,5,2,3))
vec_density <- kde2d(vec[ ,2],vec[ ,3], n=2000)
image.plot(vec_density,col = hm_col_scale,xlab = paste(variablenames[2],units[2],sep=" "),
           ylab=paste(variablenames[3],units[3],sep=" "),cex.lab=1.4,cex.axis=1.4)
text(bestestimate_full[linearindx[2]],bestestimate_full[linearindx[3]]+0.01,"The best estimate",col="black",cex=1.4)
points(bestestimate_full[linearindx[2]],bestestimate_full[linearindx[3]],col="black",pch=20,cex=2)
mtext("Probability density",side=3,at=0.55,cex=1.2)
dev.off()

# (6).divergence test
load("Metdata/proj_30yravg_2070_2099")
projyears<-30
climnum<-dim(proj_30yravg)[1]
parasamplenum<-dim(proj_30yravg)[2]
mean<-rep(NA,parasamplenum)
std<-rep(NA,parasamplenum)
for(i in 1:parasamplenum){
  mean[i]<-mean(proj_30yravg[ ,c(1:i)])
  std[i]<-sqrt(var(as.vector(proj_30yravg[ ,c(1:i)])))
}
png(paste("Plots/supplement/convergence_new.png"), width = 1000, height = 618)
par(mar=c(5,5,2,5))
plot(c(2:parasamplenum),mean[2:parasamplenum],type="l",lwd=2,log="x",
     cex.lab=2,cex.axis=2,xlab="Sample size",ylab="Mean (bushel/acre)")
par(new = T)
plot(c(2:parasamplenum),std[2:parasamplenum],type="l",lty=2,lwd=2,log="x",
     axes=F, xlab=NA, ylab=NA)
axis(side = 4,cex.axis=2)
mtext(side = 4, line = 3, 'Standard deviation (bushel/acre)',cex = 2)
legend(500,65,lty=c(1,2),lwd=c(2,2),col=c("black","black"),
       legend=c("Mean","Standard deviation"),bty="n",cex=2)
dev.off()

# (7) full uncertainty time series
#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
rowMax<-function(data){
  apply(data,1,max,na.rm=TRUE)
}
rowMin<-function(data){
  apply(data,1,min,na.rm=TRUE)
}
meanyield_anomaly<-rep(NA,hindyear)
for (i in 1:hindyear){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,w=Data$area[indx])
}
load("Metdata/hind_bestfit")
load("Metdata/proj_bestfit")
load("Metdata/step")
load("Metdata/proj_parasample")
hindyears<-length(hind_fit)
projyears<-dim(proj_parasample)[1]
climnum<-dim(proj_parasample)[2]+1
parasamplenum<-dim(proj_parasample)[3]
#Annual range of 2 uncertainty sources
annualprojmin<-rowQuantiles(matrix(proj_parasample,nrow=projyears,ncol=(climnum-1)*parasamplenum),probs=0)
annualprojmax<-rowQuantiles(matrix(proj_parasample,nrow=projyears,ncol=(climnum-1)*parasamplenum),probs=1)
#Annual range of only parametric uncertainty (linear shifted climate in two 30-year periods)
load("Metdata/proj_linearshifted_bestfit_2020_2049")
load("Metdata/proj_linearshifted_bestfit_2070_2099")
load("Metdata/proj_linearshifted_parasample_2020_2049")
load("Metdata/proj_linearshifted_parasample_2070_2099")
png("Plots/supplement/time_series_full_uncertainty_new.png", width = 1000, height = 618)
#time series plot 
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1979+3.3,2099-3.3),ylim = c(min(annualprojmin),max(annualprojmax)),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1979:2018,2018:1979),c(hind_fit+step,rev(hind_fit-step)),col="darkseagreen1",border=NA) #hindcast para
lines(c(1979:2018),hind_fit,col="forestgreen",lwd=2.5) # model best estimate
polygon(c(2018:2099,2099:2018),c(annualprojmax[13:94],rev(annualprojmin[13:94])),col="lightblue1",border=NA)#para + clim
for (i in 1:(climnum-1)){ #only climate
  lines(c(2018:2099),proj_fit[c(13:94),i],col="blue")
}
points(c(1979:2018),meanyield_anomaly,col="black",pch=20,cex=1.6)
legend(1980,-120,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield anomaly observation","Best model's hindcasts","Best model's projections"),bty="n",cex=2)
legend(1981.5,-194, fill=c("darkseagreen1","lightblue1"),
       legend = c("Hindcast window", "Yield projections uncertainty range (climate + parameter)"),bty="n",cex=2)
dev.off()

# (8) residual plot
hind<-predict(model,Data)
Data$residual<-Data$yield_anomaly-hind
Data$residual[which(Data$residual>=60)]=60
Data$residual[which(Data$residual<=-60)]=-60
png(paste("Plots/supplement/residual_new.png"), width = 1000, height = 618)
i=3
index<-which(Data$year==levels(Data$year)[i])
a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data[index,c(3,21)], values = "residual") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-60,60),name="bu/acre")+
  theme(plot.title = element_text(size=20), legend.text = element_text(size=15),legend.key.size = unit(0.4, "in"),legend.title = element_text(size=15))
plot(a)
dev.off()


# (9) predictive performances between using 32 years data and 40 years data
load("Metdata/Metdataframe/Data_Metobs")

for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i+2])
  if (i==1){
    old_indx<-indx
  }
  old_indx<-append(old_indx,indx)
}
old_Data <- Data[old_indx, ]
#Old best model
model<-lm(yield~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmin_GS+Tmin_sqr+
            VPD_GS+VPD_sqr+fips+year,data=old_Data,weights = area) 
Coef<-summary(model)$coefficients
old_Data$yield_anomaly<-rep(NA,dim(old_Data)[1])
for (i in 1:dim(old_Data)[1]){
  row_year<-which(row.names(Coef)==paste("year",old_Data$year[i],sep=""))
  if (length(row_year)==0){
    yeareffect<-0
  } else{
    yeareffect<-Coef[row_year,1]
  }
  old_Data$yield_anomaly[i]<-old_Data$yield[i]-yeareffect
}
old_Data$yield_anomaly<-old_Data$yield_anomaly-weighted.mean(old_Data$yield_anomaly,na.rm=T,w = old_Data$area)

old_fit<-rep(NA,32)
for (i in 1:32){
  indx<-which(old_Data$year==levels(old_Data$year)[i+2])
  old_fit[i]<-weighted.mean(old_Data$yield_anomaly[indx],na.rm=T,w = old_Data$area[indx])
} 

model<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmin_GS+Tmin_sqr+
            VPD_GS+VPD_sqr,data=old_Data,weights = area) 
old_proj<-predict(model,Data)
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

png(paste("Plots/supplement/comparison_new.png"), width = 1000, height = 618)
par(mar=c(4,5.1,1.6,2.1))
plot(years,hind_fit,type="l",xlab="Year",ylab="Yield anomaly (bu/acre)",
     ylim=c(-40,15),col="red",cex.axis=1.8,cex.lab=1.8,lwd=2) #best fit with 40 yrs data
points(years,meanyield_anomaly,pch=20,cex=2) #Observation
lines(c(1981:2012),old_fit,type = "l",col="blue") #best fit with 32 yrs data
points(c(1979:1980,2013:2018),old_proj_annual,col="blue",pch=1,cex=3)
legend("bottomleft",col=c("black","blue","red","blue"),lty=c(NA,NA,1,1),
       lwd=c(NA,NA,2,2),pch=c(20,1,NA,NA),bty="n",
       legend=c("Observations","Yield hindcast calculated by 32 years data",
       "Best fit with 40 years data","Best fit with 32 years data"),cex=1.8)
dev.off()