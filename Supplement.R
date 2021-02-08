# This script generates the supplement plots
rm(list = ls())
graphics.off()
library(usmap)
library(MASS)
library(fields)
load("Metdata/valid_samples")
hindyear<-32
projyear<-94
parasamplenum<-length(Para)
Para<-matrix(unlist(Para),nrow = parasamplenum,ncol=13,byrow = TRUE)
load("Metdata/Metdataframe/Data_Metobs")
Data$GDD_sqr<-Data$GDD_GS^2
Data$EDD_sqr<-Data$EDD_GS^2
Data$Tmax_sqr<-Data$Tmax_GS^2
Data$Tmin_sqr<-Data$Tmin_GS^2
Data$Pr_sqr<-Data$Pr_GS^2
Data$VPD_sqr<-Data$VPD_GS^2
model<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+Tmin_GS+Tmin_sqr+
            Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data)
bestestimate<-summary(model)$coefficient[ ,1]
bestestimatestd<-summary(model)$coefficient[ ,2]
fullmodel<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmax_GS+Tmax_sqr+Tmin_GS+Tmin_sqr+
                Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data)
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
        "2070-2099 EDD projection of MACA-METDATA"),bty="n",cex=1.5)
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
                                                               "2070-2099 temperature projection of MACA-METDATA"),bty="n",cex=1.5)
dev.off()

# (3).variable distribution plots
varnames<-c("Intercept~(bu/acre)","GDD~(bu/acre/de*gree~day)","GDD^2~(bu/acre/de*gree~day^2)",
            "EDD~(bu/acre/de*gree~day)","EDD^2~(bu/acre/de*greee~day^2)","Tmax~(bu/acre/de*gree)",
            "Tmax^2~(bu/acre/de*gree^2)","Tmin~(bu/acre/de*gree)","Tmin^2~(bu/acre/de*gree^2)",
            "Pr~(bu/acre/mm)","Pr^2~(bu/acre/mm^2)","VPD~(bu/acre/hPa)","VPD^2~(bu/acre/hPa^2)")
panel<-c("A","B","C","D","E","F","G","H","I","J","K","L","M")
png(paste("Plots/supplement/para_distribution.png"), width = 1200, height = 1200)
m <- rbind(c(1:4),c(5:8),c(9:12),c(13,14,14,14))
#m <- rbind(c(1:5),c(6:10),c(11,12,13,14,14))
layout(m)
par(mar=c(6,6,2,2))
for (i in 1:13){
  vect1<-Para[ ,i]
  if ((i<5) || (i>7)){
    if (i<5){
      vect2<-dnorm(seq(min(vect1),max(vect1),len=512),mean=bestestimate[i],sd=bestestimatestd[i])
    }
    if (i>7){
      vect2<-dnorm(seq(min(vect1),max(vect1),len=512),mean=bestestimate[i-3],sd=bestestimatestd[i-3])
    }
    plot(0,0,xlim=c(min(vect1),max(vect1)),ylim=c(0,max(max(density(vect1)$y),max(vect2))),
         xlab=parse(text=varnames[i]),ylab="",type = "n",cex.axis=1.7,cex.lab=1.8)
    lines(density(vect1),xlim=c(min(vect1,0),max(vect1,0)),
          type="l",lwd=2,xlab=variable.names(model)[i],ylab="",col="red",main="",add=TRUE)
    lines(seq(min(vect1),max(vect1),len=512),vect2,xlim=c(min(vect1,0),max(vect1,0)),
          type="l",lwd=2,xlab=variable.names(model)[i],ylab="",col="black",main="",add=TRUE)
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
vec<-Para[ ,linearindx]
variablenames<-c("Intercept","GDD","EDD","Tmax","Tmin","Pr","VPD")
units<-c("(bu/acre)","(bu/acre/degree day)","(bu/acre/degree day)",
         "(bu/acre/degree)","(bu/acre/degree)","(bu/acre/mm)","(bu/acre/hPa)")
distance_name<-c(0,17,32,46.5,61,77)
distance<-c(0,11,27,44,58,74)
height<-c(0,-13.5,-13.5,-12,-12,-12,-12)
png(paste("Plots/supplement/scatter.png"), width = 1000, height = 618)
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
for (j in 1:6){
  mtext(variablenames[j], outer=FALSE, side=4, at=23, cex=1.4, las=1, line=-99.5+distance_name[j])
  mtext(units[j], outer=FALSE, side=4, at=21.5, cex=1.2, las=1, line=-99.5+distance[j])
}
for (i in 2:7){
  mtext(variablenames[i], outer=FALSE, side=4, at=3.8*(7-i), cex=1.4, las=1, line=-8)
  mtext(units[i], outer=FALSE, side=4, at=3.8*(7-i)-1.5, cex=1.2, las=1, line=height[i])
}
dev.off()

png(paste("Plots/supplement/scatterGDDEDD.png"), width = 1000, height = 618)
par(mar=c(3,3,2,3))
vec_density <- kde2d(vec[ ,2],vec[ ,3], n=2000)
image.plot(vec_density,col = hm_col_scale,xlab = variablenames[2],ylab=variablenames[3],cex.lab=1.4,cex.axis=1.4)
text(bestestimate_full[linearindx[2]],bestestimate_full[linearindx[3]]+0.06,"The best estimate",col="black",cex=1.4)
points(bestestimate_full[linearindx[2]],bestestimate_full[linearindx[3]],col="black",pch=20,cex=2)
mtext("Probability density",side=3,at=0.575,cex=1.2)
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
png(paste("Plots/supplement/convergence.png"), width = 1000, height = 618)
par(mar=c(5,5,2,5))
plot(c(2:parasamplenum),mean[2:parasamplenum],type="l",lwd=2,log="x",
     cex.lab=2,cex.axis=2,xlab="Sample size",ylab="Mean (bu/acre)")
par(new = T)
plot(c(2:parasamplenum),std[2:parasamplenum],type="l",lty=2,lwd=2,log="x",
     axes=F, xlab=NA, ylab=NA)
axis(side = 4,cex.axis=2)
mtext(side = 4, line = 3, 'Standard deviation (bu/acre)',cex = 2)
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
meanyield_anomaly<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,weight=Data$area)
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
png("Plots/supplement/time_series_full_uncertainty.png", width = 1000, height = 618)
#time series plot 
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981+3.3,2099-3.3),ylim = c(min(c(annualprojmin)),max(c(annualprojmax))),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(hind_fit+step,rev(hind_fit-step)),col="darkseagreen1",border=NA) #hindcast para
lines(c(1981:2012),hind_fit,col="forestgreen",lwd=2.5) # model best estimate
polygon(c(2012:2099,2099:2012),c(annualprojmax[7:94],rev(annualprojmin[7:94])),col="lightblue1",border=NA)#para + clim
for (i in 1:(climnum-1)){ #only climate
  lines(c(2012:2099),proj_fit[c(7:94),i],col="blue")
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=1.6)
legend(1982,-120,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield anomaly observation","Best model's hindcasts","Best model's projections"),bty="n",cex=2)
legend(1983.5,-250, fill=c("darkseagreen1","lightblue1"),
       legend = c("Hindcast window", "Yield projections uncertainty range (climate + parameter)"),bty="n",cex=2)
dev.off()

# (8) residual plot
hind<-predict(model,Data)
Data$residual<-Data$yield_anomaly-hind
Data$residual[which(Data$residual>=50)]=-50
Data$residual[which(Data$residual<=-50)]=-50
png(paste("Plots/supplement/residual.png"), width = 1000, height = 618)
i=3
index<-which(Data$year==levels(Data$year)[i])
a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=Data[index,c(3,15)], values = "residual") + labs(title = paste("Yield residual of each county in ",i+1980,sep=""))+ 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(-50,50),name="bu/acre")+
  theme(plot.title = element_text(size=20), legend.text = element_text(size=15),legend.key.size = unit(0.4, "in"),legend.title = element_text(size=15))
plot(a)
dev.off()
