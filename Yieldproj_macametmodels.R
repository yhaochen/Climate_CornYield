#This file reads the yield projections in different model structures/climate projections and makes plots

rm(list = ls())
graphics.off()
library(matrixStats)

#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)

rowMax<-function(data){
  apply(data,1,max,na.rm=TRUE)
}

rowMin<-function(data){
  apply(data,1,min,na.rm=TRUE)
}

meanyield_anomaly<-rep(NA,40)
for (i in 1:40){
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
annualprojmin<-rowQuantiles(matrix(proj_parasample,nrow=projyears,ncol=(climnum-1)*parasamplenum),probs=0.025)
annualprojmax<-rowQuantiles(matrix(proj_parasample,nrow=projyears,ncol=(climnum-1)*parasamplenum),probs=0.975)
#Annual range of only parametric uncertainty (linear shifted climate in two 30-year periods)
load("Metdata/proj_linearshifted_bestfit_2020_2049")
load("Metdata/proj_linearshifted_bestfit_2070_2099")
load("Metdata/proj_linearshifted_parasample_2020_2049")
load("Metdata/proj_linearshifted_parasample_2070_2099")

png("Plots/time_series.png", width = 1000, height = 618)
#time series plot 
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1979+3.3,2099-3.3),ylim = c(min(c(annualprojmin)),max(c(annualprojmax))),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1979:2018,2018:1979),c(hind_fit+step,rev(hind_fit-step)),col="darkseagreen1",border=NA) #hindcast para
lines(c(1979:2018),hind_fit,col="forestgreen",lwd=2.5) # model best estimate
polygon(c(2018:2099,2099:2018),c(annualprojmax[13:94],rev(annualprojmin[13:94])),col="lightblue1",border=NA)#para + clim
for (i in 1:(climnum-1)){ #only climate
  lines(c(2018:2099),proj_fit[c(13:94),i],col="blue")
}
points(c(1979:2018),meanyield_anomaly,col="black",pch=20,cex=1.6)
legend(1980,-75,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield anomaly observation","Best model's hindcasts","Best model's projections"),bty="n",cex=2)
legend(1981.5,-110, fill=c("darkseagreen1","lightblue1"),
       legend = c("Hindcast window", "95% yield projections uncertainty range (climate + parameter)"),bty="n",cex=2)
dev.off()



 #Plot of marginal yield distributions at 2020-2049/2070-2099
 for (k in 1:2){
   if (k==1){
     load("Metdata/proj_linearshifted_bestfit_2020_2049")
     load("Metdata/proj_linearshifted_parasample_2020_2049")
     selectedyears<-c(15:44)
     load("Metdata/proj_30yravg_2020_2049")
     load("Metdata/proj_fit_30yravg_2020_2049")
     load("Metdata/cumu_uncertainty_2020_2049")
   }
   if (k==2){
     load("Metdata/proj_linearshifted_bestfit_2070_2099")
     load("Metdata/proj_linearshifted_parasample_2070_2099")
     selectedyears<-c(65:94)
     load("Metdata/proj_30yravg_2070_2099")
     load("Metdata/proj_fit_30yravg_2070_2099")
     load("Metdata/cumu_uncertainty_2070_2099")
   }
   #All
   vect_para_clim<-as.vector(proj_30yravg)
   #projection: only climate
   vect_clim<-proj_fit_30yravg
   #projection: only parameter (best structure under linear shifted climate)
   vect_para<-as.vector(proj_30yravg[climnum, ])
   #projection: none (best structure's point estimates in each year under linear shifted climate)
   annualproj_none<-proj_fit_30yravg[climnum]

   #all distributions in 2070-2099/2020-2049
   if (k==1){
     boxstep<- -0.02
     yrange<-c(-0.06,0.14)
     plotname<-"distribution_2020_2049"
   }
   if (k==2){
     boxstep<- -0.007
     yrange<-c(-0.021,0.06)
     plotname<-"distribution_2070_2099"
   }
   png(paste("Plots/",plotname,".png"), width = 1000, height = 618)
   par(mar=c(4.5,5.1,3,21.1))
   plot(0,0,xlim = c(quantile(vect_para_clim,0.001),quantile(vect_para_clim,0.999)),ylim = yrange,
        xlab="Yield anomaly (bush/acre)",ylab="Density",type = "n",yaxt="n",cex.axis=2,cex.lab=2)
   axis(2,at=seq(0,yrange[2],by=0.025),cex.axis=2) #y axis labels
   points(annualproj_none,0,pch=19)
   color<-c(rgb(1,0,0,0.6),"forestgreen",rgb(0,0,1,0.6))
   text<-c("Only parameter uncertainty","Only climate uncertainty","Parameter + Climate uncertainty")
   boxcolor<-c(rgb(1,0,0),"forestgreen",rgb(0,0,1,0.5))

   #plot in the same order as the cumulative uncertainty figure result
   sorted<-sort(cumu_uncertainty)
   names<-names(sorted)
   for (i in 1:3){
     lines(density(get(paste("vect_",names[i],sep="")),na.rm=TRUE),col=color[which(cumu_uncertainty==sorted[1,i])],lwd=5)
     boxplot(get(paste("vect_",names[i],sep="")), horizontal = TRUE,na.rm=TRUE, xaxt="n",col=boxcolor[which(cumu_uncertainty==sorted[1,i])],
             frame=F, pch=20,ylim=yrange, main="", add=TRUE,at=i*boxstep,boxwex=boxstep,cex=0.5)
     #points(mean(get(paste("vect_",names[i],sep=""))),i*boxstep,pch=5,cex=1.2)
     #points(density(get(paste("vect_",names[i],sep="")))$x[which(density(get(paste("vect_",names[i],sep="")))$y==max(density(get(paste("vect_",names[i],sep="")))$y))],i*boxstep,pch=19,cex=1.2)
     mtext(side=4, at=i*boxstep, line=1,text[which(cumu_uncertainty==sorted[1,i])],las=1,col=color[which(cumu_uncertainty==sorted[1,i])],cex=1.5)
   }
   if (k==1){
     mtext(side=3, at=quantile(vect_para_clim,0.001), line=1, las=1,"a",cex=2)
     mtext(side=3, at=quantile(vect_para_clim,0.3), line=1, las=1,"2020 - 2049",cex=2)
   }
   if (k==2){
     mtext(side=3, at=quantile(vect_para_clim,0.001), line=1, las=1,"b",cex=2)
     mtext(side=3, at=quantile(vect_para_clim,0.3), line=1, las=1,"2070 - 2099",cex=2)
   }
   legend(quantile(vect_para_clim,0.001),yrange[2],col=c("black","red","forestgreen","blue"),pch=c(19,NA,NA,NA),
          lty=c(NA,1,1,1),lwd=c(NA,2,2,2),legend=c("Point estimate","Parameter uncertainty + linear shifted climate",
                                                   "Sampled climate projections best estimates","Sampled parameter + climate uncertainty"),bty="n",cex=1.7)
   dev.off()
 }

