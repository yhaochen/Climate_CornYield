#This file reads the yield projections in different model structures/climate projections and makes plots

rm(list = ls())
graphics.off()

colMax <- function(data) {
  apply(data, 2, max,na.rm=TRUE)
}
colMin <- function(data) {
  apply(data, 2, min,na.rm=TRUE)
}
rowMax <- function(data) {
  apply(data, 1, max,na.rm=TRUE)
}
rowMin <- function(data) {
  apply(data, 1, min,na.rm=TRUE)
}


#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)

meanyield_anomaly<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,weight=Data$area)
}

load("Metdata/hind_bestfit")
load("Metdata/proj_bestfit")
load("Metdata/hind_parasample")
load("Metdata/proj_parasample")
hindyears<-length(hind_fit)
projyears<-dim(proj_parasample)[1]
climnum<-dim(proj_parasample)[2]+1
parasamplenum<-dim(proj_parasample)[3]

#Annual range of 2 uncertainty sources
annualhistmin<-rowMin(hind_parasample)
annualhistmax<-rowMax(hind_parasample)
annualprojmin<-rowMin(proj_parasample)
annualprojmax<-rowMax(proj_parasample)
#Annual range of only parametric uncertainty (linear shifted climate in two 30-year periods)
load("Metdata/proj_linearshifted_bestfit_2020_2049")
load("Metdata/proj_linearshifted_bestfit_2070_2099")
load("Metdata/proj_linearshifted_parasample_2020_2049")
load("Metdata/proj_linearshifted_parasample_2070_2099")
annualprojmax_2020_2049<-rowMax(proj_linearshifted_parasample_2020_2049)
annualprojmax_2070_2099<-rowMax(proj_linearshifted_parasample_2070_2099)
annualprojmin_2020_2049<-rowMin(proj_linearshifted_parasample_2020_2049)
annualprojmin_2070_2099<-rowMin(proj_linearshifted_parasample_2070_2099)

png("Plots/time_series.png", width = 1000, height = 618)
#time series plot 
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981+3.3,2099-3.3),ylim = c(min(c(annualhistmin,annualprojmin)),max(c(annualhistmax,annualprojmax))),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(annualhistmax,rev(annualhistmin)),col="darkseagreen1",border=NA) #hindcast para
lines(c(1981:2012),hind_fit,col="forestgreen",lwd=2.5) # model best estimate
polygon(c(2012:2099,2099:2012),c(annualprojmax[7:94],rev(annualprojmin[7:94])),col="lightblue1",border=NA)#para + clim
#polygon(c(2020:2049,2049:2020),c(annualprojmax_2020_2049,rev(annualprojmin_2020_2049)),col="skyblue",border=NA) # para
#polygon(c(2070:2099,2099:2070),c(annualprojmax_2070_2099,rev(annualprojmin_2070_2099)),col="skyblue",border=NA) # para
for (i in 1:(climnum-1)){ #only climate
  lines(c(2012:2099),proj_fit[c(7:94),i],col="blue")
  Sys.sleep(1)
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=1.6)
legend(1982,80,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield anomaly observation","Best model's hindcasts","Best model's projections"),bty="n",cex=2)
legend(1983.5,50.8, fill=c("darkseagreen1","lightblue1"),
       legend = c("Hindcast window", "Yield projections range (climate + parameter)"),bty="n",cex=2)
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
     yrange<-c(-0.06,0.23)
     plotname<-"distribution_2020_2049"
   }
   if (k==2){
     boxstep<- -0.01
     yrange<-c(-0.03,0.12)
     plotname<-"distribution_2070_2099"
   }
   png(paste("Plots/",plotname,".png"), width = 1000, height = 618)
   par(mar=c(4.5,5.1,1.6,21.1))
   plot(0,0,xlim = c(min(vect_para_clim,na.rm = TRUE),max(vect_para_clim,na.rm = TRUE)),ylim = yrange,
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
     if (i==3) { #thicker black line for the distribution that includes 3 uncertainty sources
       lines(density(get(paste("vect_",names[i],sep="")),na.rm=TRUE),col=color[which(cumu_uncertainty==sorted[1,i])],lwd=5)
     }
     boxplot(get(paste("vect_",names[i],sep="")), horizontal = TRUE,na.rm=TRUE, xaxt="n",col=boxcolor[which(cumu_uncertainty==sorted[1,i])],
             frame=F, pch=20,ylim=yrange, main="", add=TRUE,at=i*boxstep,boxwex=boxstep,cex=0.5)
     #points(mean(get(paste("vect_",names[i],sep=""))),i*boxstep,pch=5,cex=1.2)
     #points(density(get(paste("vect_",names[i],sep="")))$x[which(density(get(paste("vect_",names[i],sep="")))$y==max(density(get(paste("vect_",names[i],sep="")))$y))],i*boxstep,pch=19,cex=1.2)
     mtext(side=4, at=i*boxstep, line=1,text[which(cumu_uncertainty==sorted[1,i])],las=1,col=color[which(cumu_uncertainty==sorted[1,i])],cex=1.5)
   }
   if (k==1){
     mtext(side=4, at=yrange[2], line=1, las=1,"(a). 2020 - 2049",cex=2)
   }
   if (k==2){
     mtext(side=4, at=yrange[2], line=1, las=1,"(b). 2070 - 2099",cex=2)
   }
   legend(min(vect_para_clim),yrange[2],col=c("black","red","forestgreen","blue"),pch=c(19,NA,NA,NA),
          lty=c(NA,1,1,1),lwd=c(NA,2,2,2),legend=c("Point estimate","Parameter uncertainty + linear shifted climate",
                                                   "Sampled climate projections best estimates","Sampled parameter + climate uncertainty"),bty="n",cex=1.5)
   dev.off()
 }


#divergence test
load("Metdata/proj_parasample")
projyears<-dim(proj_parasample)[1]
climnum<-dim(proj_parasample)[2]
parasamplenum<-dim(proj_parasample)[3]
std<-rep(NA,parasamplenum)
for(i in 1:parasamplenum){
  std[i]<-sqrt(var(proj_parasample[ , ,c(1:i)]))
}
plot(c(2:parasamplenum),std[2:parasamplenum],type="l",lwd=2,
     cex.lab=1.5,cex.axis=1.5,xlab="sample size",ylab="standard error (bush/acre)")

#variable distribution plots

rm(list = ls())
graphics.off()
library(lhs)

load("Metdata/kept_sample_num")
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
parasamplenum<-5000000
stepwidth<-3
set.seed(666)
Data$GDD_sqr<-Data$GDD_GS^2
Data$EDD_sqr<-Data$EDD_GS^2
Data$Tmax_sqr<-Data$Tmax_GS^2
Data$Tmin_sqr<-Data$Tmin_GS^2
Data$Pr_sqr<-Data$Pr_GS^2
Data$VPD_sqr<-Data$VPD_GS^2
model<-lm(lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+EDD_sqr+Tmax_GS+Tmax_sqr+Tmin_GS+Tmin_sqr+
               Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data))
bestestimate<-summary(model)$coefficient[ ,1]
bestestimatestd<-summary(model)$coefficient[ ,2]
variablenames<-variable.names(model)
variablenum<-length(variablenames)
#EDD terms best estimates from the best model are outside the 3-sigma range of the full model
extrawidth<-c(0,0,0,0.15,0.0005,0,0,0,0,0,0,0,0)
#Latin hypercube sampling in a range for all parameters defined above
MCpara<-randomLHS(parasamplenum,variablenum) #range [0,1]
for (m in 1:variablenum){
  MCpara[ ,m]<-MCpara[ ,m]*2-1 #map to [-1,1]
  MCpara[ ,m]<-MCpara[ ,m]*(bestestimatestd[m]*stepwidth) + bestestimate[m]
}
bestmodel<-lm(yield_anomaly~GDD_GS+GDD_sqr+EDD_GS+Tmin_GS+Tmin_sqr+
                Pr_GS+Pr_sqr+VPD_GS+VPD_sqr,data=Data) #best model
bestestimate_1<-summary(bestmodel)$coefficient[ ,1]
bestestimatestd_1<-summary(bestmodel)$coefficient[ ,2]
m <- rbind(c(1:5),c(6:10), c(11,12,13,14,14))
layout(m)
par(mar=c(5,5,1,1))
for (i in 1:13){
  vect1<-MCpara[validsample,i]
  if ((i<5) || (i>7)){
    if (i<5){
      vect2<-rnorm(length(validsample),mean=bestestimate_1[i],sd=bestestimatestd_1[i])
    }
    if (i>7){
      vect2<-rnorm(length(validsample),mean=bestestimate_1[i-3],sd=bestestimatestd_1[i-3])
    }
    plot(density(vect2),xlim=c(min(vect1,vect2),max(vect1,vect2)),cex.axis=1.4,cex.lab=1.4,
         type="l",lwd=2,xlab=variable.names(model)[i],ylab="density",main="")
    lines(density(vect1),col="red",lwd=2)
    abline(v=bestestimate[i],col="red",lty=2)
    abline(v=mean(vect2),lty=2)
  } else{
    plot(density(vect1),xlim=c(min(vect1,0),max(vect1,0)),cex.axis=1.4,cex.lab=1.4,
         type="l",lwd=2,xlab=variable.names(model)[i],ylab="density",col="red",main="")
    abline(v=bestestimate[i],col="red",lty=2)
    abline(v=0,lty=2)
  }
}
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(0.5,1.5, lty=c(1,1,2,2),lwd=c(2,2,2,2),col=c("red","black","red","black"),bty="n",cex=1.7,
       legend=c("Precalibration samples","Best model's distribution","Full model's best estimate","Best model's best estimate"))

