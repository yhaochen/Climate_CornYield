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
plot(0,0,xlim = c(1981,2099),ylim = c(min(c(annualhistmin,annualprojmin)),max(c(annualhistmax,annualprojmax))),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(annualhistmax,rev(annualhistmin)),col="darkseagreen1",border=NA) #hindcast para
lines(c(1981:2012),hind_fit,col="forestgreen",lwd=2.5) # model best estimate
polygon(c(2012:2099,2099:2012),c(annualprojmax[7:94],rev(annualprojmin[7:94])),col="lightblue1",border=NA)#para + clim
#polygon(c(2020:2049,2049:2020),c(annualprojmax_2020_2049,rev(annualprojmin_2020_2049)),col="skyblue",border=NA) # para
#polygon(c(2070:2099,2099:2070),c(annualprojmax_2070_2099,rev(annualprojmin_2070_2099)),col="skyblue",border=NA) # para
for (i in 1:(climnum-1)){ #only climate
  lines(c(2012:2099),proj_fit[c(7:94),i],col="blue")
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=1.6)
legend(1978.5,-56.5,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield anomaly observation","Best model's hindcasts","Best model's projections"),bty="n",cex=2)
legend(1980,-79, fill=c("darkseagreen1","lightblue1"),
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
     yrange<-c(-0.03,0.1)
     plotname<-"distribution_2070_2099"
   }
   png(paste("Plots/",plotname,".png"), width = 1000, height = 618)
   par(mar=c(4.5,5.1,1.6,21.1))
   plot(0,0,xlim = c(min(vect_para_clim,na.rm = TRUE),max(vect_para_clim,na.rm = TRUE)),ylim = yrange,
        xlab="Yield anomaly (bush/acre)",ylab="Density",type = "n",yaxt="n",cex.axis=2,cex.lab=2)
   axis(2,at=seq(0,0.2,by=0.025),cex.axis=2) #y axis labels
   points(annualproj_none,0,pch=8)
   color<-c(rgb(1,0,0,0.6),"forestgreen",rgb(0,0,1,0.6))
   text<-c("Parameter","Climate","Parameter + Climate")
   boxcolor<-c(rgb(1,0,0),"forestgreen",rgb(0,0,1,0.5))

   #plot in the same order as the 7-layer cumulative uncertainty figure result
   sorted<-sort(cumu_uncertainty)
   names<-names(sorted)
   for (i in 1:3){
     lines(density(get(paste("vect_",names[i],sep="")),na.rm=TRUE),col=color[which(cumu_uncertainty==sorted[1,i])],lwd=5)
     if (i==3) { #thicker black line for the distribution that includes 3 uncertainty sources
       lines(density(get(paste("vect_",names[i],sep="")),na.rm=TRUE),col=color[which(cumu_uncertainty==sorted[1,i])],lwd=5)
     }
     boxplot(get(paste("vect_",names[i],sep="")), horizontal = TRUE,na.rm=TRUE, xaxt="n",col=boxcolor[which(cumu_uncertainty==sorted[1,i])],
             frame=F, pch=20,ylim=yrange, main="", add=TRUE,at=i*boxstep,boxwex=boxstep,cex=0.5)
     points(mean(get(paste("vect_",names[i],sep=""))),i*boxstep,pch=5,cex=1.2)
     points(density(get(paste("vect_",names[i],sep="")))$x[which(density(get(paste("vect_",names[i],sep="")))$y==max(density(get(paste("vect_",names[i],sep="")))$y))],i*boxstep,pch=19,cex=1.2)
     mtext(side=4, at=i*boxstep, line=1,text[which(cumu_uncertainty==sorted[1,i])],las=1,col=color[which(cumu_uncertainty==sorted[1,i])],cex=2.4)
   }

   legend(min(vect_para_clim),yrange[2],col=c("black","red","forestgreen","blue"),pch=c(8,NA,NA,NA),
          lty=c(NA,1,1,1),lwd=c(NA,2,2,2),legend=c("Point estimate","Parameter","Climate","Parameter + climate"),bty="n",cex=2)
   dev.off()
 }



# # #2D paremeter distributions
# parasamplenum<-10000
# i<-16
# j<-4
# model<-lm(as.formula(paste("yield_anomaly ~ ",paste(Tnames[i], Pnames[j],sep="+"), sep="") ),data=Data, weights=Data$area)
# MCpara<-parasample(model,parasamplenum,2)
# MCpara_keep<-MCpara
# hind_parasample<-matrix(NA,nrow=32,ncol=10000)
# variablenames<-variable.names(model)
# variablenum<-length(variablenames)


# par(mfrow=c(7,7),mar=c(1,1,1,1))
# for (i in 1:7){
#   for (j in 1:7){
#     if (j<i){
#       d<-data.frame(x=MCpara[ ,j],y=MCpara[ ,i])
#       kd <- ks::kde(d, compute.cont=TRUE)
#       contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
#                                           z=estimate, levels=cont["5%"])[[1]])
#       contour_95 <- data.frame(contour_95)
#       
#       plot(d,pch=20,cex=0.5,col="red")
#       points(MCpara_keep[ ,j],MCpara_keep[ ,i],pch=20,cex=0.5,col="grey")
#       lines(contour_95$x,contour_95$y,lwd=2,col="blue")
#     }else {
#       plot(0,type='n',axes=FALSE,ann=FALSE)
#     }
#   }
# }

