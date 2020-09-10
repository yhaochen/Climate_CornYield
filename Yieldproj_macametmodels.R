#This file reads the yield projections in different model structures/climate projections and makes plots

rm(list = ls())
graphics.off()

colMax <- function(data) {
  apply(data, 2, max)
}
colMin <- function(data) {
  apply(data, 2, min)
}
rowMax <- function(data) {
  apply(data, 1, max)
}
rowMin <- function(data) {
  apply(data, 1, min)
}

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


#Metdata observation
load("Metdata/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
Data$year=Data$year+1978
Data$year<-factor(Data$year)

meanyield_anomaly<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-weighted.mean(Data$yield_anomaly[indx],na.rm=T,weight=Data$area)
}

load("Metdata/hind_bestfit")
load("Metdata/proj_bestfit")
load("Metdata/hind_parasample")
load("Metdata/proj_parasample")
hindyears<-32
projyears<-94
strunum<-dim(proj_parasample)[1]
climnum<-dim(proj_parasample)[2]/projyears+1
parasamplenum<-dim(proj_parasample)[3]

#Annual range of all 3 uncertainty sources
annualhistmin<-rep(NA,hindyears)
annualhistmax<-rep(NA,hindyears)
for (i in 1:hindyears){
  annualhistmin[i]<-min(hind_parasample[ ,i, ],na.rm=TRUE)
  annualhistmax[i]<-max(hind_parasample[ ,i, ],na.rm=TRUE)
}
annualprojmin<-rep(NA,projyears)
annualprojmax<-rep(NA,projyears)
for (i in 1:projyears){
  colindx<-seq(from = i, by=projyears, length.out = climnum-1)
  annualprojmin[i]<-min(proj_parasample[ ,colindx, ],na.rm=TRUE)
  annualprojmax[i]<-max(proj_parasample[ ,colindx, ],na.rm=TRUE)
}
#Annual range of structure+climate
annualprojmin_stru_clim<-rep(NA,projyears)
annualprojmax_stru_clim<-rep(NA,projyears)
for (i in 1:projyears){
  colindx<-seq(from = i, by=projyears, length.out = climnum-1)
  annualprojmin_stru_clim[i]<-min(proj_fit[ ,colindx],na.rm=TRUE)
  annualprojmax_stru_clim[i]<-max(proj_fit[ ,colindx],na.rm=TRUE)
}

png("Plots/time_series.png", width = 1000, height = 618)
#time series plot 
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2099),ylim = c(min(c(annualhistmin,annualprojmin))-5,max(c(annualhistmax,annualprojmax))),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(annualhistmax,rev(annualhistmin)),col="darkseagreen1",border=NA) #hindcast para
lines(c(1981:2012),hind_fit[64, ],col="forestgreen",lwd=2.5) #none (64th model best estimate)
polygon(c(2012:2099,2099:2012),c(annualprojmax[7:94],rev(annualprojmin[7:94])),col="lightblue1",border=NA)#para + stru + clim
polygon(c(2012:2099,2099:2012),c(annualprojmax_stru_clim[7:94],rev(annualprojmin_stru_clim[7:94])),col="skyblue",border=NA) #stru + clim
for (i in 1:(climnum-1)){ #only climate
  indx<-c(7:94)+(i-1)*projyears
  lines(c(2012:2099),proj_fit[strunum+1,indx],col="blue")
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=1.6)
legend(1978.5,-56.5,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield observation","Best model's hindcasts","Best model's projections"),bty="n",cex=2)
legend(1980,-85, fill=c("darkseagreen1","skyblue","lightblue1"),
       legend = c("Hindcast window","Yield projections range (climate + structure)",
                  "Yield projections range (climate + structure + parameter)"),bty="n",cex=2)
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
  selectedcols<-ser(selectedyears,projyears,climnum-1)
  #All
  vect_all<-as.vector(proj_30yravg)
  #projection: structure + climate 
  vect_stru_clim<-as.vector(proj_fit_30yravg)
  #projection: parameter + climate
  vect_para_clim<-as.vector(proj_30yravg[strunum, , ])
  #projection: structure + parameter
  vect_para_stru<- as.vector(proj_30yravg[ ,climnum, ])
  #projection: only climate (best structure performance: 64th)
  vect_clim<-proj_fit_30yravg[strunum+1, ]
  #projection: only structure (best estimate of each structure under linear shifted climate)
  vect_stru<-as.vector(proj_fit_30yravg[ ,climnum])
  #projection: only parameter (best structure under linear shifted climate)
  vect_para<-as.vector(proj_30yravg[strunum,climnum, ])
  #projection: none (best structure's point estimates in each year under linear shifted climate)
  annualproj_none<-proj_fit_30yravg[strunum+1,climnum]
  
  #all distributions in 2070-2099/2020-2049
  if (k==1){
    boxstep<- -0.014
    yrange<-c(-0.1,0.2)
    plotname<-"distribution_2020_2049"
  }
  if (k==2){
    boxstep<- -0.006
    yrange<-c(-0.045,0.075)
    plotname<-"distribution_2070_2099"
  }
  png(paste("Plots/",plotname,".png"), width = 1000, height = 618)
  par(mar=c(4.5,5.1,1.6,16.1))
  plot(0,0,xlim = c(min(vect_all),max(vect_all)),ylim = yrange,
       xlab="Yield anomaly (bush/acre)",ylab="Density",type = "n",yaxt="n",cex.axis=2,cex.lab=2)
  axis(2,at=seq(0,0.1,by=0.025),cex.axis=2) #y axis labels
  points(annualproj_none,0,pch=8)
  color<-c("black",rgb(1,0,0),rgb(0,0,1),"forestgreen",rgb(1,0,1),"orange","cyan3")
  boxcolor<-c("grey",rgb(1,0,0),rgb(0,0,1,0.5),"forestgreen",rgb(1,0,1),"orange","cyan3")
  text<-c("Parameter + Structure + Climate","Parameter","Climate","Structure",
          "Parameter + Climate","Parameter + Structure","Structure + Climate")
  #plot in the same order as the 7-layer cumulative uncertainty figure result
  sorted<-sort(cumu_uncertainty)
  names<-names(sorted)
  for (i in 1:7){
    lines(density(get(paste("vect_",names[i],sep="")),na.rm=TRUE),col=color[which(cumu_uncertainty==sorted[1,i])],lwd=2)
    if (i==7) { #thicker black line for the distribution that includes 3 uncertainty sources
      lines(density(get(paste("vect_",names[i],sep="")),na.rm=TRUE),col=color[which(cumu_uncertainty==sorted[1,i])],lwd=4)
    }
    boxplot(get(paste("vect_",names[i],sep="")), horizontal = TRUE, xaxt="n",col=boxcolor[which(cumu_uncertainty==sorted[1,i])], 
            frame=F, pch=20,ylim=yrange, main="", add=TRUE,at=i*boxstep,boxwex=boxstep,cex=0.5)
    points(mean(get(paste("vect_",names[i],sep=""))),i*boxstep,pch=5,cex=1.2)
    points(density(get(paste("vect_",names[i],sep="")))$x[which(density(get(paste("vect_",names[i],sep="")))$y==max(density(get(paste("vect_",names[i],sep="")))$y))],i*boxstep,pch=19,cex=1.2)
    mtext(side=4, at=i*boxstep, line=1,text[which(cumu_uncertainty==sorted[1,i])],las=1,col=color[which(cumu_uncertainty==sorted[1,i])],cex=1.2)
  }
  
  legend(min(vect_all),yrange[2]*1.1,col=c("black","red","forestgreen","blue","orange","cyan3",rgb(1,0,1),"black"),pch=c(8,NA,NA,NA,NA,NA,NA,NA),
         lty=c(NA,1,1,1,1,1,1,1),lwd=c(NA,2,2,2,2,2,2,2),legend=c("Point estimate","Parameter","Structure","Climate","Parameter + structure ","Structure + climate","Parameter + climate","All"),bty="n",cex=2)
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






