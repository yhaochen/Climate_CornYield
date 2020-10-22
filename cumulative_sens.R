#plot of cumulative sensitivity

#cumulative uncertainty (based on Jim et al 2019)
rm(list = ls())
graphics.off()

#read 30-year mean yield data and calculate cumulative uncertainty at each stage
for (k in 1:2){
  if (k==1){
    load("Metdata/proj_30yravg_2020_2049")
  }
  if (k==2){
    load("Metdata/proj_30yravg_2070_2099")
  }
  climnum<-dim(proj_30yravg)[1]
  parasamplenum<-dim(proj_30yravg)[2]
  #Level one: 
  std_para<-rep(NA,climnum)
  for (i in 1:climnum){
    std_para[i]<-var(proj_30yravg[i, ])*(parasamplenum-1)/parasamplenum #population var instead of sample var
  }
  std_marg_para<-sqrt(mean(std_para))

  std_clim<-rep(NA,parasamplenum)
  for (i in 1:parasamplenum){
    std_clim[i]<-var(proj_30yravg[ ,i])*(climnum-1)/climnum
  }
  std_marg_clim<-sqrt(mean(std_clim))

  #Level two:
  std_para_clim<-sqrt(var(as.vector(proj_30yravg)))*(climnum*parasamplenum-1)/(climnum*parasamplenum)
  cumu_uncertainty<-data.frame(para=std_marg_para,clim=std_marg_clim, para_clim=std_para_clim)
  
  if (k==1){
    save(cumu_uncertainty,file="Metdata/cumu_uncertainty_2020_2049")
  }

  if (k==2){
    save(cumu_uncertainty,file="Metdata/cumu_uncertainty_2070_2099")
  }
}


png("Plots/3-layer.png", width = 1000, height = 720)
#7-layer std
par(mfrow=c(2,1))
par(mar=c(4.5,5.1,1.6,2.1))
for (k in 1:2){
  if (k==1){
    load("Metdata/cumu_uncertainty_2020_2049")
  }
  if (k==2){
    load("Metdata/cumu_uncertainty_2070_2099")
  }
  color<-c(rgb(1,0,0,0.3),rgb(0,1,0,0.3),rgb(0,0,1,0.3))
  text<-c("Parameter","Climate","Parameter + Climate")
  textcolor<-c(rgb(1,0,0),"forestgreen",rgb(0,0,1))
  #sort cumulative uncertainties
  sorted<-sort(cumu_uncertainty)

  plot(0,0,xlim = c(0,max(sorted)*1.4),ylim = c(0,3.5),xlab="",
       ylab=" ",type = "n",axes="FALSE",cex.lab=2)
  axis(1,at=c(0,max(sorted)),labels=c("0.0",sprintf("%.1f",max(sorted))),cex.axis=2)
  for (i in 1:3){
    polygon(c((max(sorted)-sorted[1,i])/2, (max(sorted)+sorted[1,i])/2,
              (max(sorted)+sorted[1,i])/2, (max(sorted)-sorted[1,i])/2),
            c(i,i,i-1,i-1),col=color[which(cumu_uncertainty==sorted[1,i])],border=NA)
  }
  mtext(side=1, at=max(sorted)/2,line=1.5,"Yield anomaly (standard deviation)",cex=2) #serve as x axis label
  right=max(sorted)*1.2
  text(right,7.5,"Cumulative uncertainty",cex=2)
  
  for (i in 1:3){
    text(right,i-0.5,paste(sprintf("%.1f",sorted[1,i])," (",sprintf("%.1f",sorted[1,i]/max(sorted)*100),"%)",sep=""),
         cex=2,col=textcolor[which(cumu_uncertainty==sorted[1,i])])
    text(max(sorted)/2,i-0.5,text[which(cumu_uncertainty==sorted[1,i])],cex=2,col=textcolor[which(cumu_uncertainty==sorted[1,i])])
  }
}
dev.off()





# 
# #3-layer both
# par(mar=c(4.5,5.1,1.6,2.1))
# plot(0,0,xlim = c(min(proj_30yravg_2070_2099)-30,max(proj_30yravg_2070_2099)+30),ylim = c(0,4),xlab="Yield anomaly (range: bush/acre)",
#      ylab=" ",type = "n",axes="FALSE",cex.lab=2)
# axis(1,c(-79.4,27.4),cex.axis=2)
# #axis(1,c(-40.9,22.7),cex.axis=2)
# range_2070_2099<-c(min(proj_30yravg_2070_2099),max(proj_30yravg_2070_2099))
# polygon(c(range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para[1])/2,
#           range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para[1])/2),
#         c(1,1,0,0),col=rgb(1,0,0,0.3),border=NA)
# polygon(c(range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para_stru[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para_stru[1])/2,
#           range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para_stru[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para_stru[1])/2),
#         c(2,2,1,1),col=rgb(0,1,0,0.3),border=NA)
# polygon(c(range_2070_2099[1],range_2070_2099[2],range_2070_2099[2],range_2070_2099[1]), c(3,3,2,2),col=rgb(0,0,1,0.3),border=NA)
# 
# range_2070_2099<-c(min(proj_30yravg_2070_2099),max(proj_30yravg_2070_2099))
# polygon(c(range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para[1])/2,
#           range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para[1])/2),
#         c(1,1,0,0),col=rgb(1,0,0,0.3),border=NA)
# polygon(c(range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para_stru[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para_stru[1])/2,
#           range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]+cumu_uncertainty_2070_2099$para_stru[1])/2,range_2070_2099[1]+(cumu_uncertainty_2070_2099$all[1]-cumu_uncertainty_2070_2099$para_stru[1])/2),
#         c(2,2,1,1),col=rgb(0,1,0,0.3),border=NA)
# polygon(c(range_2070_2099[1],range_2070_2099[2],range_2070_2099[2],range_2070_2099[1]), c(3,3,2,2),col=rgb(0,0,1,0.3),border=NA)
# 
# 
# 
# text(-90,3.5,"Stage uncertainty",cex=2)
# text(34,3.5,"Cumulative uncertainty",cex=2)
# text(-99,0.5,"34.5 (33.6%)",cex=2,col=rgb(1,0,0))
# text(-99,1.5,"42.5 (41.3%)",cex=2,col="forestgreen")
# text(-99,2.5,"29.8 (25.1%)",cex=2,col=rgb(0,0,1))
# text(47,0.5,"34.5 (33.6%)",cex=2,col=rgb(1,0,0))
# text(47,1.5,"77.0 (74.9%)",cex=2,col="forestgreen")
# text(47,2.5,"106.8 (100%)",cex=2,col=rgb(0,0,1))
# text((27.4-79.4)/2,0.5,"Parameter",cex=2,col=rgb(1,0,0))
# text((27.4-79.4)/2,1.5,"Parameter + Structure",cex=2,col="forestgreen")
# text((27.4-79.4)/2,2.5,"Parameter + Structure + Climate",cex=2,col=rgb(0,0,1))
# 
# 
# 
# 
# par(mar=c(4.5,5.1,1.6,3.1))
# plot(0,0,xlim = c(10,70),ylim = c(0,3),xlab="Cumulative uncertainty (range)",
#      ylab=" ",type = "n",axes="FALSE",cex.lab=2)
# axis(1,cex.axis=2)
# segments(x0=range_all,y0=-0.2,x1=range_all,y1=3,lwd=2)
# segments(x0=range_marg_para_clim,y0=-0.2,x1=range_marg_para_clim,y1=2,col=rgb(1,0,1,0.3),lwd=2)
# segments(x0=range_marg_para_stru,y0=-0.2,x1=range_marg_para_stru,y1=2,col=rgb(1,1,0,0.3),lwd=2)
# segments(x0=range_marg_stru_clim,y0=-0.2,x1=range_marg_stru_clim,y1=2,col=rgb(0,1,1,0.3),lwd=2)
# segments(x0=range_marg_para,y0=-0.2,x1=range_marg_para,y1=1,col=rgb(1,0,0,0.3),lwd=2)
# segments(x0=range_marg_stru,y0=-0.2,x1=range_marg_stru,y1=1,col=rgb(0,1,0,0.3),lwd=2)
# segments(x0=range_marg_clim,y0=-0.2,x1=range_marg_clim,y1=1,col=rgb(0,0,1,0.3),lwd=2)
# 
# # 3-layer plot
# par(mar=c(4.5,5.1,1.6,2.1))
# plot(0,0,xlim = c(min(proj_30yravg)-30,max(proj_30yravg)+30),ylim = c(0,4),xlab="Yield anomaly (range: bush/acre)",
#      ylab=" ",type = "n",axes="FALSE",cex.lab=2)
# axis(1,c(-79.4,27.4),cex.axis=2)
# #axis(1,c(-40.9,22.7),cex.axis=2)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_para)/2,min(proj_30yravg)+(range_all+range_marg_para)/2,
#           min(proj_30yravg)+(range_all+range_marg_para)/2,min(proj_30yravg)+(range_all-range_marg_para)/2),
#         c(1,1,0,0),col=rgb(1,0,0,0.3),border=NA)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_para_stru)/2,min(proj_30yravg)+(range_all+range_marg_para_stru)/2,
#           min(proj_30yravg)+(range_all+range_marg_para_stru)/2,min(proj_30yravg)+(range_all-range_marg_para_stru)/2),
#         c(2,2,1,1),col=rgb(0,1,0,0.3),border=NA)
# polygon(c(min(proj_30yravg),max(proj_30yravg),
#           max(proj_30yravg),min(proj_30yravg)), c(3,3,2,2),col=rgb(0,0,1,0.3),border=NA)
# text(-90,3.5,"Stage uncertainty",cex=2)
# text(34,3.5,"Cumulative uncertainty",cex=2)
# text(-99,0.5,"34.5 (33.6%)",cex=2,col=rgb(1,0,0))
# text(-99,1.5,"42.5 (41.3%)",cex=2,col="forestgreen")
# text(-99,2.5,"29.8 (25.1%)",cex=2,col=rgb(0,0,1))
# text(47,0.5,"34.5 (33.6%)",cex=2,col=rgb(1,0,0))
# text(47,1.5,"77.0 (74.9%)",cex=2,col="forestgreen")
# text(47,2.5,"106.8 (100%)",cex=2,col=rgb(0,0,1))
# text((27.4-79.4)/2,0.5,"Parameter",cex=2,col=rgb(1,0,0))
# text((27.4-79.4)/2,1.5,"Parameter + Structure",cex=2,col="forestgreen")
# text((27.4-79.4)/2,2.5,"Parameter + Structure + Climate",cex=2,col=rgb(0,0,1))
# 
# 
# 
# par(mar=c(4.5,5.1,1.6,2.1))
# plot(0,0,xlim = c(-5,std_all+5),ylim = c(0,4),xlab="Yield anomaly (standard deviation)",
#      ylab=" ",type = "n",axes="FALSE",cex.lab=2)
# axis(1,labels = FALSE,at=c(0,std_all))
# polygon(c((std_all-std_marg_para)/2, (std_all+std_marg_para)/2,
#           (std_all+std_marg_para)/2, (std_all-std_marg_para)/2),
#         c(1,1,0,0),col=rgb(1,0,0,0.3),border=NA)
# polygon(c((std_all-std_marg_para_stru)/2, (std_all+std_marg_para_stru)/2,
#           (std_all+std_marg_para_stru)/2, (std_all-std_marg_para_stru)/2),
#         c(2,2,1,1),col=rgb(0,1,0,0.3),border=NA)
# polygon(c(0,std_all,std_all,0), c(3,3,2,2),col=rgb(0,0,1,0.3),border=NA)
# text(-1,3.5,"Stage uncertainty",cex=2)
# text(17,3.5,"Cumulative uncertainty",cex=2)
# text(-3,0.5,"6.1 (75.6%)",cex=2,col=rgb(1,0,0))
# text(-3,1.5,"1.4 (17.5%)",cex=2,col="forestgreen")
# text(-3,2.5,"0.6 (6.9%)",cex=2,col=rgb(0,0,1))
# text(20,0.5,"6.1 (75.6%)",cex=2,col=rgb(1,0,0))
# text(20,1.5,"7.5 (93.1)",cex=2,col="forestgreen")
# text(20,2.5,"8.1 (100%)",cex=2,col=rgb(0,0,1))
# text(std_all/2,0.5,"Parameter",cex=2,col=rgb(1,0,0))
# text(std_all/2,1.5,"Parameter + Structure",cex=2,col="forestgreen")
# text(std_all/2,2.5,"Parameter + Structure + Climate",cex=2,col=rgb(0,0,1))
# 
# 
# 
# 
# #7-layer fig range
# par(mar=c(4.5,5.1,1.6,2.1))
# plot(0,0,xlim = c(min(proj_30yravg)-30,max(proj_30yravg)+30),ylim = c(0,8),xlab="Yield anomaly (range: bush/acre)",
#      ylab=" ",type = "n",axes="FALSE",cex.lab=2)
# range=c(-79.4,27.4)
# #range=c(-40.9,22.7)
# axis(1,range,cex.axis=2)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_clim)/2,min(proj_30yravg)+(range_all+range_marg_clim)/2,
#           min(proj_30yravg)+(range_all+range_marg_clim)/2,min(proj_30yravg)+(range_all-range_marg_clim)/2),
#         c(1,1,0,0),col=rgb(0,0,1,0.3),border=NA)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_para)/2,min(proj_30yravg)+(range_all+range_marg_para)/2,
#           min(proj_30yravg)+(range_all+range_marg_para)/2,min(proj_30yravg)+(range_all-range_marg_para)/2),
#         c(2,2,1,1),col=rgb(1,0,0,0.3),border=NA)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_stru)/2,min(proj_30yravg)+(range_all+range_marg_stru)/2,
#           min(proj_30yravg)+(range_all+range_marg_stru)/2,min(proj_30yravg)+(range_all-range_marg_stru)/2),
#         c(3,3,2,2),col=rgb(0,1,0,0.3),border=NA)
# 
# 
# polygon(c(min(proj_30yravg)+(range_all-range_marg_para_clim)/2,min(proj_30yravg)+(range_all+range_marg_para_clim)/2,
#           min(proj_30yravg)+(range_all+range_marg_para_clim)/2,min(proj_30yravg)+(range_all-range_marg_para_clim)/2),
#         c(4,4,3,3),col=rgb(1,0,1,0.3),border=NA)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_stru_clim)/2,min(proj_30yravg)+(range_all+range_marg_stru_clim)/2,
#           min(proj_30yravg)+(range_all+range_marg_stru_clim)/2,min(proj_30yravg)+(range_all-range_marg_stru_clim)/2),
#         c(6,6,5,5),col=rgb(0,1,1,0.3),border=NA)
# polygon(c(min(proj_30yravg)+(range_all-range_marg_para_stru)/2,min(proj_30yravg)+(range_all+range_marg_para_stru)/2,
#           min(proj_30yravg)+(range_all+range_marg_para_stru)/2,min(proj_30yravg)+(range_all-range_marg_para_stru)/2),
#         c(5,5,4,4),col=rgb(1,1,0,0.3),border=NA)
# 
# polygon(c(min(proj_30yravg),max(proj_30yravg),
#           max(proj_30yravg),min(proj_30yravg)), c(7,7,6,6),col=rgb(0,0,0,0.3),border=NA)
# 
# text(-85,7.5,"Stage uncertainty",cex=2)
# text(32,7.5,"Cumulative uncertainty",cex=2)
# left=-97
# text(left,0.5,"30.0 (28.1%)",cex=2,col=rgb(0,0,1))
# text(left,1.5,"4.5 (4.2%)",cex=2,col=rgb(1,0,0))
# text(left,2.5,"21.9 (20.6%)",cex=2,col="forestgreen")
# text(left,3.5,"5.7 (5.3%)",cex=2,col=rgb(1,0,1))
# text(left,4.5,"14.9 (13.9%)",cex=2,col="orange")
# text(left,5.5,"8.9 (8.4%)",cex=2,col="cyan3")
# 
# text(left,6.5,"20.7 (19.5%)",cex=2,col=rgb(0,0,0))
# right=45
# text(right,0.5,"30.0 (28.1%)",cex=2,col=rgb(0,0,1))
# text(right,1.5,"34.5 (32.3%)",cex=2,col=rgb(1,0,0))
# text(right,2.5,"56.4 (52.9%)",cex=2,col="forestgreen")
# text(right,3.5,"62.1 (58.2%)",cex=2,col=rgb(1,0,1))
# text(right,4.5,"77.0 (72.1%)",cex=2,col="orange")
# text(right,5.5,"85.9 (80.5%)",cex=2,col="cyan3")
# 
# text(right,6.5,"106.6 (100%)",cex=2,col=rgb(0,0,0))
# 
# text(mean(range),0.5,"Climate",cex=2,col=rgb(0,0,1))
# text(mean(range),1.5,"Parameter",cex=2,col=rgb(1,0,0))
# text(mean(range),2.5,"Structure",cex=2,col="forestgreen")
# text(mean(range),3.5,"Parameter + Climate",cex=2,col=rgb(1,0,1))
# text(mean(range),4.5,"Parameter + Structure",cex=2,col="orange")
# text(mean(range),5.5,"Structure + Climate",cex=2,col="cyan3")
# text(mean(range),6.5,"Parameter + Structure + Climate",cex=2,col=rgb(0,0,0))
