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
#3-layer std
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
  if (k==1){
    text(0,3.3,"(a).",cex=2)
  }
  if (k==2){
    text(0,3.3,"(b).",cex=2)
  }
  text(right,3.3,"Cumulative uncertainty",cex=2)
  if (k==1){
    text(max(sorted)/2,3.3,"2020 - 2049",cex=2)
  }
  if (k==2){
    text(max(sorted)/2,3.3,"2070 - 2099",cex=2)
  }
  for (i in 1:3){
    text(right,i-0.5,paste(sprintf("%.1f",sorted[1,i])," (",sprintf("%.0f",sorted[1,i]/max(sorted)*100),"%)",sep=""),
         cex=2,col=textcolor[which(cumu_uncertainty==sorted[1,i])])
    text(max(sorted)/2,i-0.5,text[which(cumu_uncertainty==sorted[1,i])],cex=2,col=textcolor[which(cumu_uncertainty==sorted[1,i])])
  }
}
dev.off()




