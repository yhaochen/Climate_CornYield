#This file reads the yield projections in different model structures/climate projections and makes plots



#Read each climate projection
rm(list = ls())
graphics.off()
library(sp)
library(maps)
library(maptools)
library(usmap)
library(ncdf4)
library(precintcon)
library(BMS)
library(housingData)
library(binaryLogic)
source("latlong2county.R")
source("GDDEDD.R")
source("plots.R")
#Determine which county each grid is in s
file<-nc_open("/gpfs/group/kzk10/default/public/METDATA/raw/tmmx_1980.nc")
grid_lon<-as.vector(file$dim$lon$vals)
grid_lat<-as.vector(file$dim$lat$vals)
dim_lon<-length(grid_lon)
dim_lat<-length(grid_lat)
grid<-data.frame(x = rep(grid_lon,each=dim_lat), y = rep(grid_lat,dim_lon))
countyname<-latlong2county(grid)
countyname<-matrix(countyname,nrow=dim_lat,ncol=dim_lon)

#which counties' (grids') data are needed?
countys <- map('county', fill=TRUE, col="transparent", plot=FALSE)
IDs <- sapply(strsplit(countys$names, ":"), function(x) x[1])
countys_sp <- map2SpatialPolygons(countys, IDs=IDs,proj4string=CRS("+proj=longlat +datum=WGS84"))
countyNames <- sapply(countys_sp@polygons, function(x) x@ID)
countytoget<-countyNames[c(1:67,83:157,288:290,359:517,562:854,960:1143,1160:1183,1198:1564,1742:1762,1796:1957,2011:2098,2212:2278,2284:2329,2396:2490,2788:2887,2927:3053)]#PA,NY,NJ,MD,DE,DC,NC,VA,SC,WV,OH,MI,GA,KY,IN,IL,AL,TN,WI,MS,MN,MO,LA,AR,IA
countynum<-length(countytoget)
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

ANSI<-load("ANSI")
CountyANSI<-load("CountyANSI")
StateANSI<-load("StateANSI")
ma <- function(arr, n){  #moving avg function
  res = rep(NA,length(arr)-n+1)
  for(i in n:length(arr)){
    res[i] = mean(arr[(i-n+1):i])
  }
  res
}

su <- function(arr, n){  #cumulative sum
  res = rep(NA,length(arr)-n+1)
  for(i in n:length(arr)){
    res[i] = sum(arr[n:i])
  }
  res
}

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

#18 models
modelnames<-c("MIROC5","MRI-CGCM3","IPSL-CM5B-LR","IPSL-CM5A-LR", 
              "HadGEM2-ES365","GFDL-ESM2M","GFDL-ESM2G","CSIRO-Mk3-6-0","bcc-csm1-1",
              "MIROC-ESM", "IPSL-CM5A-MR", "CNRM-CM5","BNU-ESM",
              "MIROC-ESM-CHEM", "inmcm4", "HadGEM2-CC365", "CanESM2", "bcc-csm1-1-m")
#Metdata observation
load("Metdata_temp/Metdataframe/Data_Metobs")
Data$StateANSI<-factor(Data$StateANSI)
Data$year=Data$year+1978
Data$year<-factor(Data$year)

meanyield_anomaly<-rep(NA,32)
for (i in 1:32){
  indx<-which(Data$year==levels(Data$year)[i])
  meanyield_anomaly[i]<-mean(Data$yield_anomaly[indx],na.rm=T)
}

load("Metdata_temp/hind_bestfit")
load("Metdata_temp/proj_bestfit")
load("Metdata_temp/hind_parasample")
load("Metdata_temp/proj_parasample")
#load("Metdata_temp/proj_linearshifted_bestfit_2020_2049")
#load("Metdata_temp/proj_linearshifted_parasample_2020_2049")
load("Metdata_temp/proj_linearshifted_bestfit_2070_2099")
load("Metdata_temp/proj_linearshifted_parasample_2070_2099")
parasamplenum<-dim(proj_parasample)[3]

#Find the variance of yield distribution in these years:
selectedyears<-c(90:94) 
#Also find the upper and lower bound in each year under different choices
#all uncertainties considered
annualhistmin<-rep(NA,32)
annualhistmax<-rep(NA,32)
for (i in 1:32){
  annualhistmin[i]<-min(hind_parasample[ ,i, ],na.rm=TRUE)
  annualhistmax[i]<-max(hind_parasample[ ,i, ],na.rm=TRUE)
}
annualprojmin<-rep(NA,94)
annualprojmax<-rep(NA,94)
for (i in 1:94){
  colindx<-seq(from = i, by=94, length.out = 18)
  annualprojmin[i]<-min(proj_parasample[ ,colindx, ],na.rm=TRUE)
  annualprojmax[i]<-max(proj_parasample[ ,colindx, ],na.rm=TRUE)
}
selectedcols<-ser(selectedyears,94,18)
vect_all<-as.vector(proj_parasample[ ,selectedcols, ])
var(vect_all,na.rm = TRUE)

#projection: structure + climate 
annualprojmin_stru_clim<-rep(NA,94)
annualprojmax_stru_clim<-rep(NA,94)
for (i in 1:94){
  colindx<-seq(from = i, by=94, length.out = 18)
  annualprojmin_stru_clim[i]<-min(proj_fit[ ,colindx],na.rm=TRUE)
  annualprojmax_stru_clim[i]<-max(proj_fit[ ,colindx],na.rm=TRUE)
}
vect_stru_clim<-as.vector(proj_fit[ ,selectedcols])
var(vect_stru_clim,na.rm = TRUE)
#projection: parameter + climate
annualprojmin_para_clim<-rep(NA,94)
annualprojmax_para_clim<-rep(NA,94)
for (i in 1:94){
  colindx<-seq(from = i, by=94, length.out = 18)
  annualprojmin_para_clim[i]<-min(proj_parasample[63,colindx, ],na.rm=TRUE)
  annualprojmax_para_clim[i]<-max(proj_parasample[63,colindx, ],na.rm=TRUE)
}
vect_para_clim<-as.vector(proj_parasample[63,selectedcols, ])
var(vect_para_clim,na.rm=TRUE)
#projection: structure + parameter
annualprojmin_para_stru<-rep(NA,30)
annualprojmax_para_stru<-rep(NA,30)
for (i in 1:30){
  annualprojmin_para_stru[i]<-min(proj_linearshifted_parasample[ ,i, ],na.rm=TRUE)
  annualprojmax_para_stru[i]<-max(proj_linearshifted_parasample[ ,i, ],na.rm=TRUE)
}
vect_para_stru<- as.vector(proj_linearshifted_parasample[  ,30, ])
var(vect_para_stru,na.rm=TRUE)

#projection: only climate (best structure performance: 64th)
annualprojmax_clim<-rep(NA,94)
annualprojmin_clim<-rep(NA,94)
for (i in 1:94){
  colindx<-seq(from = i, by=94, length.out = 18)
  annualprojmin_clim[i]<-min(proj_fit[64,colindx],na.rm=TRUE)
  annualprojmax_clim[i]<-max(proj_fit[64,colindx],na.rm=TRUE)
}
vect_clim<-proj_fit[64,selectedcols]
var(vect_clim,na.rm=TRUE)
#projection: only structure (best estimate of each structure under linear shifted climate)
annualprojmin_stru<-rowMin(proj_linearshifted_fit)
annualprojmax_stru<-rowMax(proj_linearshifted_fit)
vect_stru<-as.vector(proj_linearshifted_fit[ ,30])
var(vect_stru,na.rm=TRUE)
#projection: only parameter (best structure under linear shifted climate)
annualprojmin_para<-rep(NA,30)
annualprojmax_para<-rep(NA,30)
for (i in 1:30){
  annualprojmin_para[i]<-min(proj_linearshifted_parasample[63,i, ],na.rm=TRUE)
  annualprojmax_para[i]<-max(proj_linearshifted_parasample[63,i, ],na.rm=TRUE)
}
vect_para<-as.vector(proj_linearshifted_parasample[63,30, ])
var(vect_para,na.rm=TRUE)
#projection: none (best structure's point estimates in each year under linear shifted climate)
annualproj_none<-proj_linearshifted_fit[64, ]
var(proj_linearshifted_fit[64,30])


#All
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2099),ylim = c(-100,30),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(windowyield_up,rev(windowyield_low)),col="darkseagreen1",border=NA) #window
lines(c(1981:2012),hind_fit[64, ],col="forestgreen",lwd=2.5) #none (64th model best estimate)
polygon(c(2012:2099,2099:2012),c(annualprojmax[7:94],rev(annualprojmin[7:94])),col="lightblue1",border=NA)#parametric + structure + climate
for (i in 1:13){
  indx<-c(7:94)+(i-1)*94
  lines(c(2012:2099),proj_fit[64,indx],col="blue")
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=2)
legend(1980,30,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield observation","Best model's hindcasts","Best model's projections (climate)"),bty="n",cex=1.5)
legend(2030,18, fill=c("darkseagreen1","lightblue1"),
       legend = c("Acceptable hindcast window","Yield projections (structural + parametric + climate)"),bty="n",cex=1.5)

# parameter + climate
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2099),ylim = c(-100,30),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(windowyield_up,rev(windowyield_low)),col="darkseagreen1",border=NA) #window
lines(c(1981:2012),hind_fit[64, ],col="forestgreen",lwd=2.5) #none (64th model best estimate)
polygon(c(2012:2099,2099:2012),c(annualprojmax_para_clim[7:94],rev(annualprojmin_para_clim[7:94])),col="lightblue1",border=NA)#parametric + climate
#polygon(c(2012:2099,2099:2012),c(annualprojmax_stru_clim[7:94],rev(annualprojmin_stru_clim[7:94])),col="skyblue",border=NA)#structure + climate
for (i in 1:13){
  indx<-c(7:94)+(i-1)*94
  lines(c(2012:2099),proj_fit[64,indx],col="blue")
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=2)
legend(1980,30,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield observation","Best model's hindcasts","Best model's projections (climate)"),bty="n",cex=1.5)
legend(2030,18, fill=c("darkseagreen1","lightblue1"),
       legend = c("Acceptable hindcast window","Best model's projections (parametric + climate)"),bty="n",cex=1.5)

#structure + climate
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2099),ylim = c(-100,30),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(windowyield_up,rev(windowyield_low)),col="darkseagreen1",border=NA) #window
lines(c(1981:2012),hind_fit[64, ],col="forestgreen",lwd=2.5) #none (64th model best estimate)
#polygon(c(2012:2099,2099:2012),c(annualprojmax_para_clim[7:94],rev(annualprojmin_para_clim[7:94])),col="lightblue1",border=NA)#parametric + climate
polygon(c(2012:2099,2099:2012),c(annualprojmax_stru_clim[7:94],rev(annualprojmin_stru_clim[7:94])),col="lightblue1",border=NA)#structure + climate
for (i in 1:13){
  indx<-c(7:94)+(i-1)*94
  lines(c(2012:2099),proj_fit[64,indx],col="blue")
}
points(c(1981:2012),meanyield_anomaly,col="black",pch=20,cex=2)
legend(1980,30,pch = c(20,NA,NA),lwd=c(NA,2,2),lty=c(NA,1,1),
       col=c("black","forestgreen","blue"),
       legend = c("Yield observation","Best model's hindcasts","Best model's projections (climate)"),bty="n",cex=1.5)
legend(2030,18, fill=c("darkseagreen1","lightblue1"),
       legend = c("Acceptable hindcast window","Different model's projections (structural + climate)"),bty="n",cex=1.5)



#boxplot
par(mar=c(4,5.1,1.6,2.1))
plot(0,0,xlim = c(1981,2099),ylim = c(-90,80),xlab="Year",ylab="Yield anomaly (bush/acre)",type = "n",cex.axis=2,cex.lab=2)
polygon(c(1981:2012,2012:1981),c(annualhistmax,rev(annualhistmin)),col="darkseagreen1",border=NA) # parametirc + structure
polygon(c(1981:2012,2012:1981),c(colMax(hind_fit[c(2:64), ]),rev(colMin(hind_fit[c(2:64), ]))),col="green2",border=NA)# structure
lines(c(1981:2012),SRhind,col="forestgreen",lwd=2.5) #none (S&R best estimate)
boxplot(Yield_proj~Year,data=yieldprojdata,add=TRUE,at=seq(2015,2095,by=10),axes=FALSE,boxwex=5,outpch=20)
points(c(1981:2012),meanyield,col="black",pch=20,cex=2)
legend("topright",pch = c(20,NA,NA,NA),lwd=c(NA,2,NA,NA),lty=c(NA,1,NA,NA),
       col=c("black","forestgreen",NA,NA),fill=c("white","white","darkseagreen1","green2"),
       border=c(NA,NA,"white","white"),x.intersp = c(2,2,1.3,1.3),
       legend = c("Yield observation","Best estimate of S&R model yield hindcasts","Range of model structural + parametric uncertainty",
                  "Range of only model parametric uncertainty"),bty="n",cex=1)


Years_proj<-proj_parasample
for (i in 7:94){
  colindx<-seq(from = i, by=94, length.out = 13)
  Years_proj[ ,colindx, ]<-i
}
yieldprojdata<-data.frame(Yield_proj=as.vector(proj_parasample),Year=floor((as.vector(Years_proj)+2005)/10)*10+5)
yieldprojdata<-yieldprojdata[complete.cases(yieldprojdata), ]
yieldprojdata<-yieldprojdata[which(yieldprojdata$Year>=2012), ]

#marginal yield distribution figure

i<-c(90:94)
vect1<-rep(NA,13*parasamplenum)
for (j in 1:13){
  colindx<-i+94*(j-1)
  for (k in 1: parasamplenum){
  vect1[(j-1)*parasamplenum+k]<-mean(proj_parasample[63,colindx,k],na.rm=TRUE)
  }
}

vect2<-rep(NA,13*64)
for (j in 1:13){
  colindx<-i+94*(j-1)
  for (k in 1:64){
    vect2[(j-1)*64+k]<-mean(proj_fit[k,colindx])
  }
}

vect3<-matrix(NA,nrow=13*63,ncol = parasamplenum)
for (j in 1:13){
  colindx<-i+94*(j-1)
  for (k in 1:63){
    vect3[(j-1)*63+k, ]<-colMeans(proj_parasample[k,colindx, ],na.rm=TRUE)
  }
}
vect3<-as.vector(vect3)
vect3<-vect3[!is.na(vect3)]
boxdens3(vect1,vect2,vect3,"Yield anomaly (bush/acre)", "Parametric + climate uncertainty",
        "Structural + climate uncertainty", "Parametric + structural + climate uncertainty")
points(vect1,rep(0,13),pch=20,col="red")
legend(-100,0.045,legend = "Best estimates from 13 climate projections",pch=20,col="red",bty="n",cex=1.4)


#all distributions in 2099:
par(mar=c(4.5,5.1,1.6,2.1))
plot(0,0,xlim = c(-100,20),ylim = c(-0.035,0.07),xlab="Yield anomaly (bush/acre)",ylab="Density",type = "n",cex.axis=2,cex.lab=2)
points(proj_linearshifted_fit[64,38],0,pch=8)
lines(density(vect_para,na.rm=TRUE),col=rgb(1,0,0),lwd=2)
lines(density(vect_stru,na.rm=TRUE),col=rgb(0,1,0),lwd=2)
lines(density(vect_clim,na.rm=TRUE),col=rgb(0,0,1),lwd=2)
lines(density(vect_para_stru,na.rm=TRUE),col=rgb(1,0.75,0),lwd=2,lty=2)
lines(density(vect_stru_clim,na.rm=TRUE),col=rgb(0,1,1),lwd=2,lty=2)
lines(density(vect_para_clim,na.rm=TRUE),col=rgb(1,0,1),lwd=2,lty=2)
lines(density(vect_all,na.rm=TRUE),col="black",lwd=3)
boxplot(vect_para, horizontal = TRUE, xaxt="n",col=rgb(1,0,0), 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.005,boxwex=0.003) 
boxplot(vect_stru, horizontal = TRUE, xaxt="n",col=rgb(0,1,0), 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.01,boxwex=0.003) 
boxplot(vect_clim, horizontal = TRUE, xaxt="n",col=rgb(0,0,1), 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.015,boxwex=0.003) 
boxplot(vect_para_stru, horizontal = TRUE, xaxt="n",col=rgb(1,0.75,0), 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.02,boxwex=0.003) 
boxplot(vect_stru_clim, horizontal = TRUE, xaxt="n",col=rgb(0,1,1), 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.025,boxwex=0.003) 
boxplot(vect_para_clim, horizontal = TRUE, xaxt="n",col=rgb(1,0,1), 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.03,boxwex=0.003) 
boxplot(vect_all, horizontal = TRUE, xaxt="n",col="grey", 
        frame=F, pch=20,ylim=c(-100,20), main="", add=TRUE,at=-0.035,boxwex=0.003) 

legend(-103,0.08,col=c("black","red","blue","green"),pch=c(8,NA,NA,NA),lty=c(NA,1,1,1),legend=c("Point estimate",
        "Parametric","Climatic","Structural"),bty="n",cex=1.5)
legend(-70,0.08,col=c(rgb(1,0.78,0),rgb(0,1,1),rgb(1,0,1),"black"),lty=c(2,2,2,1),legend=c("Parametric + structural ",
        "Structural + climatic","Parametric + climatic","All"),bty="n",cex=1.5)
# #2D paremeter distributions
parasamplenum<-10000
i<-16
j<-4
model<-lm(as.formula(paste("yield_anomaly ~ ",paste(Tnames[i], Pnames[j],sep="+"), sep="") ),data=Data, weights=Data$area)
MCpara<-parasample(model,parasamplenum,2)
MCpara_keep<-MCpara
hind_parasample<-matrix(NA,nrow=32,ncol=10000)
variablenames<-variable.names(model)
variablenum<-length(variablenames)

#for hindcast of each parametric sample
col_data_hind<-rep(NA,variablenum-1) #first variable is intercept
for (m in 1:(variablenum-1)){
  col_data_hind[m]<-which(colnames(Data)==variablenames[m+1])
}

hind<-predict(model,Data)
for (k in 1:32){
  indx<-which(Data$year==levels(Data$year)[k])
  for (n in 1:parasamplenum){
    hind_parasample[k,n]<-MCpara[n,1]+MCpara[n,2:variablenum]%*%colMeans(as.matrix(Data[indx,col_data_hind]))
  }
}
#keep the models if hindcast pass the window
for (n in 1:parasamplenum){
  for (k in 1:32){
    if ((hind_parasample[k,n] > windowyield_up[k]) || (hind_parasample[k,n] < windowyield_low[k])){
      MCpara_keep[n, ] <- NA
    }
  }
}

par(mfrow=c(7,7),mar=c(1,1,1,1))
for (i in 1:7){
  for (j in 1:7){
    if (j<i){
      d<-data.frame(x=MCpara[ ,j],y=MCpara[ ,i])
      kd <- ks::kde(d, compute.cont=TRUE)
      contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                          z=estimate, levels=cont["5%"])[[1]])
      contour_95 <- data.frame(contour_95)
      
      plot(d,pch=20,cex=0.5,col="red")
      points(MCpara_keep[ ,j],MCpara_keep[ ,i],pch=20,cex=0.5,col="grey")
      lines(contour_95$x,contour_95$y,lwd=2,col="blue")
    }else {
      plot(0,type='n',axes=FALSE,ann=FALSE)
    }
  }
}




