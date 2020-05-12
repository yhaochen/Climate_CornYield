# This file uses Bayesian regression (MCMC) for PRISM data

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
load("Data_prism")
Data_prism<-Data_prism[complete.cases(Data_prism), ]
Data_prism$StateANSI<-factor(Data_prism$StateANSI)
Data_prism$year<-factor(Data_prism$year)
#test with partial data
Datatouse=Data_prism#[1:1000, ]
Yield<-Datatouse$Yield_anomaly
GDD<-Datatouse$GDD
EDD<-Datatouse$EDD
SPI<-Datatouse$SPI
year=32
#MCMC
# posterior = prior*likelihood
# input theta is the parameter space including the parameters required for likelihood function (sigma for normal distribution)




# inside the likelihood function, add: 1: temporal fixed effects, quadratic precip & yield in state
#2: allowing for spatial correlation in residuals
logp = function(theta){
  N = dim(Datatouse)[1] #data size
  # Calulate model simulations.
  ymean=theta[1]; gdd = theta[2]; edd = theta[3]; spi=theta[4]; sigma=theta[5]
  model = GDD*gdd+EDD*edd+SPI*spi+ymean # model to use
  # Estimate the residuals (i.e., the deviation of the observations from the
  # model simulation).
  resid = Yield - model
  # Get the log of the likelihood function.
  log.likelihood = -N/2*log(2*pi) - N*log(sigma) - 1/2*sum(resid^2/sigma^2)
  # Use an improper uniform "relatively uninformative" prior.
  log.prior = 0 # log(1)
  # Bayesian updating: update the probability estimates based on the observations.
  log.posterior = log.likelihood + log.prior
  # Return the unnormalized log posterior value.
  return(log.posterior)
}

#adaptive mcmc
accept.mcmc = 0.234
gamma.mcmc = 0.6
NI <- 1e5
step <- c(0.02, 0.002, 0.01, 0.1, 0.02) #GDD,EDD,SPI,yieldmean,sigma
p<-lm(Yield_anomaly~GDD+EDD+SPI,data=Data_prism)
p0 <- c(unname(p$coefficients), 16) #initial values: linear regression
# Run MCMC calibration. Note: Runs for several minutes.
mcmc.out <- MCMC(logp, NI, p0, scale = step, adapt = TRUE, acc.rate = accept.mcmc,
                 gamma = gamma.mcmc, list = TRUE, n.start = round(0.01*NI))
mcmc.chains <- mcmc.out$samples

par(mfrow = c(2,3))
plot(c(1:NI),mcmc.chains[ ,1],type="l")
plot(c(1:NI),mcmc.chains[ ,2],type="l")
plot(c(1:NI),mcmc.chains[ ,3],type="l")
plot(c(1:NI),mcmc.chains[ ,4],type="l")
plot(c(1:NI),mcmc.chains[ ,5],type="l")
par(mfrow = c(2,3))
plot(density(mcmc.chains[c(10000:NI),2]),xlab="",ylab="Density",main="GDD",cex.lab=1.3,cex.axis=1.3)
plot(density(mcmc.chains[c(10000:NI),3]),xlab="",ylab="",main="EDD",cex.lab=1.3,cex.axis=1.3)
plot(density(mcmc.chains[c(10000:NI),4]),xlab="",ylab="",main="SPI",cex.lab=1.3,cex.axis=1.3)
plot(density(mcmc.chains[c(10000:NI),1]),xlab="",ylab="",main="Yield mean",cex.lab=1.3,cex.axis=1.3)
plot(density(mcmc.chains[c(10000:NI),5]),xlab="",ylab="",main="error variance",cex.lab=1.3,cex.axis=1.3)

theta_best<-c(mean(mcmc.chains[c(10000:NI),1]),mean(mcmc.chains[c(10000:NI),2]),mean(mcmc.chains[c(10000:NI),3]),mean(mcmc.chains[c(10000:NI),4]),mean(mcmc.chains[c(10000:NI),5]))

#hindcast and res
hind<-GDD*theta_best[2]+EDD*theta_best[3]+SPI*theta_best[4]+theta_best[1]
res<-Yield-hind
Data_prism$hindcast<-hind
Data_prism$residual<-res
plot(density(res))
save(Data_prism,file="Data_prism_complete")
indx<-which((Data_prism$StateANSI==1)&(Data_prism$countyANSI==5))
acf(res[indx],main="")

#task1: plot temporal residuals, spatial residuals
startingindex<-which(Data_prism$year==1) # want to find the counties with complete record (32 years), first find all years with "1"
loc<-which(diff(startingindex)==32)      # if two "1"s are 32 years away, it means a complete record
loc<-c(loc,length(startingindex))        # also test the last "1" mannually, an it is a complete record so add it to the end of the location vector
startingindex<-startingindex[loc]        # now this is the starting index of the counties we want
complete_countynum<-length(startingindex)
ACF<-rep(NA,complete_countynum)
PACF<-rep(NA,complete_countynum)
for (i in 1:complete_countynum){
  restoget<-Data_prism$resfix[startingindex[i]:(startingindex[i]+31)]
  a<-acf(restoget)
  ACF[i]<-a$acf[2]
  b<-pacf(restoget)
  PACF[i]<-b$acf[1]
}
par(mar=c(5.1,5.1,3.1,2.1))
plot(density(ACF),xlab="AR(1) parameter",ylab="Density",main="",cex.lab=1.5,cex.axis=1.5)  #ACF distribution for all counties (with 32 years record)
abline(v=mean(ACF),col="blue",lty=2)

#task2: heat map showing spatial autocorrelation: plot residuals at each county with a color
yeartoget<-1
index<-which(Data_prism$year==yeartoget)
newdata$resfix[which(newdata$resfix< -40)]=-40
newdata$resfix[which(newdata$resfix>40)]=40

plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"), #for res
           data=newdata[index ,c(15,12)], values = "resfix") + labs(title = paste("Yield residual of each county in ",yeartoget+1980,sep="")) + 
  scale_fill_gradient2(low = "blue",mid="white", high ="red",midpoint = 0, limits=c(-40,40),name="bush/acre")+ 
  theme(plot.title = element_text(size=14))
plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),  # for Pr
           data=newdata[index ,c(7,12)], values = "Pr") + labs(title = paste("Precipitation of each county in ",yeartoget+1980,sep="")) + 
  scale_fill_gradient2(low = "blue",mid="white", high ="red",midpoint = 100, limits=c(0,200),name="mm")+ 
  theme(plot.title = element_text(size=14))
plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),  # for GDD
           data=newdata[index ,c(3,12)], values = "GDD") + labs(title = paste("GDD of each county in ",yeartoget+1980,sep="")) + 
  scale_fill_gradient2(low = "blue",mid="white", high ="red",midpoint = 1700, limits=c(700,2700),name="degrees Celcius")+ 
  theme(plot.title = element_text(size=14))
plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),  # for EDD
           data=newdata[index ,c(5,12)], values = "EDD") + labs(title = paste("EDD of each county in ",yeartoget+1980,sep="")) + 
  scale_fill_gradient2(low = "blue",mid="white", high ="red",midpoint = 150, limits=c(0,300),name="degrees Celcius")+ 
  theme(plot.title = element_text(size=14))
#variance of residuals
resvar<-rep(NA,length(levels(newdata$fips)))
resmean<-rep(NA,length(levels(newdata$fips)))
for (i in 1:length(levels(newdata$fips))) {
  index<-which(newdata$fips==levels(newdata$fips)[i])
  resvar[i]<-var(newdata$resfix[index])
  resmean[i]<-mean(newdata$resfix[index],na.rm=TRUE)
}
resframe<-data.frame(fips=levels(newdata$fips),var=resvar,mean=resmean,std=sqrt(resvar))
plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
           data=resframe, values = "std") + labs(title = "Corn yield residuals standard error") + 
  scale_fill_gradient2(low = "blue",mid="white", high ="red", limits=c(0,20),name="Residual standard error (bu/acre)") +theme(plot.title = element_text(size=24))


#task3: variogram test exponential vs iid
geoCounty
hind<-newdata$hindfix
Data_prism$residual<-newdata$resfix
Data_prism$lat<-rep(NA,length(hind))
Data_prism$lon<-rep(NA,length(hind))
levelind<-as.integer(Data_prism$fips)
for (i in 1:length(hind)){
  c<-levels(Data_prism$fips)[levelind[i]] #which level (county)
  index<-which(geoCounty==c)
  Data_prism$lat[i]<-geoCounty[index,5]
  Data_prism$lon[i]<-geoCounty[index,4]
}
index<-which(Data_prism$year==yeartoget)
geoobj<-matrix(NA,nrow=length(index),ncol=3)
geoobj[ ,1]<-Data_prism$residual[index]
geoobj[ ,2]<-Data_prism$lon[index]
geoobj[ ,3]<-Data_prism$lat[index]
#PM <- as.geodata(geoobj, coords.col = 2:3, data.col = 1)  
#v1=variog(PM,option = "smooth")

PM<-data.frame(resmean=Data_prism$residual[index],lon=111*(geoobj[ ,2]+82)*cos(geoobj[ ,3]),lat=111*(geoobj[ ,3]-37))
coordinates(PM)<-~lon+lat
V<-variogram(resmean~1,locations=PM$lon+PM$lat,data=PM,width=10,cutoff=600,alpha=c(0,45,90,135))
Vfit<-vgm(psill=165, model="Sph", nugget=350, range=100)
plot(V,model=Vfit,pch=20,col="black",xlab="Distance (km)",ylab="Semivariance",cex.lab=1.4,cex.axis=1.4)


par(mar=c(5.1,5.1,3.1,2.1))
plot(v1$u,v1$v,type="l",xlab="Distance (deg)",ylab="semivariance",xlim=c(0,23),cex.lab=1.6,cex.axis=1.6)
abline(h=var(Data_prism$residual[index]),lty=2)
legend("topleft",col=c("black","black","red"),lty = c(1,2,1),lwd=c(1.9,1,1),legend=c("Semivariance of resifuals","Theoretical semivariance for iid residuals","Fitted logarithm function"),bty="n",cex=1.4)


#Bayes factor
#sample points from variogram result, smooth variogram gives 760000 points, sample every 100
dist<-v1$u[seq(1, length(v1$u), 100)]
semiva<-v1$v[seq(1,length(v1$v), 100)]
variogframe<-data.frame(x=dist,y=semiva)
exponential.model <- lm(semiva ~ log(dist))
newdata=dist
a=predict(exponential.model,dist=newdata)
lines(dist,a,type="l",col="red")
res_iidcase<-semiva-mean(semiva)
res_logcase<-as.vector(semiva-a)


logp_iid = function(theta){
  N = length(semiva) #data size
  # Calulate model simulations.
  ymean=theta[1]; sigma=theta[2]
  model = ymean # model to use
  # Estimate the residuals (i.e., the deviation of the observations from the
  # model simulation).
  resid = semiva-model
  # Get the log of the likelihood function.
  log.likelihood = -N/2*log(2*pi) - N*log(sigma) - 1/2*sum(resid^2/sigma^2)
  # Use an improper uniform "relatively uninformative" prior.
  log.prior = 0 # log(1)
  # Bayesian updating: update the probability estimates based on the observations.
  log.posterior = log.likelihood + log.prior
  # Return the unnormalized log posterior value.
  return(log.posterior)
}
logp_log = function(theta){
  N = length(semiva) #data size
  # Calulate model simulations.
  k=theta[1];b=theta[2];sigma=theta[3]
  model = k*log(dist)+b# model to use
  # Estimate the residuals (i.e., the deviation of the observations from the
  # model simulation).
  resid = semiva-model
  # Get the log of the likelihood function.
  log.likelihood = -N/2*log(2*pi) - N*log(sigma) - 1/2*sum(resid^2/sigma^2)
  # Use an improper uniform "relatively uninformative" prior.
  log.prior = 0 # log(1)
  # Bayesian updating: update the probability estimates based on the observations.
  log.posterior = log.likelihood + log.prior
  # Return the unnormalized log posterior value.
  return(log.posterior)
}
accept.mcmc = 0.234
gamma.mcmc = 0.6
NI <- 20000
step <- c(1,0.6) #mean, sigma
p0 <- c(mean(semiva),sqrt(var(semiva))) #initial values: linear regression
# Run MCMC calibration. Note: Runs for several minutes.
mcmc.out <- MCMC(logp_iid, NI, p0, scale = step, adapt = TRUE, acc.rate = accept.mcmc,
                 gamma = gamma.mcmc, list = TRUE, n.start = round(0.01*NI))
mcmc.chains <- mcmc.out$samples

accept.mcmc = 0.234
gamma.mcmc = 0.6
NI <- 20000
step <- c(0.1,0.1,1) #k,b,sigma
p0 <- c(52,224,91) #initial values: linear regression
# Run MCMC calibration. Note: Runs for several minutes.
mcmc.out <- MCMC(logp_log, NI, p0, scale = step, adapt = TRUE, acc.rate = accept.mcmc,
                 gamma = gamma.mcmc, list = TRUE, n.start = round(0.01*NI))
mcmc.chains <- mcmc.out$samples

#plot ar1 parameters for each county
index<-which(Data_prism$fips=="42027")
resseries<-Data_prism$resfix[index]
arima(resseries, order=c(1,0,0))
par(mar=c(5.1,5.1,3.1,2.1))
plot(c(1981:2012),Data_prism$Yield[index],xlab="Year",ylab="Yield (bu/acre)",cex.lab=1.5,cex.axis=1.5,pch=20)
lines(c(1981:2012),Data_prism$hindfix[index],col="blue")
legend("topleft",pch = c(20,NA),lty = c(NA, 1),col=c("black","blue"),legend=c("Yield observation in Center county, PA","Yield hindcast"),bty="n",cex=1.2)
par(mar=c(5.1,5.1,3.1,2.1))
plot(c(1981:2012),Data_prism$resfix[index],xlab="Year",ylab="Yield residuals (Bu/acre)",cex.lab=1.5,cex.axis=1.5,pch=20)
plot(density(Data_prism$resfix[index]),xlab="Year",ylab="Density",main="",cex.axis=1.4,cex.lab=1.5)
acf(Data_prism$resfix[index],main="")
pacf(Data_prism$resfix[index],xlab="Lag (year)",main="")
qqnorm(Data_prism$resfix[index],pch=20,cex.axis=1.4,cex.lab=1.4,main="Normal Q-Q plot of yield residuals in Center county, PA",cex.main=1.4)
qqline(Data_prism$resfix[index],col="blue",lwd=2)