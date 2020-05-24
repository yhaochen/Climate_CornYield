rm(list = ls())
graphics.off()
source("plots.R")
source("GDDEDD.R")
#maca-met projections 2066-2099 34 tot
load("Data/MACAv2-METDATA/GDD_all")
load("Data/MACAv2-METDATA/EDD_all")
load("Data/MACAv2-METDATA/Tmax_all")
load("Data/MACAv2-METDATA/Tmin_all")
load("Data/MACAv2-METDATA/Pr_all")
#maca-met historical 1970-2004 35 tot
load("Data/MACAv2-METDATA/GDD_historical_all")
load("Data/MACAv2-METDATA/EDD_historical_all")
load("Data/MACAv2-METDATA/Tmax_historical_all")
load("Data/MACAv2-METDATA/Tmin_historical_all")
load("Data/MACAv2-METDATA/Pr_historical_all")
#met observation 1979-2016, 38 tot
load("Data/MET/GDD_all")
load("Data/MET/EDD_all")
load("Data/MET/Tmax_all")
load("Data/MET/Tmin_all")
load("Data/MET/Pr_all")

#take 26 yrs: 1979-2004 2074-2099
idx_proj<-c(9:34)
idx_hist<-c(10:35)
idx_obs<-c(1:26)

EDD_proj<-EDD_proj[ ,idx_proj]
EDD_hist<-EDD_hist[ ,idx_hist]
EDD_obs<-EDD_obs[ ,idx_obs]

GDD_proj<-GDD_proj[ ,idx_proj]
GDD_hist<-GDD_hist[ ,idx_hist]
GDD_obs<-GDD_obs[ ,idx_obs]

Pr_proj<-Pr_proj[ ,idx_proj]
Pr_hist<-Pr_hist[ ,idx_hist]
Pr_obs<-Pr_obs[ ,idx_obs]

Tmax_proj<-Tmax_proj[ ,idx_proj]
Tmax_hist<-Tmax_hist[ ,idx_hist]
Tmax_obs<-Tmax_obs[ ,idx_obs]

Tmin_proj<-Tmin_proj[ ,idx_proj]
Tmin_hist<-Tmin_hist[ ,idx_hist]
Tmin_obs<-Tmin_obs[ ,idx_obs]

boxhist3(EDD_obs,EDD_hist,EDD_proj,"EDD (degree days)","METDATA EDD observation (1979-2004)","MACA EDD hindcast (1979-2004)","MACA EDD projection (2074-2099)")
boxhist3(GDD_obs,GDD_hist,GDD_proj,"GDD (degree days)","METDATA GDD observation (1979-2004)","MACA GDD hindcast (1979-2004)","MACA GDD projection (2074-2099)")
boxhist3(Pr_obs,Pr_hist,Pr_proj,"Precipitation (mm)","METDATA Pr observation (1979-2004)","MACA Pr hindcast (1979-2004)","MACA Pr projection (2074-2099)")
boxhist3(Tmax_obs,Tmax_hist,Tmax_proj,"Tmax (degree C)","METDATA Tmax observation (1979-2004)","MACA Tmax hindcast (1979-2004)","MACA Tmax projection (2074-2099)")
boxhist3(Tmin_obs,Tmin_hist,Tmin_proj,"Tmin (degree C)","METDATA Tmin observation (1979-2004)","MACA Tmin hindcast (1979-2004)","MACA Tmin projection (2074-2099)")

Temp_proj<-(Tmax_proj+Tmin_proj)/2
Temp_hist<-(Tmax_hist+Tmin_hist)/2
Temp_obs<-(Tmax_obs+Tmin_obs)/2
boxhist2(Temp_hist-mean(Temp_hist),Temp_proj-mean(Temp_proj),"Mean temperature anomaly (degree C)","MACA temperature hindcast (1979-2004)","MACA temperature projection (2074-2099)")
cdf2(Temp_hist-mean(Temp_hist),Temp_proj-mean(Temp_proj),"Mean temperature anomaly (degree C)","MACA temperature hindcast (1979-2004)","MACA temperature projection (2074-2099)")
boxhist2(Pr_hist-mean(Pr_hist),Pr_proj-mean(Pr_proj),"Precipitation anomaly (mm)","MACA precipitation hindcast (1979-2004)","MACA precipitation projection (2074-2099)")
cdf2(Pr_hist-mean(Pr_hist),Pr_proj-mean(Pr_proj),"Precipitation anomaly (mm)","MACA Pr hindcast (1979-2004)","MACA Pr projection (2074-2099)")

ks.test(Temp_hist-mean(Temp_hist),Temp_proj-mean(Temp_proj))
ks.test(Pr_hist-mean(Pr_hist),Pr_proj-mean(Pr_proj))


#map plot of these mean variables
load("ANSI")
change<-data.frame(fips=ANSI,GDD=rowMeans(GDD_proj)-rowMeans(GDD_hist),EDD=rowMeans(EDD_proj)-rowMeans(EDD_hist),Pr=rowMeans(Pr_proj)-rowMeans(Pr_hist), 
                   Tmax=rowMeans(Tmax_proj)-rowMeans(Tmax_hist),Tmin=rowMeans(Tmin_proj)-rowMeans(Tmin_hist),Tmean=rowMeans(Temp_proj)-rowMeans(Temp_hist))
a=plot_usmap(regions = "counties",include=c("PA","NY","NJ","MD","DE","DC","NC","VA","SC","WV","OH","MI","GA","KY","IN","IL","AL","TN","WI","MS","MN","MO","LA","AR","IA"),
             data=change, values = "Tmean") + labs(title = "Change of mean temperature between 2074-2099 and 1979-2004")+ 
  scale_fill_gradient2(low = "white", high ="red", mid="lightpink",limits=c(2.5,4.5),midpoint=3.5,name="T (degree C)")+theme(plot.title = element_text(size=14))
plot(a)







#test for EDD
load("Data/MACAv2-METDATA/EDD_proj_linearshift")
boxhist2(EDD,EDD_proj,"EDD (degree ays)","EDD projection under linear shifted temperature (2074-2099)","MACA-METDATA EDD projection (2074-2099)")


