#plot of boxplot,pdf,cdf,survival function
combinedplot<-function(vector_data,Main_title,variable_name){ #input should be a vector
  layout(mat = matrix(c(1,2,3,4),4,1, byrow=TRUE),  height = c(1,3,3,3))
  vector_data<-vector_data[!is.na(vector_data)]
  Max = max(vector_data, na.rm = TRUE)
  Min = min(vector_data, na.rm = TRUE)
  par(mar=c(0, 5.1, 2.1, 2.1))
  boxplot(vector_data, horizontal = TRUE, xaxt="n", frame=F, pch=20,ylim=c(Min,Max), main=Main_title, cex.main=2.0) #boxplot
  par(mar=c(4, 5.1, 0, 2.1))
  plot(density(vector_data),xlim=c(Min,Max),xaxt="n",yaxt="n",xlab = "", ylab = "Density",main = "",cex.lab=2) #pdf
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  par(mar=c(4, 5.1, 1.6, 2.1))
  plot(ecdf(vector_data),xlim=c(Min,Max),xaxt="n",yaxt="n",xlab = "",ylab = "CDF",main = "",cex.lab=2)
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  par(mar=c(4,5.1,1.6,2.1))
  e<-ecdf(vector_data)
  plot(vector_data,1-e(vector_data),xlim=c(Min,Max),xaxt="n",yaxt="n",xlab = "",ylab="Survival function",main="",cex.lab=2,pch=20)#survival function
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
}

eccdf <- function(x) 
{ return( 1 - ecdf(x)(x) ) } 

combinedplot_2<-function(vector_data1,vector_data2,Main_title,variable_name1,variable_name2){ # 2 variables case
  layout(mat = matrix(c(1,2,3,4,5),5,1, byrow=TRUE),  height = c(1,1,3,3,3))
  vector_data1<-vector_data1[!is.na(vector_data1)]
  vector_data2<-vector_data2[!is.na(vector_data2)]
  #for EDD, get only positive values
  #vector_data1<-vector_data1[which(vector_data1>0)]
  #vector_data2<-vector_data2[which(vector_data2>0)]
  Max = max(max(vector_data1, na.rm = TRUE),max(vector_data2,na.rm=TRUE))
  #Max=0.2
  Min = min(min(vector_data1, na.rm = TRUE),min(vector_data2,na.rm=TRUE))
  #print(c("max=",Max,"min=",Min))
  par(mar=c(0, 5.1, 2.1, 2.1))
  boxplot(vector_data1, horizontal = TRUE, xaxt="n", frame=F, pch=20,ylim=c(Min,Max), main=Main_title, cex.main=2.0) #boxplot
  points(mean(vector_data1),1,col="blue",pch=16)
  par(mar=c(0, 5.1 ,2.1, 2.1))
  boxplot(vector_data2, horizontal = TRUE, xaxt="n", frame=F, pch=20,ylim=c(Min,Max),col="red")
  points(mean(vector_data2),1,col="blue",pch=16)
  print(mean(vector_data1))
  print(mean(vector_data2))
  par(mar=c(4, 5.1, 0, 2.1))
  plot(density(vector_data1),xlim=c(Min,Max),ylim=c(0,0.2),xaxt="n",yaxt="n",xlab = "", ylab = "Density",main = "",cex.lab=2) #pdf
  lines(density(vector_data2),xlim=c(Min,Max),xlab = "", ylab = "",main = "",cex.lab=2,col="red")
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  legend("topright",lty=c(1,1),col=c("black","red"),legend=c(variable_name1,variable_name2),cex=2)
  par(mar=c(4, 5.1, 1.6, 2.1))
  plot(ecdf(vector_data1),xlim=c(Min,Max),xaxt="n",yaxt="n",xlab = "",ylab = "CDF",main = "",cex.lab=2,col="white")
  e1<-ecdf(vector_data1)
  e2<-ecdf(vector_data2)
  curve(e1(x),add=TRUE)
  curve(e2(x),add=TRUE,col="red")
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  par(mar=c(4,5.1,1.6,2.1))
  plot(vector_data1,1-e1(vector_data1),xaxt="n",yaxt="n",xlim=c(Min,Max),ylim=c(0,1),xlab = "GDD",ylab="Survival function",main="",cex.lab=2,col="white")#survival function
  curve(1-e1(x),add=TRUE)
  curve(1-e2(x),add=TRUE,col="red")
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
}

combinedplot_3<-function(vector_data,Main_title,variable_name){ #historgram + boxplots
  layout(mat = matrix(c(1,2),4,1, byrow=TRUE),  height = c(6,1))
  vector_data<-vector_data[!is.na(vector_data)]
  Max = max(vector_data, na.rm = TRUE)
  Min = min(vector_data, na.rm = TRUE)
  par(mar=c(0, 5.1, 2.1, 2.1))
  boxplot(vector_data, horizontal = TRUE, xaxt="n", frame=F, pch=20,ylim=c(Min,Max), main=Main_title, cex.main=2.0) #boxplot
  par(mar=c(4, 5.1, 0, 2.1))
  hist(vector_data,xlim=c(Min,Max),xaxt="n",yaxt="n",xlab = variable_name , freq=FALSE, main = "",cex.lab=2) #histogram
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
}