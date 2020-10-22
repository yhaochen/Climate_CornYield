#Sampling parametric uncertainty by Latin Hypercube Sampling process

#function input: model, sample number, sample range 
library(lhs) #Latin Hypercube Sampling
parasample<-function(model,n,range){
  paranum<-length(model$coefficients) #number of parameters in the model
  paramean<-rep(NA,paranum) #parameter means
  parastd<-rep(NA,paranum)  #parameter standard errors
  output<-summary(model) #summary of model, including the parameter means/stds
  for (i in 1:paranum){
    paramean[i]<-output$coefficients[i,1]
    parastd[i]<-output$coefficients[i,2]
  }
  lhs<-randomLHS(n,paranum)
  for (i in 1:n){
    for (j in 1:paranum){
      samplingrange<-range[which(colnames(range)==names(model$coefficients)[j])]
      lhs[i,j]<-paramean[j]-sig*parastd[j]+lhs[i,j]*2*sig*parastd[j]
    }
  }
  lhs
}

parasample_log<-function(model,n,sig){ #logyield model, with linear and quadratic time in the model, do not sampling them
  paranum<-length(model$coefficients)-2 #number of parameters in the model
  paramean<-rep(NA,paranum) #parameter means
  parastd<-rep(NA,paranum)  #parameter standard errors
  output<-summary(model) #summary of model, including the parameter means/stds
  for (i in 1:paranum){
    paramean[i]<-output$coefficients[i,1]
    parastd[i]<-output$coefficients[i,2]
  }
  lhs<-randomLHS(n,paranum)
  for (i in 1:n){
    for (j in 1:paranum){
      lhs[i,j]<-paramean[j]-sig*parastd[j]+lhs[i,j]*2*sig*parastd[j]
    }
  }
  lhs
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



frac_between <- function (vect1,vect2,vect3){ #vect2 is lower, vect3 is higher
  totlen<-length(vect1)
  out<-0
  for (i in 1:totlen){
    if ((vect1[i]<vect2[i]) | (vect1[i]>vect3[i])){
      out<-out+1
    }
  }
  return((totlen-out)/totlen)
}