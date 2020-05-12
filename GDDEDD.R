# Input: two vectors of the same length: Tmax, Tmin in celcius degree
# Output: total Growing degree days and extreme degree days
#Now add a "frost degree days?"
GDDEDD <- function(Tmax,Tmin) {
  if (length(Tmax)!=length(Tmin)){
    print("Tmax and Tmin must have the same length")
    exit
  }
  n<-length(Tmax)
  Tmax_1=Tmax
  Tmin_1=Tmin
  for (i in 1:n){
    if (Tmax[i]>29){
      Tmax[i]=29
    }
    if (Tmax[i]<10){
      Tmax[i]=10
    }
    if (Tmin[i]>29){
      Tmin[i]=29
    }
    if (Tmin[i]<10){
      Tmin[i]=10
    }
  }
  GDD<-sum((Tmax+Tmin)/2-10)
  for (i in 1:n){
    if (Tmax_1[i]<29){
      Tmax_1[i]=29
    }
    if (Tmin_1[i]<29){
      Tmin_1[i]=29
    }
  }
  EDD<-sum((Tmax_1+Tmin_1)/2-29)
  c(GDD,EDD)
}