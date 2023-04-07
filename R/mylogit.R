mylogit <- function(x){
  x[x==0] = 0.0001
  x[x==1] = 0.999
  log(x/(1-x))
}
