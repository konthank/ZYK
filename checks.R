
warnings<-function(){
  
  if(k <= 0 || k%%1 != 0{
    stop("k must be non-negative, non-zero integer")
  }
  if(n <= 0 || n%%1 != 0){
    stop("n must be non-negative, non-zero integer")
  }
  if(g == dunif || g == dPareto || g == dt || g == dF || g == dchisq){
    stop("input density is not log-concave")
  }
  if(hPrime[1]<=0 || hPrime[2]>=0){
    stop{"Upper and Lower Bounds are not valid"}
  }

     
}
