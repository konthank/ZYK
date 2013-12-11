# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.


simulation <-function(g, lb=-100, ub=100, n, k=2){
# input function is g(x) and D, where D the domain D=(lower,upper)
# n is the number of samples to generate from g(x)

h <- function (x) log((g(x))) # g will be input density

 # derivative function of h:
hPrime <- function(x){
  d = sqrt(.Machine$double.eps)
  return((h(x+d/2)-h(x-d/2))/d)
}

z_j = rep(NA,k-1) # initialize a vector to store values for z.

# initialize the abscissae in Tk.
x1 = -1
xk = x1 + 2*(k-1)
Tk <- seq(x1,xk,length=k)

# z_generate is a function that returns the z vector based on equation 1.
z_generate <- function(Tk){
 for (j in 1:(k-1)){
  xj1=Tk[j+1]
  xj=Tk[j]
  z_j[j] = (h(xj1)-h(xj)-xj1*hPrime(xj1)+xj*hPrime(xj))/(hPrime(xj)-hPrime(xj1)) # x starts from 1 and z starts from 0 so need to adjust index.
  z=c(z_j,100)
 }
 return (z)
}

z=z_generate(Tk)

# compute uk(x)
uk <- function(x, Tk){
  i = 1
  while (z[i]<x){
    i = i+1
  }
  return (h(Tk[i-1])+(x-Tk[i-1])*hPrime(Tk[i-1]))
}

library(stats)

# compute sk(x)
sk <- function(x) { #not sure if we need the input u
  integrand = function(y) exp(uk(y))
  return(integrand(x)/integrate(integrand, lower=lb, upper=ub))
}

# compute lk(x):
lk <- function(x,Tk) {
  # if x<x[1] or x>x[k] define l=-Inf
  if ((x < Tk[1]) || (x > Tk[k])){
    return (-Inf)
  }
  j=1
  # find j such that x in (x[j],x[j+1])
  while (Tk[j+1] < x){
    j = j+1
  }
  return (((Tk[j+1]-x)*h(Tk[j])+(x-Tk[j])*h(Tk[j+1]))/(Tk[j+1]-Tk[j]))
}


yZ <- function(k,Tk){
  yZ <- rep(NA, k+1) # corresponding y values for Tk
  yZ[1] <- hPrime(Tk[1])*z[1] + h(Tk[1]) - hPrime(Tk[1])*Tk[1]
  for (i in 1:k) {
    slope = hPrime(Tk[i])
    x = Tk[i]
    y = h(x)
    yZ[i+1] = slope*z[i+1] + (y - slope*x)
  }
  return(yZ)
}

#area <- function(k, Tk){
#  yZ <-yZ(k, Tk)
#  area <- rep(NA, k)
#  for (i in 1:length(z)-1) {
#   x = Tk[i]
#   y = yZ[i]
#   x2 = z[i+1]
#   y2 = yZ[i+1]
#   slope = (y2-y) / (x2-x)
#   int = y - slope*x
#  
#   area[i] = exp(y2)/slope - exp(y)/slope
#  }
#  return(area)
# }

area <- function(k, Tk){
  tmp <- c(Tk[1], z)
  areas <- (exp(h(Tk)-hPrime(Tk)*Tk)/hPrime(Tk))*(exp(hPrime(Tk)*tmp[2:length(tmp)])-exp(hPrime(Tk)*tmp[1:(length(tmp)-1)]))
  return(areas)
}


## Keep sampling until we accept n values
size = 0
drawxStar <- function(k, Tk){
  ## Calculate xStar, out sample candidate
  sumArea<- sum(area(k, Tk))
  ## Let sumA be the total area under u(x) and ratioA be the cumulative percent area of each piece under u(x)
  pct = area(k, Tk)/sumArea
  cumsumArea <- cumsum(pct)
  u = runif(1)
  j = which(pct <= cumsumArea)[1]
  if(j > 1){
    hpj <- hPrime(Tk[j])
    xStar <- 1/hpj*log(exp(hpj*z[j-1]) + sumArea*(u-cumsumArea[j-1])*hpj*exp(hpj*Tk[j])/exp(h(Tk[1])))
  }else{
    hp1 <- hPrime(Tk[1])
    xStar <- 1/hp1*log(exp(hp1*lb) + sumArea*u*hp1*exp(hp1*Tk[1])/exp(hp1))
  }
  return(xStar)
}
#############################################################################################

samples = NULL
while(size < n){
  w = runif(1)
  xStar <- drawxStar(k, Tk)
  if (w <= exp(lk(xStar, Tk))-uk(xStar, Tk)){
    size = size+1
    samples[size] = xStar
  }
  else{
    hxStar = h(xStar)
    hPxStar = hPrime(xStar)
    Tk = sort(c(Tk, xStar))
    k = k + 1
    z = z_generate(Tk)
    if(w <= exp(hxStar-uk(xStar, Tk))){
      size = size + 1
      samples[size] = xStar
    }
  }
}





####################################################################



# Log-concave checks:
# Some densities will not work with this kind of sampling; we need to check for those.

  # Works: Normal, Exponential, Gamma(r>=1), Beta
  # Doesn't work: Gamma(r<1), Pareto, t, F, chi-squared

  # We will also have to check that the bounds are appropriate for the given density.
