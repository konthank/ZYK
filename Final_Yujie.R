# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.


lb = -Inf
ub = Inf
n= 100
k=2

ars(g=dnorm, lb=-100, ub=100, n=100, k=2)

ars <-function(g, lb=-100, ub=100, n, k=2){
# input function is g(x), lower bound and upper bound, k is the number of initial points
# n is the number of samples to generate from g(x)
 
  
# initialize the abscissae in Tk.
x1 = -1
xk = 1        #x1 + 2*(k-1)

Tk <- seq(x1,xk,length=k)

hk <- h(Tk)

if (min(hk) == -Inf) stop("log(g(x)) returns -Inf, try other initial values")

dhk <- hPrime(Tk)

size = 0
## Keep sampling until we accept n values
samples = NULL
z=z_generate(k,Tk, lb, ub)

while(size < n){
  w = runif(1)

  areas <- area(k, Tk)
  xStar <- drawxStar(k, Tk)
  if (w <= exp(lk(xStar, Tk))-uk(xStar, Tk)){
    size = size+1
    samples[size] = xStar
  }
  else{
    hxStar = h(xStar)
    hPxStar = hPrime(xStar)
    TkOld = Tk # for debugging purpose
    Tk = sort(c(Tk, xStar))
    k = k + 1
    if(w <= exp(hxStar-uk(xStar, Tk))){
      size = size + 1
      samples[size] = xStar
    }
  }
}


}


##Helper functions:

# Calculate h(x)
h <- function (x) {
  log((g(x))) }

# derivative function of h(x):
hPrime <- function(x){
  d = sqrt(.Machine$double.eps)
  return((h(x+d/2)-h(x-d/2))/d)
}

# z_generate is a function that returns the z vector based on equation 1.
z_generate <- function(k,Tk, lb, ub){
  z_j = rep(NA,k-1) # initialize a vector to store values for z.
  z_j=(hk[-1] - hk[-k] - Tk[-1]*dhk[-1] + Tk[-k]*dhk[-k]) / (dhk[-k]-dhk[-1]) # equation 1 in the paper
  z=c(z_j,ub)
  return (z)
}

# compute uk(x)
uk <- function(x, Tk){
  i = sapply(x, function(x){which(x<=z)[1]})
  return (hk[i]+(x-Tk[i])*dhk[i])
}

# compute sk(x)
sk <- function(x, Tk, lb, ub) {
  library(stats)
  uk <- uk(x, Tk)
  fx <- function() exp(uk)
  return(exp(uk)/integrate(fx(x), lower=lb, upper=ub))
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
  z <- z_generate(k, Tk, lb, ub)
  yZ <- rep(NA, k+1) # corresponding y values for Tk
  yZ[1] <- hPrime(Tk[1])*z[1] + h(Tk[1]) - hPrime(Tk[1])*Tk[1]
  for (i in 1:k) {
    slope = hPrime(Tk[i])
    x = Tk[i]
    y = h(x)
    yZ[i+1] = slope*z[i+1] + (y - slope*x)
  }
  return(yZ)
} ####Need to modify


area <- function(k, Tk){
  z <- z_generate(k, Tk, lb, ub)
  yZ <-yZ(k+1, Tk)
  area <- rep(NA, k)
  y1 <- yZ[1:(length(yZ)-1)]
  y2 <- yZ[2:length(yZ)]
  x1 <- z[1:(length(z)-1)]
  x2 <- z[2:length(z)]
  for (i in 1:k) {
    slope = (y2[i]- y1[i])/ (x2[i]-x1[i])
    int = y1[i] - slope*x1[i]
    area[i] = exp(y2[i])/slope - exp(y1[i])/slope
  }
  return(area)
}


area <- function(k, Tk){
  tmp <- c(Tk[1], z)
  areas <- (exp(hk)-dhk*hk)/dhk)*(exp(dhk*tmp[-1])-exp(dhk*tmp[-length(tmp)]))
  return(areas)
}


drawxStar <- function(k, Tk){
  ## Calculate xStar, out sample candidate
  sumArea<- sum(areas)
  ## sumArea is the total area under u(x), and pct is the percentage of each piece of area over the total area
  pct = areas/sumArea
  cumsumCDF <- cumsum(pct)
  #cumsumCDF gives the piece-wise CDF values under u(x)
  u = runif(1)
  j = which(pct <= cumsumCDF)[1]
  if(j > 1){
    hpj <- hPrime(Tk[j])
    xStar <- 1/hpj*log(exp(hpj*z[j-1]) + sumArea*(u-cumsumCDF[j-1])*hpj*exp(hpj*Tk[j])/exp(h(Tk[1])))
  }else{
    hp1 <- hPrime(Tk[1])
    xStar <- 1/hp1*log(exp(hp1*lb) + sumArea*u*hp1*exp(hp1*Tk[1])/exp(hp1))
  }
  return(xStar)
}


####################################################################

# Log-concave checks:
# Some densities will not work with this kind of sampling; we need to check for those.

  # Works: Normal, Exponential, Gamma(r>=1), Beta
  # Doesn't work: Gamma(r<1), Pareto, t, F, chi-squared

  # We will also have to check that the bounds are appropriate for the given density.
