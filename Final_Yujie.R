# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.


lb = -Inf
ub = Inf
n= 100
k=2

fx <- function(x){ dnorm(x, mean=1, sd=0.5)}
#test functions 
ars(g=fx, lb=-100, ub=100, n=100, xvec=c(-2,3)) #runs forever!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ahhhh
ars(g=dnorm, lb=-100, ub=100, n=100, xvec=c(-2,3)) # the plot looks weird


ars <-function(g, lb=-100, ub=100, n, xvec=c(-1,1)){
# input function is g(x), lower bound and upper bound, k is the number of initial points
# n is the number of samples to generate from g(x)


# initialize the abscissae in Tk.

k <- length(xvec)
x1 <- xvec[1]
xk <- xvec[k]        #x1 + 2*(k-1)
Tk <- sort(xvec)

hk <<- h(Tk, g)
if (min(hk) == -Inf) stop("log(g(x)) returns -Inf, try other initial values")
if (k < 2) stop("Must have at least 2 inital values")
if (min(Tk) < lb | max(Tk) > ub) stop("Initial values have to be inside the bounds.")
dhk <<- hPrime(Tk, g)

if ((ub == Inf & dhk[k] > 0) | (lb == -Inf & dhk[1] < 0)) stop("Initial values should be on both sides of the mode")



size <- 0
## Keep sampling until we accept n values
samples <- NULL


while(size < n){
  hk <<- h(Tk, g)
  dhk <<- hPrime(Tk, g)
  z <<- z_generate(k,Tk, lb, ub)
  w = runif(1)
  xStar <- drawxStar(z,hk,dhk,Tk)
  if (w <= exp(lk(xStar, Tk)-uk(xStar, Tk))){
    size = size+1
    samples[size] = xStar
  }
  else{
    hxStar <- h(xStar, g)
    hPxStar <- hPrime(xStar, g)
    TkOld <- Tk # for debugging purpose
    Tk <<- sort(c(Tk, xStar))
    k <<- k + 1
    if(w <= exp(hxStar-uk(xStar, Tk))){
      size <<- size + 1
      samples[size] = xStar
    }
  }
}

return(samples)
}


##Helper functions:

# Calculate h(x)
h <- function (x, g) {
  log((g(x))) }

# derivative function of h(x):
hPrime <- function(x, g){
  d = sqrt(.Machine$double.eps)
  return((h(x+d/2, g)-h(x-d/2, g))/d)
}

# z_generate is a function that returns the z vector based on equation 1.
z_generate <- function(k,Tk, lb, ub){
  z_j=(hk[-1] - hk[-k] - Tk[-1]*dhk[-1] + Tk[-k]*dhk[-k]) / (dhk[-k]-dhk[-1]) # equation 1 in the paper
  z=c(lb,z_j,ub)
  z=sort(z)
  return (z)
}

# compute uk(x)
uk <- function(x, Tk){
  #i = sapply(x, function(x){which(x<=z)[1]})
  j <- which(x<=z)[1]
  i <- j-1
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
  k <- length(Tk)
  j <- which(x<=Tk)[1]
  j <- j-1
  # if x<x[1] or x>x[k] define l=-Inf
  if ((x < Tk[1]) || (x > Tk[k])){
    return (-Inf)
  }
  # find j such that x in (x[j],x[j+1])
  
  return (((Tk[j+1]-x)*hk[j]+(x-Tk[j])*hk[j+1])/(Tk[j+1]-Tk[j]))
}


#Calculate y values of z under uk(x)


# Updated
area <- function(yZ){
  y1 <- yZ[1:(length(yZ)-1)]
  y2 <- yZ[2:length(yZ)]
  x1 <- z[1:(length(z)-1)]
  x2 <- z[2:length(z)]
  slope = (y2 - y1) / (x2 - x1)
  area = (exp(y2)-exp(y1))/slope
  return(area)
}


#area <- function(k, Tk){
#  tmp <- c(Tk[1], z)
#  areas <- (exp(hk)-dhk*hk)/dhk)*(exp(dhk*tmp[-1])-exp(dhk*tmp[-length(tmp)]))
#  return(areas)
#}


drawxStar <- function(z, hk, dhk, Tk){
  k <- length(Tk)
  yZ <- hk - (Tk - z[-(k+1)])*dhk
  yZ[k+1] <- hk[k] + (z[k+1] - Tk[k])*dhk[k]
  areas <- area(yZ)
  ## Calculate xStar, out sample candidate
  sumArea<- sum(areas)
  ## sumArea is the total area under s(x), and pct is the percentage of each piece of area over the total area
  pct = areas/sumArea
  cumsumCDF <- cumsum(pct)
  #cumsumCDF gives the CDF values under s(x)
  u = runif(1)
  j = which(u <= cumsumCDF)[1]
  x1 = z[j]
  x2 = z[j+1]
  y1 = yZ[j]
  y2 = yZ[j+1]
  slope = (y2-y1)/(x2-x1)
  xStar = 1/slope*log(u*exp(y2)+(1-u)*exp(y1))-y1/slope+x1
    #if(j > 1){
  #  hpj <- hPrime(Tk[j])
  #  xStar <- 1/hpj*log(exp(hpj*z[j-1]) + sumArea*(u-cumsumCDF[j-1])*hpj*exp(hpj*Tk[j])/exp(h(Tk[1])))
  #}else{
  #  hp1 <- hPrime(Tk[1])
  #  xStar <- 1/hp1*log(exp(hp1*lb) + sumArea*u*hp1*exp(hp1*Tk[1])/exp(hp1))
  #}
  return(xStar)
}


####################################################################

# Log-concave checks:
# Some densities will not work with this kind of sampling; we need to check for those.

  # Works: Normal, Exponential, Gamma(r>=1), Beta
  # Doesn't work: Gamma(r<1), Pareto, t, F, chi-squared

  # We will also have to check that the bounds are appropriate for the given density.
