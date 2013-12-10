# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.



simulation <-function(g, lb= -Inf, ub= Inf, n){
# input function is g(x) and D, where D the domain D=(lower,upper)
# n is the number of samples to generate from g(x)
k = 2  # two starting abscissae

h <- function (x) log((g(x)))  # g will be input density
hPrime <- function(x){
  d = sqrt(.Machine$double.eps)
  return((h(x+d/2)-h(x-d/2))/d)
} # derivative function of h

z_j = rep(NA,k-1) # initialize a vector to store values for z.

# initialize the abscissae in Tk.
if (lb==-Inf & ub==Inf){
  x1 = 0
}
if( lb==-Inf & ub!=Inf){
    x1 = ub-1
}
if (lb!=-Inf & ub==Inf){
  x1 = lb+1
}
if( lb!=-Inf & ub!=Inf){
  x1 = (lb+ub)/2
}

xk = x1 + 1*(k-1)

Tk <- seq(x1,xk,length=k)

# z_generate is a function that returns the z vector based on equation 1.
z_generate <- function(Tk){
 for (j in 1:(length(Tk)-1)){
  xj1=Tk[j+1]
  xj=Tk[j]
  z_j[j] = (h(xj1)-h(xj)-xj1*hPrime(xj1)+xj*hPrime(xj))/(hPrime(xj)-hPrime(xj1)) # x starts from 1 and z starts from 0 so need to adjust index.
  z=c(lb,z_j,ub)
 }
 return (z)
}

z=z_generate(Tk)

# compute uk(x)
uk <- function(x){
  i = 1
  while (z[i]<x){
    i = i+1
  }
  return (h(Tk[i-1])+(x-Tk[i-1])*hPrime(Tk[i-1]))
}


library(stats)

sk <- function(x) { #not sure if we need the input u
  integrand = function(y) exp(uk(y))
  return(integrand(x)/integrate(integrand, lower=lb, upper=ub))
}

###Need to Modify#################################################
#
yZ <- rep(NA, k+1) # corresponding y values for Tk
yZ[1] <- hPrime(TK[1])*z[1] + h(Tk[1]) - hPrime(Tk[1])*Tk[1]
for (i in 2:k) {
  slope = hPrime[i]
  x = Tk[i]
  y = h(x)
  yZ[i+1] = slope*z[i+1] + (y - slope*x)
}

area <- rep(NA, k)
for (i in 1:length(Tk)) {
  x = z[i]
  y = yZ[i]
  x2 = z[i+1]
  y2 = yZ[i+1]
  slope = (y2-y) / (x2-x1)
  int = y - slope*x
  
  area[i] = exp(y2)/slope - exp(y)/slope
}

sumArea<- sum(area)
## Let sumA be the total area under u(x) and ratioA be the cumulative percent area of each piece under u(x)
pct = area/sumArea
cumsumArea <- cumsum(pct)

## Keep sampling until we accept n values
size = 0
while (size < n) {
  ## Calculate xStar, out sample candidate
  u = runif(1)
  j = which(pct <= cumsumArea)[1]
  if(j > 1){
    hpj <- hPrime(Tk[j])
    xStar <- 1/hpj*log(exp(hpj*z[j-1]) + sumArea*(u-cumsumArea[j-1])*hpj*exp(hpj*Tk[j])/exp(h(Tk[1])))  
  }else{
    hp1 <- hPrime(Tk[1])
    xStar <- 1/hp1*log(exp(hp1*lb) + sumArea*u*hp1*exp(hp1*Tk[1])/exp(hp1))
  }
  
#############################################################################################
  
  ##???? Calculate xStar evaluated at u(x), uxStar
  um = (hz[indexA+1] - hz[indexA]) / (z[indexA+1] - z[indexA])
  ub = hz[indexA] - um*z[indexA]
  uxStar = um*xStar + ub
  
  ##???? Calculate xStar evaluated at l(x), lxStar
  lxStar = 0
  indexTk = which(xStar<Tk)[1]
  if (indexTk == 1 | length(indexTx) == 0) lxStar = 0
  else {
    lm = (hx(Tk[indexTk]) - hx(Tk[indexTk-1])) / (Tk[indexTk] - Tk[indexTk-1])
    lb = hx[indexTk] - lm*Tk[indexTk]
    lxStar = lm*xStar + lb
  }
  
  ## Automatically accept if w <= exp(lxStar - uxStar)
  w = runif(1)
  LUratio = exp(lxStar - uxStar)
  if (w <= LUratio) {
    samps[size+1] = xStar
    size = size+1
  }
  ## Otherwise, evaluate h(x), update our vectors, and check to see if we accept xStar
  else {
    hxStar = hx(xStar)
    Tk = sort(c(Tk, xStar))
    xIndex = which(Tk == xStar)
    if (xIndex == 1) {
      hk = c(hx(xStar), hk)
      dhk = c(dhx(xStar), dhk)
    }
    else if (xIndex == length(Tk)) {
      hk = c(hk, hx(xStar))
      dhk = c(dhk, dhx(xStar))
    }
    else {
      hk = c(hk[1:(xIndex-1)], hx(xStar), hk[xIndex:length(hk)])
      dhk = c(dhk[1:(xIndex-1)], dhx(xStar), dhk[xIndex:length(dhk)])
    }
    zNew1 = (hk[xIndex] - hk[xIndex] - Tk[i+1]*dhk[i+1] + Tk[i]*dhk[i]) / (dhk[i] - dhk[i+1])
    zNew2 =
      z[xIndex] = znew1
    z = sort(c(z, znew2)) ####
    
    HUratio = exp(hxStar - uxStar)
    if (w <= HUratio) {
      samps[size+1] = xStar
      size = size+1
    }
  }
}
}





#####################################################################
lk <- function(x,Tk) {
  # if  x<x[1] or x>x[k] define l=-Inf
  if ((x < Tk[1]) || (x > Tk[length(x)])){
    return (-Inf)
  }
  j=1
  # find j such that x in (x[j],x[j+1])
  while (Tk[j+1] < x){ 
    j = j+1
  }
  return (((Tk[j+1]-x)*h(Tk[j])+(x-Tk[j])*h(Tk[j+1]))/(Tk[j+1]-Tk[j]))
}


# Sampling step and updating step

output=c() # store accepted points

while (length(output)<n){


 # generate xStar from s(x) (not sure how, maybe inversCDF(uniform)?)
 # could use sample():
  
 w = runif(1)
 
 # squeezing test
 if (w<=exp(lk(xStar)-uk(xStar))){
   output <- append(output, xStar)
 }
 # updating below
 else{
   T <- sort(append(T,xStar)) #include xStar in T
   # rejection test
    if (w<=exp(h(xStar)-u(xStar))){
      output <- append(output,xStar)
    }
   z <- z_generate(T)
   u <- function(x,z,T)
   s <- function(x,u)
   l <- function(x,T)
  }
}

return(output)
}


# Log-concave checks:
# Some densities will not work with this kind of sampling; we need to check for those.

  # Works: Normal, Exponential, Gamma(r>=1), Beta
  # Doesn't work: Gamma(r<1), Pareto, t, F, chi-squared

  # We will also have to check that the bounds are appropriate for the given density.
