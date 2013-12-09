# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.



simulation <-function(g, D = c(-Inf,Inf), n, k = 10){
# input function is g(x) and D, where D the domain D=(lower,upper)
# n is the number of samples to generate from g(x)
h <- function (x) log((g(x)))  # g will be input density


hPrime <- function(x){
  d = 10^-10   # machine epsilon??
  return((h(x+d)-h(x))/d)
} # derivative function of h

z = rep (0,k+1) # initialize a vector to store values for z.

# initialize the abscissae in Tk.
if (D[1]==-Inf){
  # The following few lines choose x1 such that hPrime(x1)>0, not sure if this is the best way to do it.
  # Start x1 at large negative number, then add 1.
  x1 = 10^-8
  while (hPrime(x1) <= 0){
    x1 = x1 + 1
  }
  if(D[2]== Inf){
    xk = 10^8
    while (hPrime(xk) >= 0){
    xk = xk - 1
  }
  }
  else{
    xk = D[2]
  }
}
else{
  x1 = D[1]
  if(D[2]== Inf){
    xk = 10^8
    while (hPrime(xk) >= 0){
    xk = xk - 1
  }else{
    xk = D[2]
  }
}
}
Tk <- seq(x1,xk,length=k)


# z_generate is a function that returns the z vector based on equation 1.
z_generate <- function(Tk){

 for (j in 2:length(Tk)){
  z[j] = (h(Tk[j])-h(Tk[j-1])-Tk[j]*hPrime(Tk[j])+Tk[j-1]*hPrime(Tk[j-1]))/(hPrime(Tk[j-1])-hPrime(Tk[j])) # x starts from 1 and z starts from 0 so need to adjust index.
  z[1] = D[1]
  z[k+1] = D[2]
 }
return (z)
}

z=z_generate(Tk)
# compute u(x)
u_generate <- function(x){
  i = 1
  while (z[i]<x){
    i = i+1
  }
  return (h(Tk[i-1])+(x-Tk[i-1])*hPrime(Tk[i-1]))
}


library(stats)

s <- function(x,u) {
  integrand = function(y) exp(u_generate(y))
  return(integrand(x)/integrate(integrand, lower=D[1], upper=D[2]))
}

l <- function(x,Tk) {
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
 if (w<=exp(l(xStar)-u(xStar))){
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