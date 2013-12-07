# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.

# Formula calculations:
# These will be inside of a function and most likely a for loop; we need to figure out how to 
# calculate them, given inputs.

upperBound
lowerBound
x = 1:k  # k will be input number of values
j = 1:(k-1)

# input function is g(x) and D, where D the domain D=(lower,upper)
h = function (x) log((g(x)))  # g will be input density

deriv <- function(f,x){ 
  # computes derivative of a given function f at x
  h = 10^-10
  return((f(x+h)-f(x))/h)
}

hPrime <- function(x) deriv(h,x) # derivative function of h

z = rep (0,k+1) # initialize a vector to store values for z.

# initialize the abscissae in Tk.
k=10  # initial value of k, we may change this.
if (D[1]=-Inf){
  # The following few lines choose x1 such that hPrime(x1)>0, not sure if this is the best way to do it.
  x1 = 0
  while (hPrime(x1) <= 0){
    x1 = x1 - 1
  }
}
if (D[2]==Inf){
  xk = 0
  while (hPrime(xk) >= 0){
    xk = xk + 1
  }
}
else {
  x1 = D[1]
  xk = D[2]
}
T = seq(x1,xk,length=k)


for (j in 2:k){
  z[j] = (h(T[j])-h(T(j-1))-T[j]*hPrime(T[j])+T[j-1]*hPrime(T[j-1]))/(hPrime(T[j-1])-hPrime(T[j])) # x starts from 1 and z starts from 0 so need to adjust index.
  z[0] = D[1]
  z[k+1] = D[2]
}



u = h[x] + (x-x[j])*hPrime[x]

s = exp{u[x]}/integrate(u[x], lowerBound, upperBound)

l = ((x[j+1] - x)*h[x[j]] + (x - x[j])*h[x[j+1]])/(x[j+1]-x[j])


f<-function(x) x^2


# Log-concave checks:
# Some densities will not work with this kind of sampling; we need to check for those.

  # Works: Normal, Exponential, Gamma(r>=1), Beta
  # Doesn't work: Gamma(r<1), Pareto, t, F, chi-squared

  # We will also have to check that the bounds are appropriate for the given density.