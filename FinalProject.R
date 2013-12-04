# For consistency, we can use camel-case throughout the project to name variables, functions, etc. This is
# when the word begins lower-case, and any other words added on are upper-cased, e.g. meanFunction, xStar.

# Formula calculations:
# These will be inside of a function and most likely a for loop; we need to figure out how to 
# calculate them, given inputs.

upperBound
lowerBound
x = 1:k  # k will be input number of values
j = 1:(k-1)

h = log(f)  # f will be input density
hPrime

u = h[x] + (x-x[j])*hPrime[x]

s = exp{u[x]}/integrate(u[x], lowerBound, upperBound)

l = ((x[j+1] - x)*h[x[j]] + (x - x[j])*h[x[j+1]])/(x[j+1]-x[j])


# Log-concave checks:
# Some densities will not work with this kind of sampling; we need to check for those.

  # Works: Normal, Exponential, Gamma(r>=1), Beta
  # Doesn't work: Gamma(r<1), Pareto, t, F, chi-squared

  # We will also have to check that the bounds are appropriate for the given density.