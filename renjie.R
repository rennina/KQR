source("source.R") # R code for generating the FHT data
library(kernlab)
nobs = 100
p = 20
rho = 0.1

tau = 0.5
x=genx2(nobs,p,rho)
y=genjerry(x,3)

# z_sigma = 10^(-1.2)
ulam <- 0.05
sigma = 0.01
## initialize kernel function
kern <- rbfdot(sigma=sigma)
K = kernelMatrix(kern, x)

