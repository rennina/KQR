## six observations from the simple example in KQR paper ##
source("./kqr.r")

n <- 6
p <- 1

tau <- 0.5

tr.x <- matrix(runif(n*p, -2, 2), nrow=n, ncol=p)
tr.f <- sin(pi * tr.x[,1])/(pi * tr.x[,1])
e.n <- rnorm(n, 0, 0.2)
tr.y <- tr.f + e.n
tr.k <- spline.kernel(tr.x, tr.x)
  
res <- kqrpath(tr.y, tau, K=tr.k)
