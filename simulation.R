## This is the simulation data from KQR paper section 5.2 ##
source("./kqr.r")

## simulation data from (40)
n <- 200  # training
N1 <- 10000   # validation
N2 <- 10000   # test 
p <- 2

#tau <- 0.1
#tau <- 0.3
tau <- 0.5  # set tau = 50%

# 200 training data
tr.x <- matrix(runif(n*p, 0, 1), nrow = n, ncol = p)
tr.f <- 40 * exp(8*((tr.x[, 1]-0.5)^2 + (tr.x[, 2]-0.5)^2)) * (exp(8*((tr.x[, 1]-0.2)^2 + (tr.x[, 2]-0.7)^2)) + exp(8*((tr.x[, 1]-0.7)^2 + (tr.x[, 2]-0.2)^2)))^(-1)
## Normal
e.n <- rnorm(n, 0, 1)   # normal error
tr.y <- tr.f + e.n

# 10000 validation data
va.x <- matrix(runif(N1*p, 0, 1), nrow = N1, ncol = p)
va.f <- 40 * exp(8*((va.x[, 1]-0.5)^2 + (va.x[, 2]-0.5)^2)) * (exp(8*((va.x[, 1]-0.2)^2 + (va.x[, 2]-0.7)^2)) + exp((8*(va.x[, 1]-0.7)^2 + (va.x[, 2]-0.2)^2)))^(-1)
e.n <- rnorm(N1, 0, 1)   # normal error
va.y <- va.f + e.n

# 10000 test data
te.x <- matrix(runif(N2*p, 0, 1), nrow = N2, ncol = p)
te.f <- 40 * exp(8*((te.x[, 1]-0.5)^2 + (te.x[, 2]-0.5)^2)) * (exp(8*((te.x[, 1]-0.2)^2 + (te.x[, 2]-0.7)^2)) + exp(8*((te.x[, 1]-0.7)^2 + (te.x[, 2]-0.2)^2)))^(-1)
e.n <- rnorm(N2, 0, 1)  # normal error
te.y <- te.f + e.n

tr.k <- spline.kernel(tr.x, tr.x)
va.k <- spline.kernel(tr.x, va.x)
te.k <- spline.kernel(tr.x, te.x)

res <- kqrpath(tr.y, tau, K=tr.k)

## SIC
ind <- order(res$Csic)[1]
fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
sic.n.df[j] <- res$df[ind]
sic.n.mcd[j] <- pf(fit, te.f.n, tau) / N2