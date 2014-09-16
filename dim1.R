source("./kqr.r")
## 1-Dimensional
n <- 200
N1 <- 10000
N2 <- 10000
p <- 1

## define kernel: radial basis
radial.kernel <- function(x, y = x, sigma = 0.2) {
  n <- nrow(x)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- -2 * a + normx + outer(rep(1, n), normy, "*")
  exp(- a / (2 * sigma^2))
}

tau <- 0.5

#### Validation Data
va.x <- matrix(runif(N1, 0, 1), nrow=N1, ncol=p)
va.f <- 2 * (exp( -30 * (va.x[,1]-0.25)^2 ) + sin(pi * va.x[,1]^2))
## Normal
e.n <- rnorm(N1, 0, 1)
va.y.n <- va.f + e.n

#### Test Quantiles
te.x <- matrix(runif(N2*p, 0, 1), nrow=N2, ncol=p)
te.f <- 2 * (exp( -30 * (te.x[,1]-0.25)^2 ) + sin(pi * te.x[,1]^2))
## Normal
q.n <- qnorm(tau, 0, 1)
te.f.n <- te.f + q.n

############# Method in KQR paper ##############
#### Degrees of Freedom
sic.n.df <- NULL
gacv.n.df <- NULL
gold.n.df <- NULL

## Mean Check Deviation
sic.n.mcd <- NULL
gacv.n.mcd <- NULL
gold.n.mcd <- NULL

### Training data
  tr.x <- matrix(runif(n*p, 0, 1), nrow=n, ncol=p)
  tr.k <- radial.kernel(tr.x, tr.x)
  va.k <- radial.kernel(tr.x, va.x)
  te.k <- radial.kernel(tr.x, te.x)
  tr.f <- 2 * (exp( -30 * (tr.x[,1]-0.25)^2 ) + sin(pi * tr.x[,1]^2))
  e.n <- rnorm(n, 0, 1)
  tr.y <- tr.f + e.n

  res <- kqrpath(tr.y, tau, K=tr.k)

j <- 1

  ## Gold Standard
  va.g <- (res$alpha0.g0 + t(t(va.k) %*% res$alpha.g)) / res$lambda
  err.g <- apply(va.g, 1, pf, va.y.n, tau) / N1
  ind <- order(err.g)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gold.n.df[j] <- res$df[ind]
  gold.n.mcd[j] <- pf(fit, te.f.n, tau) / N2
##################################################
lambda <- res$lambda[ind] ### The gold standard lambda chosen by KQR

## Using the kqr function (but will not allow adding weights):

library(kernlab)
r.kqr <- kqr(tr.x, tr.y, C = res$lambda[ind], kern = "rbfdot", kpar = list(sigma = 1/(2*0.2^2)), scale = FALSE)

## Using the ipop function (the chosen function to use with prefixed lambda):

H = tr.k
c = - tr.y
A = rep(1, n)
b = 0
r = 0
l = matrix(lambda * (tau - 1), n, 1)
u = matrix(lambda * tau, n, 1)
 
r.ipop <- ipop(c, H, A, b, l, u, r)


###### Weighted version #########

## Centered data:
percent <- 0.05 ## 5% of the data points are censored
cen <- sample(1:n, n*percent)  ## Indices of censored
tr.y[cen] <- 0
tr.d <- as.numeric(tr.y != 0)  ## The indicator (= 1 uncensored)
x <- as.matrix(tr.x[-cen])  ## New data without censored
y <- as.matrix(tr.y[-cen])

## Estimating G:
library(survival)
sfit <- survfit(Surv(tr.y, 1 - tr.d) ~ 1)
G <- summary(sfit, time = tr.x)$surv[-cen]
w <- 1 / G

K <- radial.kernel(x, x)
H = (2*w - 1) * K
l = lambda * (tau - 1) * w
u = lambda * tau * w

c = - y
A = rep(1, 190)
b = 0
r = 0

r.ipop <- ipop(c, H, A, b, l, u, r)


########## These codes are not used ################
lambda <- res$lambda[ind]
Dmat <- cbind(matrix(0, 3*n+1, 2*n+1), rbind(matrix(0, 2*n+1, n), radial.kernel(tr.x)/lambda))
dvec <- c(rep(tau, n), rep(1-tau, n), rep(0, n+1))
A1 <- cbind(matrix(0, n, n), diag(1, n), matrix(-1, n, 1), -radial.kernel(tr.x)/lambda)
A2 <- cbind(diag(1, n), matrix(0, n, n), matrix(1, n, 1), radial.kernel(tr.x)/lambda)
A3 <- cbind(diag(1, n), matrix(0, n, 2*n+1))
A4 <- cbind(matrix(0, n, n), diag(1, n), matrix(0, n, n+1))
#Amat <- t(rbind(A1, A2, A3, A4))
Amat <- rbind(A1, A2, A3, A4)
bvec <- -c(tr.y, tr.y, rep(0, n), rep(0, n))
l <- c(10^(-8), rep(tau - 1, n), rep(10^(-8), 2*n))
u <- c(10^8, rep(tau, n), rep(10^8, 2*n))
r <- rep(10^8, 4*n)

r.ipop <- ipop(c = dvec, H = Dmat, A = Amat, b = bvec, l = l, u = u, r = r, maxiter = 1000, margin = 10^(-2))
