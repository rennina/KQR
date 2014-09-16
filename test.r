source("./kqr.r")
## 1-Dimensional
n <- 200
N1 <- 10000
N2 <- 10000
p <- 1

#### Validation Data
va.x <- matrix(runif(N1, 0, 1), nrow=N1, ncol=p)
va.f <- 2 * (exp( -30 * (va.x[,1]-0.25)^2 ) + sin(pi * va.x[,1]^2))
## Normal
e.n <- rnorm(N1, 0, 1)
va.y.n <- va.f + e.n
## Double Exponential
e.de <- rbinom(N1, 1, 0.5)
e.de[e.de==0] <- -1
e.de <- e.de * rexp(N1)
va.y.de <- va.f + e.de
## T3
e.t <- rt(N1, 3)
va.y.t <- va.f + e.t
## Mixture
e.m <- rbinom(N1, 1, 0.1)
e.m <- e.m * rnorm(N1, 0, 5) + (1-e.m) * rnorm(N1, 0, 1)
va.y.m <- va.f + e.m

#### Test Quantiles
te.x <- matrix(runif(N2*p, 0, 1), nrow=N2, ncol=p)
te.f <- 2 * (exp( -30 * (te.x[,1]-0.25)^2 ) + sin(pi * te.x[,1]^2))
tau <- 0.5
## Normal
q.n <- qnorm(tau, 0, 1)
te.f.n <- te.f + q.n
## Double Exponential
q.de <- log(2*tau)
te.f.de <- te.f + q.de
## T3
q.t <- qt(tau, 3)
te.f.t <- te.f + q.t
## Mixture
q.m <- uniroot(function(x) 0.1*pnorm(x, 0, 5) + 0.9*pnorm(x, 0, 1) - tau, low=-3, up=3, tol=1e-10)$root
te.f.m <- te.f + q.m


#### Degrees of Freedom
sic.n.df <- NULL
sic.de.df <- NULL
sic.t.df <- NULL
sic.m.df <- NULL
gacv.n.df <- NULL
gacv.de.df <- NULL
gacv.t.df <- NULL
gacv.m.df <- NULL
gold.n.df <- NULL
gold.de.df <- NULL
gold.t.df <- NULL
gold.m.df <- NULL

## Mean Check Deviation
sic.n.mcd <- NULL
sic.de.mcd <- NULL
sic.t.mcd <- NULL
sic.m.mcd <- NULL
gacv.n.mcd <- NULL
gacv.de.mcd <- NULL
gacv.t.mcd <- NULL
gacv.m.mcd <- NULL
gold.n.mcd <- NULL
gold.de.mcd <- NULL
gold.t.mcd <- NULL
gold.m.mcd <- NULL

for (j in 1:50) {
  cat("j:", j, "\n")
  tr.x <- matrix(runif(n*p, 0, 1), nrow=n, ncol=p)
  tr.k <- spline.kernel(tr.x, tr.x)
  va.k <- spline.kernel(tr.x, va.x)
  te.k <- spline.kernel(tr.x, te.x)

  tr.f <- 2 * (exp( -30 * (tr.x[,1]-0.25)^2 ) + sin(pi * tr.x[,1]^2))
  
### Standard Normal
  e.n <- rnorm(n, 0, 1)
  tr.y <- tr.f + e.n
  res <- kqrpath(tr.y, tau, K=tr.k)
  ## SIC
  ind <- order(res$Csic)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  sic.n.df[j] <- res$df[ind]
  sic.n.mcd[j] <- pf(fit, te.f.n, tau) / N2
  ## GACV
  ind <- order(res$Cgacv)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gacv.n.df[j] <- res$df[ind]
  gacv.n.mcd[j] <- pf(fit, te.f.n, tau) / N2
  ## Gold Standard
  va.g <- (res$alpha0.g0 + t(t(va.k) %*% res$alpha.g)) / res$lambda
  err.g <- apply(va.g, 1, pf, va.y.n, tau) / N1
  ind <- order(err.g)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gold.n.df[j] <- res$df[ind]
  gold.n.mcd[j] <- pf(fit, te.f.n, tau) / N2
  
### Double Exponential
  e.de <- rbinom(n, 1, 0.5)
  e.de[e.de==0] <- -1
  e.de <- e.de * rexp(n)
  tr.y <- tr.f + e.de
  res <- kqrpath(tr.y, tau, K=tr.k)
  ## SIC
  ind <- order(res$Csic)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  sic.de.df[j] <- res$df[ind]
  sic.de.mcd[j] <- pf(fit, te.f.de, tau) / N2
  ## GACV
  ind <- order(res$Cgacv)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gacv.de.df[j] <- res$df[ind]
  gacv.de.mcd[j] <- pf(fit, te.f.de, tau) / N2
  ## Gold Standard
  va.g <- (res$alpha0.g0 + t(t(va.k) %*% res$alpha.g)) / res$lambda
  err.g <- apply(va.g, 1, pf, va.y.de, tau) / N1
  ind <- order(err.g)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gold.de.df[j] <- res$df[ind]
  gold.de.mcd[j] <- pf(fit, te.f.de, tau) / N2
  
### T3
  e.t <- rt(n, 3)
  tr.y <- tr.f + e.t
  res <- kqrpath(tr.y, tau, K=tr.k)
  ## SIC
  ind <- order(res$Csic)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  sic.t.df[j] <- res$df[ind]
  sic.t.mcd[j] <- pf(fit, te.f.t, tau) / N2
  ## GACV
  ind <- order(res$Cgacv)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gacv.t.df[j] <- res$df[ind]
  gacv.t.mcd[j] <- pf(fit, te.f.t, tau) / N2
  ## Gold Standard
  va.g <- (res$alpha0.g0 + t(t(va.k) %*% res$alpha.g)) / res$lambda
  err.g <- apply(va.g, 1, pf, va.y.t, tau) / N1
  ind <- order(err.g)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gold.t.df[j] <- res$df[ind]
  gold.t.mcd[j] <- pf(fit, te.f.t, tau) / N2
  
### Mixture
  e.m <- rbinom(n, 1, 0.1)
  e.m <- e.m * rnorm(n, 0, 5) + (1-e.m) * rnorm(n, 0, 1)
  tr.y <- tr.f + e.m
  res <- kqrpath(tr.y, tau, K=tr.k)
  ## SIC
  ind <- order(res$Csic)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  sic.m.df[j] <- res$df[ind]
  sic.m.mcd[j] <- pf(fit, te.f.m, tau) / N2
  ## GACV
  ind <- order(res$Cgacv)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gacv.m.df[j] <- res$df[ind]
  gacv.m.mcd[j] <- pf(fit, te.f.m, tau) / N2
  ## Gold Standard
  va.g <- (res$alpha0.g0 + t(t(va.k) %*% res$alpha.g)) / res$lambda
  err.g <- apply(va.g, 1, pf, va.y.m, tau) / N1
  ind <- order(err.g)[1]
  fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind]
  gold.m.df[j] <- res$df[ind]
  gold.m.mcd[j] <- pf(fit, te.f.m, tau) / N2
  
  gc()
  gc()
  gc()
  gc()
  gc()
  gc()

}

