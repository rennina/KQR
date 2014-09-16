### This is to compare with the kqr function in the KQR paper ###
source("./kqr.r")
source("./ipnw.R")

############# Method in KQR paper ##############
tr.k <- radial.kernel(tr.x, tr.x)
va.k <- radial.kernel(tr.x, va.x)
te.k <- radial.kernel(tr.x, te.x)
res <- kqrpath(tr.y, tau, K=tr.k)

va.g <- (res$alpha0.g0 + t(t(va.k) %*% res$alpha.g)) / res$lambda
err.g <- apply(va.g, 1, pf, va.y, tau) / N1
ind <- order(err.g)[1]
fit <- (res$alpha0.g0[ind] + t(te.k) %*% res$alpha.g[,ind]) / res$lambda[ind] 
  
lambda <- res$lambda[ind] ### The gold standard lambda chosen by KQR

############# Interior point ##############
res.ip <- ip.nw(tr.x, tr.y, lambda, tau)
fit.ip <- primal(res.ip)[1] + t(te.k) %*% primal(res.ip)[2:201] / lambda

###### Weighted version #########

## Censored data:
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
