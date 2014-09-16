source("./ipnw.R")

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
tr.f <- 40 * exp(8*(tr.x[, 1]-0.5)^2 + (tr.x[, 2]-0.5)^2) * (exp(8*(tr.x[, 1]-0.2)^2 + (tr.x[, 2]-0.7)^2) + exp(8*(tr.x[, 1]-0.7)^2 + (tr.x[, 2]-0.2)^2))^(-1)
e.n <- rnorm(n, 0, 1)   # normal error
tr.y <- tr.f + e.n

# 10000 validation data
va.x <- matrix(runif(N1*p, 0, 1), nrow = N1, ncol = p)
va.f <- 40 * exp(8*(va.x[, 1]-0.5)^2 + (va.x[, 2]-0.5)^2) * (exp(8*(va.x[, 1]-0.2)^2 + (va.x[, 2]-0.7)^2) + exp(8*(va.x[, 1]-0.7)^2 + (va.x[, 2]-0.2)^2))^(-1)
e.n <- rnorm(N1, 0, 1)   # normal error
va.y <- va.f + e.n

# 10000 test data
te.x <- matrix(runif(N2*p, 0, 1), nrow = N2, ncol = p)
te.f <- 40 * exp(8*(te.x[, 1]-0.5)^2 + (te.x[, 2]-0.5)^2) * (exp(8*(te.x[, 1]-0.2)^2 + (te.x[, 2]-0.7)^2) + exp(8*(te.x[, 1]-0.7)^2 + (te.x[, 2]-0.2)^2))^(-1)
## testing quantiles
#q.n <- qnorm(tau, 0, 1)   
#te.y <- te.f + q.n
e.n <- rnorm(N2, 0, 1)  # normal error
te.y <- te.f + e.n

### choose lambda using gold standard ###
lambdac <- seq(0.01, 3, by = 0.2) ## candidate lambda values
err.g <- NULL # prediction error

## define the check function
pf <- function(beta0, y0, tau) {
  tmp <- y0 - beta0
  sum(tau*tmp[tmp>0]) - sum((1-tau)*tmp[tmp<=0])
}

va.k <- radial.kernel(tr.x, va.x)
for (i in 1:length(lambdac)) {
	res <- ip.nw(tr.x, tr.y, lambdac[i], tau)
	va.g <- primal(res)[1] + t(va.k) %*% primal(res)[2:201] / lambdac[i]
	err.g[i] <- pf(va.g, va.y, tau) / N1
}	
ind <- order(err.g)[1]
lambda <- lambdac[ind]
##########################
