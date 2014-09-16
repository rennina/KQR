
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
va.y.n <- va.f + e.n

# 10000 test data
te.x <- matrix(runif(N2*p, 0, 1), nrow = N2, ncol = p)
te.f <- 40 * exp(8*(te.x[, 1]-0.5)^2 + (te.x[, 2]-0.5)^2) * (exp(8*(te.x[, 1]-0.2)^2 + (te.x[, 2]-0.7)^2) + exp(8*(te.x[, 1]-0.7)^2 + (te.x[, 2]-0.2)^2))^(-1)
q.n <- qnorm(tau, 0, 1)   # normal error
te.f.n <- te.f + q.n

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


##### interior point algorithm: using ipop function #####
## non-weighted version (to compare with KQR paper):
ip.nw <- function(x, y, lambda, tau) {
	require(kernlab)
	H = radial.kernel(x, x)
	c = - y
	A = rep(1, n)
	b = 0
	r = 0
	l = matrix(lambda * (tau - 1), n, 1)
	u = matrix(lambda * tau, n, 1)
	ipop(c, H, A, b, l, u, r)
}

### choose lambda using gold standard ###
lambdac <- seq(0.1, 3, by = 0.2) ## candidate lambda values
err.g <- NULL # prediction error
## check function
pf <- function(beta0, y0, tau) {
  tmp <- y0 - beta0
  sum(tau*tmp[tmp>0]) - sum((1-tau)*tmp[tmp<=0])
}
va.k <- radial.kernel(tr.x, va.x)
for (i in 1:length(lambdac)) {
	res <- ip.nw(tr.x, tr.y, lambdac[i], tau)
	va.g <- t(t(va.k) %*% primal(res)) / lambdac[i]
	err.g[i] <- apply(va.g, 1, pf, va.y.n, tau) / N1
}	
ind <- order(err.g)[1]
lambda <- lambdac[ind]
##########################
