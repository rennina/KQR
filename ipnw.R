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
	K <- radial.kernel(tr.x, tr.x)
	H <- matrix(0, 3*n+1, 3*n+1)
	H[2:(n+1), 2:(n+1)] <- K/lambda
	c <- c(rep(0, n+1), rep(tau, n), rep(1-tau, n))
	b <- tr.y
	r <- rep(0, n)
	k <- matrix(0, n, n)
	for (i in 1:n) {k[i,i] <- sum(K[i,])/lambda}
	A <- cbind(rep(1, n), k, diag(n), diag(-1, n, n))
	l <- c(rep(-10^5, (n+1)), rep(0, 2*n))
	u <- rep(10^5, 3*n+1)
	ipop(c, H, A, b, l, u, r)
}