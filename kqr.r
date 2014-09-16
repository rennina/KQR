#######################################################################
### All rights are reserved by the authors.
### Authors: Youjuan Li, Yufeng Liu, Ji Zhu, University of Michigan (jizhu@umich.edu)
### Date:    01/01/2006
#######################################################################
###this defines the check function
pf <- function(beta0, y0, tau) {
  tmp <- y0 - beta0
  sum(tau*tmp[tmp>0]) - sum((1-tau)*tmp[tmp<=0])
}

###this calculate the initialization of data with y_i= a_th quantile
N.Initialization <- function(K, y, n, a, index, eps=1e-8) {
### K is n by n kernel matrix
### y is n by 1 response vector
### a is the quantile
### index is the initial point at the elbow
  Right <- seq(y)[y > y[index]]
  Left <- seq(y)[y < y[index]]
  nplus <- length(Right)
  nminus <- length(Left)
  Elbow <- index
  alpha.gamma <- rep(a, n)
  alpha.gamma[Elbow] <- nminus*(1-a) - nplus*a
  alpha.gamma[Left] <- -(1-a)
  f <- K %*% alpha.gamma
###find lambda at which new point enters the elbow, we can assume one point at a time
  lambdai <- double(n)
  lambdai[Elbow] <- -1
  lambdai[-Elbow] <- (f[Elbow] - f[-Elbow])/(y[Elbow] - y[-Elbow])
  lambda.entry <- max(lambdai)
  i <- match(lambda.entry, lambdai, 0)
  if (i == 0)
    stop("No match with lambda\n")
  else if (length(i) > 1)
    stop("More than two points arrive at the elbow\n")
  obs <- i
###don't need to consider the situation of point leaving the elbow
###because the point must stay on the elbow before the new point enters
  alpha0.gamma0 <- lambda.entry*y[Elbow] - f[Elbow] 
  Elbow <- union(Elbow, obs)
  list(Elbow=Elbow, Right=setdiff(Right,Elbow), Left=setdiff(Left,Elbow), lambda=lambda.entry, alpha0.gamma0=alpha0.gamma0, alpha.gamma=alpha.gamma)
}

###calculate the matrix of A after the event that point(s) leave the elbow 
DowndateKstar <- function(Kstar, index) {
  index <- index+1
  Kstar[-index, -index, drop=FALSE]
}

poly.kernel <- function(x, y=x, param.kernel=1) {
  if (param.kernel == 1)
    x %*% t(y)
  else (x %*% t(y) + 1)^param.kernel
}

radial.kernel <- function(x, y=x, param.kernel = 1) {
###Note param.kernel is now gamma 
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  exp(-a / (2 * param.kernel^2))
}

spline.kernel <- function(x, y=x) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  p <- ncol(x)
  n <- nrow(x)
  m <- nrow(y)
  kk <- matrix(1,n,m)  
  for (j in 1:p) {
    k1 <- (x[,j]-0.5) %o% (y[,j]-0.5)
    k2 <- 1/4 * ((x[,j]-0.5)^2 - 1/12) %o% ((y[,j]-0.5)^2 - 1/12)
    k3 <- abs(outer(x[,j], y[,j], "-"))
    k3 <- ((k3 - 0.5)^4 - 1/2 * (k3 - 0.5)^2 + 7/240) / 24
    kk <- kk * (1 + k1 + k2 - k3)
  }
  return (kk)
}

eta.kernel <- function(x, y=x,param.kernel=1) {
  ###Note there's no param.kernel 
  n <- nrow(x)
  m <- nrow(y)
  p <- ncol(x)
  normx <- drop((x^2) %*% rep(1, p))
  normy <- drop((y^2) %*% rep(1, p))
  a <- x %*% t(y)
  a <- (-2 * a + normx) + outer(rep(1, n), normy, "*")
  for (i in 1:n) {
  for (j in 1:m){
  if (a[i,j]<=1e-9) {a[i,j] <- 0}
  else a[i,j] <- {(1/(8*pi))*a[i,j]*0.5*log(a[i,j])}
  }
  }
  a
}

SolveKstar <- function(Kstar, eps=.0001) {
  onestar <- rep(1, ncol(Kstar))
  onestar[1] <- 0
  ##solve(Kstar+diag(length(onestar))*eps,onestar)
  ##solve(Kstar, onestar, tol=.Machine$double.eps)
  tmpqr <- qr(Kstar, tol=.Machine$double.eps)
  if (tmpqr$rank < dim(Kstar)[1]) {
    return( NULL )
  } else {
    return( qr.solve(tmpqr, onestar) )
  }
}

kqrpath <- function(y, a, K, eps=1e-10, Nmoves = 5000, lambda.min=1e-6) {
  require(boot)

	### Add jitter to y
  tmpy <- unique(y)
  if (length(tmpy) < length(y)) {
    dif <- min(diff(sort(tmpy)))
    y <- y + runif(length(y), -dif/10000, dif/10000)
  }
  
	### Initializations
  n <- length(y) 
  I <- rep(1,n) 
  Right <- Elbow <- NULL
  Left <- seq(n)
  Kscript <- K * outer(1/y, I)   

	### We start with a maximum of 3*n moves, but it can be increased
	### Initializations of counters
  alpha.g <- matrix(1, n, Nmoves)
  alpha0.g0 <- double(Nmoves)
  Cgacv <- double(Nmoves)
  Csic <- double(Nmoves)
  checkf <- double(Nmoves)
  df <- integer(Nmoves)
  fits <- matrix(1, n, Nmoves)
  Elbow.list <- as.list(seq(Nmoves))
  Size.Elbow <- integer(Nmoves)
  lambda <- double(Nmoves)
  Kstar <- matrix(0, 1, 1)

	###Initialization of path
	### Two cases 
  yr <- sort(y)
  beta0 <- yr[ floor(n*a)+1 ]
  index <- match(beta0, y)
  if ( (n*a) != floor(n*a) ) {
    init <- N.Initialization(K, y, n, a, index)
  }else {
    ## beta0 is the right end of the interval
    init <- I.Initialization(K, y, n, a, index)
  }

  Elbow <- init$Elbow
  Right <- init$Right
  Left <- init$Left
  Elbow.list[[1]] <- Elbow
  df[1] <- length(Elbow)
  lambda0 <- init$lambda
  alpha0.g0[1] <- init$alpha0.gamma0
  alpha.g[,1] <- init$alpha.gamma
  Kstar <- UpdateKstar(Kstar, Kscript[Elbow,Elbow], NULL, NULL, y[Elbow])
  lambda[1] <- lambda0
  fl <- (K %*% alpha.g[, 1] + alpha0.g0[1])/lambda0
  k <- 1

  seqN <- 1:length(tr.y)

  ## Compute the path
  while ((k < Nmoves) && (lambda[k] > lambda.min) ) {
  
	cat("k:", k, "\n")
	###Implement the model selection critera SIC and GACV, Validation:  
    	fits[,k] <- (alpha0.g0[k] + K %*% alpha.g[,k]) / lambda[k]
    	checkf[k] <- pf(fits[,k], y, a)
    	trHat <- length(Elbow.list[[k]])
    	Cgacv[k] <- checkf[k] / (n-trHat)
    	Csic[k] <- log(checkf[k]/n) + (log(n)/(2*n))*trHat

	### Consider the situation when the elbow becomes empty
    	if ((length(Elbow) == 0) && ((a*n) != floor(a*n)))
      	stop("Elbow cann't be empty for this data\n")
   
 	if ((length(Elbow) == 0) && ((a*n) == floor(a*n))) {
	
		### The elbow has become empty; need to resort to an initial condition
	      	
		cat('Problematic section\n')
		
		# Problem found in this section
		beta0 <- alpha0.g0[k] / lambda[k]
      	lamb <- lambda[k]
     	 	L <- Left
      	init <- Init.Initialization(K, y, n, a, beta0, L, lamb)
      	lambda0 <- init$lambda
      	alpha0.g0[k+1] <- init$alpha0.gamma0
      	Elbow <- init$Elbow
      	Left <- init$Left
      	Right <- init$Right
      	lambda[k+1] <- lambda0
      	alpha.g[,k+1] <- init$alpha.gamma
      	Kstar <- UpdateKstar(Kstar, Kscript[Elbow,Elbow], NULL, NULL, y[Elbow])
      	fl <- (K %*% alpha.g[, k+1] + alpha0.g0[k+1]) / lambda0
    
  	} else {
      	bstar <- SolveKstar(Kstar)
      	if (is.null(bstar)) {
        		cat("Singular Kstar\n")
        		break
  		}

      	b0 <- bstar[1]
      	b <- bstar[-1]
		
		### Now find the first event
		### Check for immobile margin
      	gl <- as.matrix(K[, Elbow]) %*% b + b0
      	dl <- fl - gl
      	immobile <- (sum(abs(dl))/n) < eps
		
		### now check for exits from Elbow
      	temp <-  - alpha.g[Elbow, k] + lambda[k] * b
      	lambda.right <- (a + temp) / b
      	lambda.right[abs(b) < eps] <- -1 #anything negative
      	lambda.left <- (-(1-a) + temp) / b
      	lambda.left[abs(b) < eps] <- -1
      	lambda01 <- c(lambda.right, lambda.left)
      	lambda.exit <- max(0, lambda01[lambda01 < lambda[k] - eps])


		### Check to see if we leave the margin when it is immobile
      	if (immobile & (lambda.exit < eps)) break

		### Now check for entries
      	if (!immobile) {
        		lambdai <- (lambda[k] * dl) / (y - gl)
        		lambdai[abs(y-gl) < eps] <- -Inf
        		lambdai[Elbow] <- -Inf  
        		lambda.entry <- max(-1, lambdai[lambdai < lambda[k] - eps])
      
		} else {
        		lambda.entry <- -1  #any negative will do
      
		}
      	lambda.max <- max(lambda.entry, lambda.exit)  

		### update lambda, alphas and fit
      
		lambda[k + 1] <- lambda.max
      	alpha.g[, k + 1] <- alpha.g[, k]
      	alpha.g[Elbow, k + 1] <- alpha.g[Elbow, k] - (lambda[k] - lambda.max) * b
      	alpha0.g0[k + 1] <- alpha0.g0[k] - (lambda[k] - lambda[k + 1]) * b0
      	fl <- (lambda[k] / lambda[k + 1]) * dl + gl

		### update active sets
      	if (lambda.entry > lambda.exit) {
	
			###point joins the elbow
        		i <- match(lambda.entry, lambdai, 0)
        		if (length(i) > 1)
          			stop("more than one point enters the elbow\n")
	
			###assumes for now there is only 1
        		if (match(i, Left, FALSE)) {
          			Left <- setdiff(Left, i)
        		} else {
          			Right <- setdiff(Right, i)
        		}
        		Kstar <- UpdateKstar(Kstar, Kscript[i, i], drop(Kscript[i, Elbow]), drop(Kscript[Elbow,i]), y[i])
        		Elbow <- c(Elbow, i)
      
		} else {
	
			###point(s) leaves the elbow; can be more than one
        		idrop <- Leaveright <- NULL
        		i <- Elbow[abs(lambda.right-lambda.exit) < eps]
        		if (length(i) > 0) {
          			Leaveright <- rep(TRUE, length(i))
          			idrop <- i
        		}
        		i <- Elbow[abs(lambda.left-lambda.exit) < eps]
        		if (length(i) > 0) {
          			Leaveright <- c(Leaveright, rep(FALSE, length(i)))
          			idrop <- c(idrop,i)
        		}
        
			cat("idrop:", length(idrop), "\n")
        		cat("leaveright:", Leaveright, "\n")
        
			for (j in seq(along=idrop)) {
          			if (Leaveright[j]) {
            		Right <- c(Right, idrop[j])
          			} else {
            		Left <- c(Left, idrop[j])
          			}

          			mi <- match(idrop[j], Elbow)
          			Kstar <- DowndateKstar(Kstar,mi)
          			Elbow <- Elbow[-mi]
        		} # for
      	} # else
    	} # else
    	k <- k + 1
    	Elbow.list[[k]] <- Elbow
    	df[k] <- length(Elbow)

  } # while

  obj <- list(alpha.g=alpha.g[,seq(k-1)], alpha0.g0=alpha0.g0[seq(k-1)], lambda=lambda[seq(k-1)], Elbow=Elbow.list[seq(k-1)], Cgacv=Cgacv[seq(k-1)], Csic=Csic[seq(k-1)], checkf=checkf[seq(k-1)], df=df[seq(k-1)], y=y, fits = fits[,1:(k-1)])
  class(obj)<-"kqrpath"
  obj
}

###this calculate the initialization of data with y_i!= a_th quantile
I.Initialization <- function(K, y, n, a, index, eps=1e-8, Nmoves=3*n){
### beta0 is the right end of the interval
  lambda.temp <- double(Nmoves)
  i.temp <- integer(Nmoves)
  beta0 <- y[index]
  Right <- seq(y)[y >= beta0]
  Left <- seq(y)[y < beta0]
  Elbow <- NULL
  lambdai <- double(n)
  alpha.gamma <- rep(a,n)
  alpha.gamma[Left]<- -(1-a)
  f <- K %*% alpha.gamma 
  
  seqN <- 1:length(y)

  # A1 %*% x <= b1;  A2 %*% x >= b2  and x >=0
  A1 <- cbind(rep(1,length(Right)), rep(-1, length(Right)), f[Right])
  A2 <- cbind(rep(1,length(Left)), rep(-1, length(Left)), f[Left])
  b1 <- y[Right]
  b2 <- y[Left]

  cost <- c(0,0,1)
  g <- simplex(cost, A1, b1, A2,b2, maxi=T)

  beta0 <- g$soln[1] - g$soln[2]
  lambda <- 1/g$soln[3] 
  Elbow <- seqN[abs(y - beta0 - g$soln[3]*f)<1e-10]
  Left <- setdiff(Left,Elbow)
  Right <- setdiff(Right, Elbow)

  alpha0.gamma0 <- lambda*beta0

  list(Elbow=Elbow, Right=Right, Left=Left, lambda=lambda, alpha.gamma=alpha.gamma, beta0=beta0, alpha0.gamma0=alpha0.gamma0)
}

###this is to resort the initialization of data with y_i!= a_th quantile(with elbow empty)
Init.Initialization<-function(K, y, n, a, beta0, L, lamb, eps=1e-8, Nmoves=3*n){
  lambda.temp <- double(Nmoves)
  i.temp <- integer(Nmoves)
  if (length(L) != a*n) stop("error")
  Left <- L
  Elbow <- NULL
  Right <- setdiff(seq(y), L)
  alpha.gamma <- rep(a,n)
  alpha.gamma[L] <- -(1-a)
  f <- K %*% alpha.gamma
  lambdai <- f / (y - beta0)

  seqN <- 1:length(y)

  # A1 %*% x <= b1;  A2 %*% x >= b2  and x >=0
  A1 <- cbind(rep(1,length(Right)), rep(-1, length(Right)), f[Right])
  A2 <- cbind(rep(1,length(Left)), rep(-1, length(Left)), f[Left])
  b1 <- y[Right]
  b2 <- y[Left]

  cost <- c(0,0,1)
  g <- simplex(cost, A1, b1, A2,b2, maxi=T)

  beta0 <- g$soln[1] - g$soln[2]
  lambda <- 1/g$soln[3] 
  Elbow <- seqN[abs(y - beta0 - g$soln[3]*f)<1e-10]
  Left <- setdiff(Left,Elbow)
  Right <- setdiff(Right, Elbow)

  alpha0.gamma0 <- lambda*beta0

  list(Elbow=Elbow, Right=Right, Left=Left, lambda=lambda, alpha.gamma=alpha.gamma, alpha0.gamma0=alpha0.gamma0)
}

###calculate the matrix of A after the event that point(s) enter the elbow    
UpdateKstar <- function(Kstar, Kell, Krow, Kcol, y) {
  yreverse <- 1/y
  I <- rep(1, length(y))
  if (length(y) == 1) {
    advec <- c(yreverse, Krow)
    aI <- c(I, Kcol)
    Kstar <- rbind(Kstar, advec)
    cbind(Kstar, c(aI, Kell))
  }
  else {
    adrect <- cbind(yreverse, Krow)
    adI <- rbind(I, Kcol)
    Kstar <- rbind(Kstar, adrect)
    cbind(Kstar, rbind(adI, Kell))
  }
}
