######################################################################
## These functions generate the FHT data introduced in:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
##        Regularization Paths for Generalized Linear Models via Coordinate Descent. 
##        Journal of Statistical Software, 33(1), 1-22. 
##        URL http://www.jstatsoft.org/v33/i01/.
##
## Copyright Jerome Friedman, Trevor Hastie, Robert Tibshirani
##
## The original comments and header are preserved.

genjerry=
function(x,snr){
# generate data according to Friedman's setup
n=nrow(x)
p=ncol(x)
b=((-1)^(1:p))*exp(-(2*(1:p)-1)/20)
f=x%*%b
e=rnorm(n)
k=sqrt(var(f)/(snr*var(e)))
y=f+k*e
return(y)
}

genx2=function(n,p,rho){
#    generate x's multivariate normal with equal corr rho
 if(abs(rho)<1){
  beta=sqrt(rho/(1-rho))
 x0=matrix(rnorm(n*p),ncol=p)
 z=rnorm(n)
 x=beta*matrix(z,nrow=n,ncol=p,byrow=F)+x0
 }
 if(abs(rho)==1){ x=matrix(z,nrow=n,ncol=p,byrow=F)}

return(x)
}
