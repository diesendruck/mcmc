library(MASS)
library(mvtnorm)

## read data & hyperparameters

setwd("~/Google Drive/2. SPRING 2015/MCMC - Prof Mueller/M2")

## set the hyper parameters
Kmx <- 10               # max number of frequencies
b0 <- rep(0,2*Kmx+1)    # p(beta|K) = N(b0, A0^-1)
## (A0 precision matrix)
A0 <- diag(0.1,2*Kmx+1)
asig <- 1               # 1/sig2 ~ Ga(asig, bsig)
bsig <- 1
lambda <- 5

## read data
dta <- read.table("cepheid.dta",header=T,skip=15)
o <- order(dta$phase)   # sort the data by phase - for easier printing
dta <- dta[o,]
n <- nrow(dta)
x <- dta$phase
y <- dta$velocity

rTruncPois=function(lambda){
  r=rpois(1,lambda);
  while(r==0){
    r=rpois(1,lambda);
  }
  return(r);
}

make.x <- function(xx,k=1) { 
  ## make columns in the design matrix for
  ## k-th harmonic
  tm <- xx*(k*2*pi)
  return(cbind(sin(tm),cos(tm)))
}

## create the design matrix X for up to Kmx trig polynomials
##    use it by selecting the first (1+2K) columns for current K
X <- rep(1,n)
for(k in 1:Kmx) {
  X <- cbind(X,make.x(x,k))
}

## design matrix for plotting fits on a grid
n0 <- 100
X0 <- rep(1,n0)
x0 <- seq(from=0,to=1.1,length=n0)
for(k in 1:Kmx){
  X0 <- cbind(X0,make.x(x0,k))
}

## Model
## prior and likelihood evaluation

lpbeta <- function(b) { # evaluate log p(b)
  idx  <- 1:length(b)
  sd <- 1/sqrt(diag(A0)[idx])  # rely on A0 being diagonal!
  ## and recall A0 is precision matrix
  # lp <- sum(dnorm(b, m=b0[idx], sd=sd, log=TRUE))   # HIS VERSION
  lp <- sum(dnorm(b, m=b0[idx], sd=sqrt(10), log=TRUE))
  return(lp)
}

lpsig2 <- function(sig2) { ## evluate log p(1/sig2)=Ga(a0,b0)
  ## using parametrization of Ga such that
  ##       E[X]=a/b for X~Ga(a,b)
  ##       pdf p(x) = c*x^(a-1) exp(-x*b)
  lp <- log(dgamma(1/sig2, asig, bsig))
  return (lp)
}

loglik <- function(K, b, sig2) { # evaluate log likelihood 
  p <- length(b)
  n <- length(y)
  yhat <- X[,1:p] %*% b
  lp <- sum(dnorm(y, m=yhat, sd=sqrt(sig2), log=TRUE))
  return (lp)
}

lmarg <- function(K) { 
  # evaluates marginal, as in candidate's formula -- spell it out!!
  sig2 <- 0.5
  ## (i) find the posterior moments
  idx <- 1:(2*K+1)                  # select cols of design matrix 
  XK <- X[,idx]                     #     under K
  H <- t(XK) %*% XK
  bhat <- solve(H) %*% t(XK) %*% y  # least squares estimate
  Shatinv <- H/sig2                 #               inverse var
  S <- solve( A0[idx,idx]+Shatinv ) # posterior var (9.4) in Hoff
  bbar <- S %*% (1/sig2*t(XK)%*%y)  #          mean (9.5) in Hoff
  ## (ii) evaluate the marginal, using the formula from 1f
  lp <- 0.5*log(det(S)) +
    loglik(K,bbar,sig2) + lpbeta(bbar)
  return(lp)
}

ljointpost=function(K, b, sig2){
  #print("calculating posterior with param:");
  n = length(y);
  p <- length(b)
  yhat <- X[,1:p] %*% b
  
  # Posterior. Cancel out terms without K, yhat, b.
#   lpo  <- (-n/2)*log(sig2) + (-1/2*(1/sig2)*sum((y-yhat)^2)) - 1/sig2 
#   lpo <- lpo - 1/2*sum(b^2)/10 - lambda + K*log(lambda) - log(factorial(K))
  lpo  <- (-1/2)*(1/sig2)*sum((y-yhat)^2) - 1/2*sum(b^2)/10 + K*log(lambda) - 
            log(factorial(K))
  
  return(lpo);
}

## plots

plt.dta <- function(plt.spline=F) { 
  # plots data and adds smoothing spline (if plt.spine=T)
  plot(x,y, pch=19,bty="l",xlab="PHASE",ylab="VELOCITY",xlim=c(0,1.1),
       ylim=c(-25,50))
  if (plt.spline){ # add smoothing spline
    fit <- smooth.spline(x,y)
    lines(fit,col=3,type="l",lwd=2,lty=3)
  }
}

## MCMC

niter <- 100
mcmc <- function(niter=500) { # main loop for mcmc
  K <- 10     # initial values
  K <- rTruncPois(lambda)
  sig2 <- 1
  ## initialize lists to save imputed parameter values and mean function
  sig2list <- NULL
  Klist <- NULL
  flist <- NULL
  for(iter in 1:niter) {
    b <- sample.b(K,sig2)
    check.b(b)
    sig2 <- sample.sig2(K,b)
    th <- rj(K,b,sig2)
    K <- th$K
    b <- th$b
    
    ## update lists
    sig2list <- c(sig2list,sig2)
    Klist <- c(Klist,K)
    p <- length(b)
    idx <- 1:p
    f <- X0[,idx]%*%b
    flist <- rbind(flist,t(f))
    if (iter %% 50 == 0)# just for info as it runs..
      cat("\n it=",iter, "K=",K,"sig2=",format(sig2)) 
  }
  cat("\n")
  return(list(K=Klist,f=flist,sig2=sig2list))
}

# Gibbs transition probs

sample.b <- function(K,sig2) { # generate b ~ p(b | K, sig2, y)
  idx <- 1:(2*K+1)   # select columns (elements) for K harmonics
  Xk <- X[,idx]      # Subset of design matrix, with columns for setting of K.
  A0k <- A0[idx,idx] # Precision matrix of beta prior.
  b0k <- b0[idx]     # Mean vector of zeros from beta prior.
  
  # Full conditionals from Question 2
  V <- solve(t(Xk)%*%Xk/sig2+A0k)
  mm <- V%*%(t(Xk)%*%y/sig2)
  L <- t( chol(V))             # LL' = V
  b <- mm + L %*% rnorm(2*K+1) # b ~ N(m,V)
  
  return (b)
}

check.b <- function(b) { # for debugging -- quick plausibility check
  p <- length(b)
  yhat <- X[,1:p] %*% b
  z <- y-yhat
  s2 <- sum(z*z)/n
  if (s2 > 100){
    cat("\n *** large residual error = ",format(s2),"\n")
    info <- FALSE
  } else
    info <- TRUE
  return(info)
}

sample.sig2 <- function(K,b) { # generate 1/sig2 ~ p(1/sig2 | K,b,y)
  p <- length(b)
  idx <- 1:(2*K+1)   # select columns (elements) for K harmonics
  Xk <- X[,idx]      # Subset of design matrix, with columns for setting of K.
  
  a1 <- (n/2)+1
  b1 <- t(y-Xk%*%b)%*%(y-Xk%*%b)/2 + 1
  
  sig2inv <- rgamma(1,shape=a1,rate=b1)
  sig2 <- 1/sig2inv
  return (sig2)
}

## RJ

rj <- function(K, b, sig2) {
  ## prob of move up ("birth") or down ("death")
  for(j in 1:5) { # just to get more moves up & down
    q <- qbirth(K) # Returns q=1/2, or if K=1 then q=1.
    if (runif(1) < q) {
      th <- rj.birth(K,b,sig2)
    } else {
      th <- rj.death(K,b,sig2)
    }
    K <- th$K
    b <- th$b
  }
  return(th)
}

rj.birth <- function(K, b, sig2) { 
  ## birth move -- add one harmonic
  
  ## 1. generate auxiliaries (u1,u2) for the new regression coefficients
  ##    we use a normal linear regression of the residual on the
  ##    (K+1)-st harmonics
  u <- rnorm(2)
  ## 2. Transfer
  map.output <- Tmap(K, b, u, sig2) # Generate b1 ("b.tilde") with T(b,u) = b1
  b1 <- map.output$b1
  #dif <- sum((y-X[,1:length(b)]%*%b)^2) - sum((y-X[,1:length(b1)]%*%b1)^2)

  logJ <- map.output$logJ
  ## 2. acc prob
  r <- rho(K, b1, b, u, logJ, sig2)
  ## 3. accept (or not)
  if (runif(1) < r){ # accept with pr min(1,rho)
    b <- b1
    K <- K+1
    cat(" + ")
  }
  ## else reject (do nothing :-)
  return(th=list(K=K,b=b))
}

rj.death <- function(K1, b1, sig2) { 
  ## death move -- drop last harmonic
  ## it is convenient for notation to label now
  ## the current pars K1,b1
  ## proposed pars    K, b
  # 1. T inv mapping (b,u) = Tinv(b1)
  Tb <- Tinv(K1, b1, sig2)
  b <- Tb$b
  u <- Tb$u
  logJ <- Tb$logJ
  K <- K1-1
  ## 2. acc ratio (for opposite birth move)
  r <- rho(K, b1, b, u, logJ, sig2)
  ## 3. accept (or not)
  if (runif(1) < 1/r) {  # accept with pr min(1, 1/rho)
    b1 <- b
    K1 <- K1-1
    cat(" - ")
  }
  ## else reject (do nothing :-)
  return(th=list(K=K1,b=b1))
}    

qbirth <- function(K) {
  ifelse(K==1,1,
    ifelse(K==Kmx,0,0.5))
}

Tmap <- function(K, b, u, sig2) { 
  ## propose longer regr vector
  ## bnew = m + Lu or u = L^-1 (bnew-m)
  fit <- qu(K, b, sig2)
  bnew <- fit$m + fit$L %*% u
  ## compute Jacobian - save on log scale
  logJ <- sum(log(diag(fit$L)))
  b1 <- c(b,bnew)
  
  return(list(b1=b1,logJ=logJ))
}

qu <- function(K, b, sig2) { 
  ## find m,L for mapping T: (b,u) -> b1,  below
  ## bnew = m + L*u, 
  ## Use a regression of residuals on (K+1)-st harmonic
  ##     to determine m and L
  idx <- 1:(2*K+1)   # select columns (elements) for K harmonics
  Xk <- X[,idx]      # Subset of design matrix, with columns for setting of K.
  eps <- y-Xk%*%b
  regression <- lm(eps ~ sin((K+1)*2*pi*x) + cos((K+1)*2*pi*x))
  mk <- regression$coefficients[2:3]
  Vk <- vcov(regression)[2:3,2:3]
  Lk <- t(chol(Vk))
  return (list(m=mk, V=Vk, L=Lk))
}

Tinv <- function(K1,b1,sig2) {
  ## proposed (shorter) par vector
  ## bnew = m + Lu or u = L^-1 (bnew-m)
  K <- K1-1
  p <- 2*K+1
  b <- b1[1:p]  
  bnew <- b1[c(p+1,p+2)]
  
  ## back out auxiliary u, and logJ
  fit <- qu(K, b, sig2)
  u <- solve(fit$L)%*%(bnew-fit$m)
  logJ <- sum(log(diag(fit$L)))
  
  return (list(b=b, u=u, logJ=logJ))
}

rho <- function(K, b1, b, u, logJ, sig2) { 
  ##  acceptance ratio for birth move,
  ##  moving from b -> (b,u)
  # current parameter: b
  # proposal:          b1
  if (length(b) != 2*K+1) {# check
    cat("\n *** Error 1: ",
        "rho(.) should be called with rho(K,b1,b,u,J).\n")
  }
  if (length(b1) != length(b)+2) { # check
    cat("\n *** Error 2: ",
        "rho(.) should be called with rho(K,b1,b,u,J).\n")
  }
  K1 <- K+1
  
  lqu <- sum(pnorm(u, m=0, sd=1, log=TRUE)) # TODO: Shouldn't this be dnorm?

  ljointpostratio <- ljointpost(K1, b1, sig2) - ljointpost(K, b, sig2)
  
  # Priors and likelihood
  rho <- exp(ljointpostratio)
  
  # Transition probabilities
  if (K==1) {
    rho <- rho*2
  } else if (K==Kmx-1) {
    rho <- rho/2
  } else {
    birth.weight <- qbeta((Kmx-K+0.01)/Kmx, 1, 5)  # New rule: K close to Kmx, rho gets smaller.
    rho <- rho*birth.weight    #           K close to 1, rho gets bigger.
  }
  rho <- rho/exp(lqu)         # Auxiliary
  rho <- rho*exp(logJ)        # Jacobian

  return (rho)
}

# 1f
   
plt.marg <- function() {
   klist <- 1:Kmx
   lm <- rep(0,Kmx) # initialize
   for(k in klist) {
     lm[k] <- lmarg(k)
     lm <- lm-max(lm)
     lpk <- lm + dpois(klist,lambda,log=TRUE)
     pk <- exp(lpk)/(sum(exp(lpk)))
     barplot(pk,names=1:Kmx,xlab="K",ylab="p(K | y)",bty="l")
   }
}

# run it all
   
do.it <- function() {
   ## do the RJ MCMC
   M <- 10000
   mc <- mcmc(M)
   
   ## plots
   it <- 1:M
   plot(it,mc$K,bty="l",xlab="ITERATION",ylab="K",type="l",ylim=c(0,Kmx))
   
   idx <- it
   idx <- it[-(1:1000)] # drop first 1000 as initial transient
   hist(mc$K[idx], xlab="K", main="", bty="l",
     breaks=(4:10)+0.5, xlim=c(3,11), ylim=c(0,1), prob=T, ylab="p(K | y)")
   
   plot(it,mc$sig2,bty="l",xlab="ITERATION",ylab="SIG2",type="l")
   hist(mc$sig2[idx],xlab="SIG2",ylab="p(SIG2 | y)",main="", prob=T)
   
   plt.dta()
   mlist <- 1000+(1:18)*500
   for(m in mlist)
   lines(x0,mc$f[m,],type="l",lty=2)
   plt.marg()
   
   x0 <- seq(length=100,from=0,to=1)
   n0 <- length(x0)
   sig2 <- 1
   fbar <- mean(mc$f)
   y0 <- rnorm(n0,m=fbar,s=1)
   X0 <- rep(1,n0)
   for(k in 1:Kmx)
   X0 <- cbind(X0,make.x(x0,k))
   write.table(format(cbind(y0,x0),digits=2),file="yx.dta",
   row.names=F,quote=F,sep="\t")
   round(t(X0) %*% X0, digits=1)/100*2
}



do.it()



