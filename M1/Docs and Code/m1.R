
require(gtools)  # I found a Dir r.v. generator in this package
## if needed install it with
## install.packages("gtools")  

## DATA
setwd("~/Google Drive/2. SPRING 2015/MCMC - Prof Mueller/M1/Docs and Code")
x <- scan("sage.dta")                # raw data
y <- table(x)                        # counts
N <- length(y)       
names(y) <- 1:N
n <- length(x)

## HYPERPARS:
rho <- 0.1
as <- 0.9
bs <- 0.1
a1 <- .1
a0 <- 10

## initialize
## this function creates a list with
##   z=(z1,.. zN); pistar=pi*, r=(r~[1],..r~[M0])
##   q=(q~[1],..q~[M1])
## You can use it to initialize the state of the MC
init <- function() { # initialize parameters
  ## z
  z <- ifelse(y<10,0,1)
  ## pi*: empirical frquency
  A0 <- which(z==0);     A1 <- which(z==1) # Set of indices where z=0 / z=1.
  M0 <- length(A0);      M1 <- length(A1) # Num of indices where z=0 / z=1.
  Y0 <- sum(y[A0]);      Y1 <- sum(y[A1]) # Num of X's with z=0 / z=1.
  pistar <- sum(Y1)/n
  ## q and r: empirical fequencies
  q <- y[A1]/Y1  # this is the q~ of the text
  r <- y[A0]/Y0  # this is the r~ of the text
  return(th=list(z=z,pistar=pistar,q=q,r=r,A0=A0,A1=A1,
                 Y0=Y0,Y1=Y1,M0=M0,M1=M1))
}

## main function for MCMC
gibbs <- function(n.iter, verbose) {
  TH <- NULL # initialize - will save pi*,z here
             ##             for each iteration
  PI <- NULL # similar - will save (pi1,.., piN) here
  th <-  init()  
  pistar <- th$pistar        # initialize pistar = pi*
  z <- th$z                  # initialize z
  # TODO: Do the above need to be available in every iteration?
  
  for (i in 1:n.iter) { # loop over iterations
    A0 <- which(z==0);     A1 <- which(z==1) # Set of indices where z=0 / z=1.
    M0 <- length(A0);      M1 <- length(A1) # Num of indices where z=0 / z=1.
    Y0 <- sum(y[A0]);      Y1 <- sum(y[A1]) # Num of X's with z=0 / z=1.
    z <- sample.z(pistar, z, Y1, Y0, M1, M0)       # 1. z ~ p(z | pistar, y)
    q <- sample.q(pistar, z, A1, M1)       # 2. q ~ p(q | pistar,z,y)
    r <- sample.r(pistar, z, A0, M0)       # 3. r ~ p(r | pistar,z,y)
    pistar <- sample.pistar(Y1, Y0)    # 4. pi
    if (verbose > 0){
      if (i %% 10 == 0) {
        # print short summary, for every 10th iteration
        print(paste("Iteration: ", i))
        print(round(z, 2))
        print(round(q, 2))
        print(round(r, 2))
        print(round(pistar, 2))
      }  
    }
    ## save iteration
    TH <- rbind(TH, c(pistar,z))

    pi <- rep(0,N)
    pi[z==1] <- pistar*q[z==1]
    pi[z==0] <- (1-pistar)*r[z==0]
    PI <- rbind(PI, pi)
  }
  
  return(list(TH=TH, PI=PI))
}

ex <- function() { 
  ## RUN these lines to get the plots
  n.iter <- 1000; verbose=0
  gbs <- gibbs(n.iter, verbose)
  ## assume gbs returns a list with elements
  ## TH = (niter x p) matrix with each row being the
  ##      state (pi, z)
  ## PI = (niter x 1) vector with pi
  TH <- gbs$TH
  PI <- gbs$PI
  its <- 1:n.iter

  ## trajectory plot of Z
  plot(its, TH[,1],xlab="ITER",ylab="PI*",bty="l",type="l",
       main="Trajectories of PI*")

  ## boxplot
  boxplot(log(PI), main="Marginal Posterior Distributions of P(PIj|y)",
          xlab="PIj", ylab="Log(PI)")

  ## plotting posterior means vs. mle's
  pibar <- apply(PI,2,mean) # posterior means
  pihat <- as.numeric(y)/n
  
  plot(pihat, pibar, type="p", 
       pch=19, bty="l",xlab="MLE pihat", ylab="E(pi | y)",
       main="Posterior Means Versus MLEs")
  abline(0,1)

  ## same thing, zoom in to left lower corner
  plot(pihat, pibar, type="p", xlim=c(0,0.03), ylim=c(0,0.03),
       pch=19, bty="l",xlab="MLE pihat", ylab="E(pi | y)",
       main="Posterior Means Versus MLEs")
  abline(0,1)

}

#######################################################
## aux functions

sample.z <- function(pistar, z, Y1, Y0, M1, M0) {
  LogPostZ <- function(z) {
    A0 <- which(z==0);     A1 <- which(z==1) # Set of indices where z=0 / z=1.
    M0 <- length(A0);      M1 <- length(A1) # Num of indices where z=0 / z=1.
    Y0 <- sum(y[A0]);      Y1 <- sum(y[A1]) # Num of X's with z=0 / z=1.
    
    log(pistar^(Y1+as-1) * (1-pistar)^(Y0+bs-1) * 
                      rho^M1*(1-rho)^(M0)) +
    lgamma(M1*a1) - lgamma(a1)*M1 + lgamma(M0*a0) - lgamma(a0)*M0 +
    sum(lgamma((a1+y)^z)) - lgamma(sum(z*(a1+y))) +
    sum(lgamma((a0+y)^(1-z))) - lgamma(sum((1-z)*(a0+y)))
  }
  
  len.z <- length(z)
  for (i in 1:len.z) {
    z[i] <- 1
    lpz.1 <- LogPostZ(z)
    z[i] <- 0
    lpz.0 <- LogPostZ(z)
    # Posterior probability (actual, not proportional):
    # post(z=1)/(post(z=1)+post(z=0))
    thing <- lpz.1 - lpz.0
    pr <- exp(thing)/(exp(thing)+1)
    z[i] <- runif(1) < pr
  }
  
  return (z)
}

sample.q <- function(pistar, z, A1, M1) {
  dir.params <- rep(a1, M1)
  # Add the y1,...,yM1 to the Dirichlet params.
  for (i in 1:M1) {
    dir.params[i] <- dir.params[i]+y[A1][i]
  }
  new.q <- rdirichlet(1, dir.params)
  q <- new.q
  return (q)
}
sample.r <- function(pistar, z, A0, M0) {
  dir.params <- rep(a0, M0)
  # Add the yM1+1,...,yN to the Dirichlet params.
  for (i in 1:M0) {
    dir.params[i] <- dir.params[i]+y[A0][i]
  }
  new.r <- rdirichlet(1, dir.params)
  r <- new.r
  return (r)
}
sample.pistar <- function(Y1, Y0) {
  new.pistar <- rbeta(1, Y1+as, Y0+bs)
  pistar <- new.pistar
  return (pistar)
}


plot( cbind( y/sum(y), c(apply(TH,2,mean))[2:78] ) )
