# MCMC - HW3.2

library(truncnorm)
library(MASS)

# RESOURCES
# http://en.wikipedia.org/wiki/Probit_model
# http://book.isito.kg/%D0%98%D0%BD%D1%84%D0%BE%D1%80%D0%BC%D0%B0%D1%82%D0%B8%D0%BA%D0%B0/%D0%9F%D1%80%D0%BE%D0%B3%D1%80%D0%B0%D0%BC%D0%BC%D0%B8%D1%80%D0%BE%D0%B2%D0%B0%D0%BD%D0%B8%D0%B5/R/Albert%20J.%20-%20Bayesian%20Computation%20with%20R,%202e%20-%202009.pdf
# http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_mcmc_sect053.htm
# http://fisher.osu.edu/~schroeder.9/AMIS900/ech7.pdf
# http://www.stat.cmu.edu/~brian/905-2009/all-papers/albert-chib-1993.pdf

Main <- function() {
  # Set up data.
  setwd("~/Google Drive/2. SPRING 2015/MCMC - Prof Mueller/hw3")
  contents <- read.csv("hw3.2.csv", header=T, sep=",")
  x <- as.matrix(contents[,c(1,2)])
  x <- cbind(1, x)
  y <- as.matrix(contents[,3])

  # Run Gibbs Sampler.
  iter <- 1000
  chains <- GibbsSampler(x, y, iter)
  b.chain <- chains[[1]]
  z.chain <- chains[[2]]
  # Thin each chain by extracting every 5th element.
  b.chain <- b.chain[, seq(1, dim(b.chain)[2], 5)]
  z.chain <- z.chain[, seq(1, dim(z.chain)[2], 5)]
  
  # Do traceplots and diagnostics for B.
  burnin <- 10
  len <- dim(b.chain)[2]
  num.betas <- dim(x)[2]
  par(mfrow=c(num.betas,2))
  estimate.b <- matrix(NA, nrow=num.betas, ncol=1)
  for (i in 1:num.betas) {
    hist(b.chain[i, ], 40, xlab="Beta Value", ylab="Frequency",
         main=paste("Histogram of Beta",i,"'s"))
    plot(ts(b.chain[i, burnin:len]), xlab="Iteration", ylab="Beta Value",
         main=paste("Mixing for Beta",i))
    estimate <- mean(b.chain[i, burnin:len])
    estimate.b[i, ] <- estimate
  }
  
  # Do traceplots and diagnostics for Z.
  n <- length(y)
  for (i in 1:n) {
    burnin <- 10
    #plot(ts(z.chain[i, burnin:len]), main=c("z", i), ylim=c(-3, 3))
  }
  estimate.z <- matrix(NA, nrow=n, ncol=1)
  for (j in 1:n) {
    estimate.z[j, ] <- mean(z.chain[j, burnin:len])
  }
  
  # Final estimates for B and Z.
  estimate.b
  estimate.z
}
  
GibbsSampler <- function(x, y, iter) {
  # Performs Gibbs sample to simulate joint posterior P(B,z|y).
  #
  # Args:
  #   x: Design matrix.
  #   y: Response vector.
  #   iter: Number of iterations.
  #
  # Returns:
  #   chains: Chains of estimated B and z random variables.
  num.betas <- dim(x)[2]
  n <- dim(x)[1]
  
  # Set starting values for beta (1's) and z (NA's) chains.
  b.chain <- matrix(1, nrow=num.betas, ncol=iter)
  z.chain <- matrix(0, nrow=n, ncol=iter)  
  
  # Initialize current variables with first elements of chains.
  b.i <- b.chain[,1]
  z.i <- z.chain[,1]
  
  for (i in 1:iter) {
    new.vars <- OneGS(b.i, z.i, n, x, y)
    b.i <- new.vars[[1]]
    z.i <- new.vars[[2]]
    b.chain[, i] <- b.i
    z.chain[, i] <- z.i
  }
  
  chains <- list(b.chain, z.chain)
  return (chains)
}

OneGS <- function(b.i, z.i, n, x, y) {
  # Runs one loop of Gibbs Sampler and returns new versions of B and Z.
  #
  # Args:
  #   b.i: Current beta vector.
  #   z.i: Current Z vector.
  #   n: Number of iterations.
  #   x: Design matrix.
  #   y: Response vector.
  #
  # Returns:
  #   output: List of updated beta and Z.
  new.z <- matrix(0, nrow=n, ncol=1)
  
  # Sample for Z's, given current betas.
  for (i in 1:n) {
    if (y[i]==1) {
      new.z[i,] <- rtruncnorm(1, a=0, b=Inf, mean=x[i,]%*%b.i, sd=1)
    } else if (y[i]==0) {
      new.z[i,] <- rtruncnorm(1, a=-Inf, b=0, mean=x[i,]%*%b.i, sd=1)
    }
  }
  z.i <- new.z
  
  # Sample for B, given new Z's.
  b.i <- mvrnorm(1, mu=solve(t(x)%*%x)%*%t(x)%*%z.i, Sigma=solve(t(x)%*%x))
  
  output <- list(b.i, z.i)
  return (output)
}

Main()
