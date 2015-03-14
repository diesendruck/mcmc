# MCMC.3.7.R


n <- 15
y <- 2
theta <- seq(0, 1, 0.001)
r <- choose(n, y)*theta^y*(1-theta)^(n-y)
plot(r)
abline(v=1000*2/15)

y <- 1:10
r <- choose(278, y)*factorial(2+y)*factorial(291-y)
plot(r)