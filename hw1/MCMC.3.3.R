# MCMC Question 3.3 (p. 228)

################
#### PART A ####  Get the posterior distributions, means, variances, and
################  95% Confidence Intervals.

# Use this given data.
data.a = c(12, 9, 12, 14, 13, 13, 15, 8, 15, 6)
data.b = c(11, 11, 10, 9, 9, 8, 7, 10, 6, 8, 8, 9, 7)

# Find quick summaries.
sum.data.a <- sum(data.a)
sum.data.b <- sum(data.b)
count.data.a <- length(data.a)
count.data.b <- length(data.b)

# Set up parameters for priors.
param.a.prior.a <- 120
param.b.prior.a <- 10
param.a.prior.b <- 12
param.b.prior.b <- 1

# Draw samples of posteriors: gamma(param.a + sum.data, param.b + count.data)
# e.g.
# p(theta.a | y) = gamma(120 + 117, 10 + 10)
# p(theta.b | y) = gamma(12 + 113, 1 + 13)
posterior.a.params <- c(param.a.prior.a + sum.data.a,
                        param.b.prior.a + count.data.a)
posterior.b.params <- c(param.a.prior.b + sum.data.b,
                        param.b.prior.b + count.data.b)
posterior.sample.a <- rgamma(1000000, shape=posterior.a.params[1],
                                      rate=posterior.a.params[2])
posterior.sample.b <- rgamma(1000000, shape=posterior.b.params[1],
                                      rate=posterior.b.params[2])

# Calculate means and variances for each posterior.
mu.a = posterior.a.params[1] / posterior.a.params[2]
mu.b = posterior.b.params[1] / posterior.b.params[2]
var.a = posterior.a.params[1] / (posterior.a.params[2])^2
var.b = posterior.b.params[1] / (posterior.b.params[2])^2

# Find 95% confidence intervals.
quantiles.a <- quantile(posterior.sample.a, probs=c(0.025, 0.975))
quantiles.b <- quantile(posterior.sample.b, probs=c(0.025, 0.975))
plot(hist(posterior.sample.a))
plot(hist(posterior.sample.b))


##############
### PART B ###  Compute and plot expectation of theta.b under varios priors.
##############

# Set up parameter that will change in prior for theta.b.
changing.values <- 1:400

# Set up vector to store updated posterior means.
changing.means.b <- vector()

# Run new posterior samples using changing prior parameters.
for (n in 1:269) {
  param.a.prior.b <- 12*n
  param.b.prior.b <- n
  new.posterior.params <- c()
  new.posterior.params <- c(param.a.prior.b + sum.data.b,
                            param.b.prior.b + count.data.b)
  posterior.sample.b <- rgamma(10000, shape=new.posterior.params[1],
                                      rate=new.posterior.params[2])
  changing.means.b[n] <- mean(posterior.sample.b)
}

# Plot results. This shows the prior params of theta.b at which the
# expectation of theta.b to be equal to that of theta.a.
plot(changing.means.b, xlab="n, where gamma prior has a=12*n, b=n",
                       main=expression(paste("Variation of E(", theta[b],
                      ") given different gamma priors.")))
which(changing.means.b > 11.85)


