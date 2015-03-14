# MCMC.3.1

n <- 100
y <- 57
step.length <- 0.001
theta <- seq(0, 1, step.length)
len <- length(theta)
results <- vector()
for.posterior <- vector()

# Flat prior on "steps+1" fixed choices for theta.
steps = 1/step.length
prior <- 1/(steps + 1)

for (i in 1:len) {
  likelihood <- choose(n, y)*theta[i]^y*(1-theta[i])^(n-y)
  results[i] <- likelihood
  for.posterior[i] <- likelihood*prior
}

plot(results)
sum(for.posterior)

plot(dbeta(seq(0, 1, 0.001), 58, 44))

