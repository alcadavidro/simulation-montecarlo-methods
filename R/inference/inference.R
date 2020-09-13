##################################################
####### MONTE CARLO METHODS FOR ESTIMATION #######
##################################################

## Let's begin with simply estimating a probability.
## Sometimes this is referred to as computing an
## expectation

## If you have a random variable X with a density function
## f(x) and we want to compute the expectation of f(x)

## Sam and Annie from 'sleepless in Seattle'
## Let A and S represent Sam's and Annie's arrival times
## at the Empire State Building, where we measure the 
## arrival time as the number of hours after noon.

## We assume that A and S are independent and uniformly
## distributed and that Annie arrives somewhere between 10:30
## and midnight and Sam arrives somewhere between 10:00 and
## 11:30 pm.

## Questions:
# 1) What is the probabilty that Annie arrives before Sam?
# 2) What is expected difference in arrival times?

## First question
# We simulate a large number of values from distribution of
# (A, S) say, 100, where A and S are independent:
obs <- 1000
sam <- runif(obs, 10, 11.5)
annie <- runif(obs, 10.5, 12)

(probabilty <- sum(annie < sam) / length(annie))
cat("Probability of Annie arrives before that Sam is", probabilty)

plot(sam, annie)
polygon(
  c(10.5, 11.5, 11.5, 10.5),
  c(10.5, 10.5, 11.5, 10.5),
  density = 10, angle = 135, col = 'red')

# Shaded region shows area A < S

# Standard error

(standard_error <- sqrt(probabilty * (1 - probabilty) / obs))

## Estimating an expectation: Question # 2
# What is the expected difference in the arrival times?

# Annie more likely to arrive later, so we model E(A - S)

difference <- annie - sam

# Monte carlo estimate is mean of these differences
mc.est <- mean(difference)
mc.est

# How is this estimation reliability? Let's calculate the standard error

se.est <- sd(difference) / sqrt(1000)
c(mc.est, se.est)

# So we estimate that Annie will arrive 0.485 hours
# after Sam arrives. SInce standard error is only 0.02,
# we can be 95% confident that the true difference is 
# within 0.04 hours of this estimate

#### GENERAL CASE WITH
#### STANDARD NORMAL DISTRIBUTION

# We generate a large number of random samples of size
# 2 from standard normal distribution, then compute the replicate
# pairs' mean differences, and then the mean of those differences

m <- 1000
g <- numeric(m)

for (i in 1:m){
  x <- rnorm(2)
  g[i] <- abs(x[1] - x[2])
}

(est <- mean(g))

# The first time we got 1.107158
# the second time got 1.142562

# can prove by integration that it equals 2/sqrt(pi) or 1.128379 and
# that the variance is 2-4/(pi) and the estandar error of estimate is

sqrt(sum((g - mean(g))^2)) / m

### Example of estimating the MSE of a trimmed mean

# A trimmed mean can be used to estimate the center of a continuous
# symmetric distribution that is not normal

# You get a trimmed mean by averaging all but the largest and smallest
# samples observations

# HEre the center is 0, we implement by writing a for loop (could also
# use the function replicate())

n <- 10   # vector is 20 long
m <- 1000 # go for 1000 replication
tmean <- numeric(m)

for (i in 1:m){
  x <- sort(rnorm(n))
  tmean[i] <- sum(x[2:(n-1)]) / (n -2)
}

(mse <- mean(tmean^2)) # determine the mse

sqrt(sum((tmean - mean(tmean))^2)) / m # standard error

