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

### Estimate MSE of the trimmed median

# actually, median is a trimmed mean
# So here we repeat for the median

n <- 20
m <- 1000
tmean <- numeric(m)

for (i in 1:m){
  x <- sort(rnorm(n))
  tmean[i] <- median(x)
}

(mse <- mean(tmean^2))
sqrt(sum((tmean - mean(tmean))^2)) / m

##### ESTIMATING A CONFIDENCE LEVEL

# Ofeten need evaluate cdf of sampling distribution
# of a statistic, when density function is unkwnon

# For example, often assume that sampled population
# is normally distributed. If population non-normal,
# true distribution of the estimator may be unkwown.

# This is a problem of integration, when you can
# generate the distribution g(x) but the true function
# of g(x) is unknown

### Example of condfience interval for variance

# consider confidence interval for variance, is
# sensitive to mild departures from normality. Can use
# MC methods to estimate true confidence level when
# normal theory confidence interval for variance
# is applied to non-normal data

# Here we calculate 95% upper confidence limit (UCL)
# for random sample size n=20 from N(0, var=4)
# distribution (sigma-squared is 4 in this case)

n <- 20
alpha <- 0.05

x <- rnorm(n, mean=0, sd=2)
# qchisq() is quantile function for chi-sq distr w/ df
(UCL <- (n-1) * var(x) / qchisq(alpha, df=n-1))

# Do it several times, upper confidence limit (UCL) all
# contain sigma-squared = 4

# If we do this large number of times, approximately
# 95% of the intervals should contain sigma-squared (=4)
# assuming sampled population is normal with a variance
# sigma-squared

### Example of MC estimate of confidence evel

# Simulation: repeat large number of times

n <- 20
alpha <- .05
UCL <- replicate(1000,
                 expr = {
                   x <- rnorm(n, mean=0, sd=2)
                   (n - 1)*var(x) / qchisq(alpha, df=n-1)
                 })

sum(UCL > 4)

# compute the mean to get the confidence level
mean(UCL > 4)

# Results is that 942 intervals satisfied (UCL > 4) so
# empirical confidence level is 94.2% very close to theorical
# level of 95%

# Could also accomplish call above like this

calcCI <- function(n, alpha){
  y <- rnorm(n, mean=0, sd=2)
  return((n-1)*var(y) / qchisq(alpha, df=n-1))
}
UCL <- replicate(1000, expr = calcCI(20, alpha = 0.05))
# count number of intervals that containes sigma² = 4

sum(UCL > 4)

### What if the population is non-normal ?
# Suppose have sampled pop of chi-sq(2) which has a variance of 4
# but is distincly non-normal

# We simply repeat the simulation, replacing N(0,4)
# samples with chi-sq(2) samples

n <- 20
alpha <- 0.05
UCL <- replicate(1000,
                 expr = {
                   x <- rchisq(n, df=2)
                   (n-1)*var(x) / qchisq(alpha, df=n-1)
                 })
sum(UCL>4)
mean(UCL > 4)

# In this experiment, only 794 if the intervals contained the
# population variance, which is far less than the 95%
# coverage

## THIS IS PARAMETRIC BOOTSTRAP

# The ordinary bootstrap is where samples are generated from observed
# samples where distribution is NOT specified, is sometimes called 
# 'non-parametric' and includes bootstrapping and jackknifing, which
# can be used to estimate the bias and the standard error of an estimate


## MORE on simulating sampling distribution with MC
# Suppose observe a random sample y1, ..., yn
# from exponential density function with unknown location
# parameter theta and you only observe sample median M

# simulatie median for exponential function
sim.median <- function(n) median(rexp(n))

M <- replicate(10000, sim.median(21))
hist(M, probability = T, main="")

# to confirm simulated draws are from
# sampling density of M, we define samp.med()
# to compute exact density

samp.med = function(m, n){
  con = factorial(n) / factorial((n - 1) / 2)^2
  con * pexp(m)^((n - 1) / 2) * (1 - pexp(m))^((n - 1) / 2) * dexp(m)
}

curve(samp.med(x, 21), add=TRUE, col="red")

## CONSTRUCTING A PERCENTILE CONF INTERVAL

quantile(M, c(0.05, 0.95))
#       5%      95% 
# 0.398545 1.102625 

# so 90% CI is (M-1.102625, M-0.398545)


##### COMPARING ESTIMATORS: THE TAXI PROBLEM
# A person is wandering the streets of a city
# and notices the following numbers of 5 taxis
# that pass by: 34, 100, 65, 81, 120

# Can she make an intelligent guess at the number of taxis in the city?
# Is a problem of statistical inference where population is collection
# of taxis driven in city and one wishes to know unknown number of
# taxis N

# Assume taxis are numbered 1 to N, each equally likely to be observed

# Consider two possible estimates: (1) the largest taxi number observed;
# and (2) twice the sample mean.

# The population is unknown

# We compare these two estimators using a simulation

# Simulate taxi number from a uniform distribution with
# known number of taxis N and compute the two estimates.
# Repeat many times and obtain two empirical sampling
# distributions

# Then can compare the two estimators by examining various
# properties of their respective sampling distribution

# taxi() function implements a single simulation.
# Two arguments: actual number of taxis N and
# sample size n.

taxi <- function(N, n){
  y <- sample(N, size=n, replace=T)
  estimate1 <- max(y)
  estimate2 <- mean(y) * 2
  c(estimate1=estimate1, estimate2=estimate2)
}

# Let's say actual number of taxis in city is 100.
# With n=5
taxi(100, 5)

# Which rule of thumb is more reliable?
# Simulation with sampling process of 1000 times

EST <- replicate(1000, taxi(100, 5))
str(EST)

# We got a matrix of 2 rows (N1, N2), 1000 columns

## ESTIMATING BIAS
# Want "unbiasedness" which means average value of
# estimator will be equal to parameter.

# We know number of taxis to be 100 and know std err

# Here we estimate their bias and standard errors
# Bias for first estimator
c(mean(EST["estimate1",]) - 100, sd(EST["estimate1",]) / sqrt(1000))
# Bias for second estimator
c(mean(EST["estimate2",]) - 100, sd(EST["estimate2",]) / sqrt(1000))

# It seems that N2 is less biased, but we can also compare them with
# respect to the mean distance from the parameter N (MAE)

### ESTIMATING MEAN DISTANCE FROM TARGET
# Compute mean absolute error and draw boxplot

absolute.error <- abs(EST - 100)
boxplot(t(absolute.error))

# using the errors
errors <- EST - 100
boxplot(t(errors))

# Use apply to find sample mean of the absolute
apply(t(absolute.error), 2, mean)

# std error
apply(t(absolute.error), 2, sd) / sqrt(1000)

### LATE TO CLASS AGAIN??

# Suppose the travel times for a particular student from home to school
# are normally distributed with mean 20 minutes and standard deviation
# 4 minutes. Each day during a five-day school week he leaves home 30 
# minutes before class. For each of the following problems, write a short
# Monte Carlo simulation function to compute the probability or expectation
# of interes

# (1) Find the Expected total traveling time of the student to school for
# five-day week. Find the simulation etimate and give the standard error
# for the simulation estimate.

travel.time <- function(n = 5){
  sum(rnorm(n=n, mean=20, sd = 4))
}

M <- replicate(1000, travel.time(n=5))
(mean(M))
# std error
(sd(M) / sqrt(1000))


# (2) Find the probability that the student is late for at least one class
# in the five-day week. Find the simulation estimate of the probability and the
# corresponding standard error

late.class <- function(){
  times <- rnorm(5, mean=20, sd=4)
  ifelse(all(times <= 30), "on.time", "late")
}

estimate.late.class <- replicate(1000, late.class())
prop.table(table(estimate.late.class))

# Expectation
p.hat <- mean(estimate.late.class == "late")
se.phat <- sqrt(p.hat * (1 - p.hat) / 1000)


# (3) On average, what will be the longest travel time to school during the 
# five-day week? Again find the simulation estimate with the standard error

max.time <- function(){
  max(rnorm(5, 20, 4))
}
longest.estimate <- replicate(1000, max.time())
mean(longest.estimate)
sd(longest.estimate) / sqrt(1000)
