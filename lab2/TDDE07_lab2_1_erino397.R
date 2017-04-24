require(MASS)
require(geoR)
require(reshape2)

# Lab 2 - Assignment 1

df = read.table("data/TempLinkoping2016.txt", header=TRUE)

# (a)

# Linear model (test)
fit = lm(temp ~ time + I(time^2), data=df)
summary(fit)
#plot(fit)

n_params = 3

# Prior hyper parameters
mu0 = c(-10, 93, -85)
covar0 = diag(c(0.64, 2.98, 2.88))
v0 = 366 - n_params
sigma_sq0 = 1

# (b)

t = df[["time"]]
temp = df[["temp"]]

plot(t)
n_draws = 10
plot(t, temp, type="p", col="lightgray", ylim=c(-30, 30), xlim=c(0,1), xlab="Time", ylab="Temperature")
for (iter in 1:n_draws) {
  sigma_sq = rinvchisq(n=1, df=v0, scale=sigma_sq0)
  beta = mvrnorm(n=1, mu=mu0, Sigma=sigma_sq*ginv(covar0))
  error = rnorm(n=1, 0, 1)
  pred = beta[1] + beta[2]*t + beta[3]*I(t)^2 + error
  lines(t, pred)
}


