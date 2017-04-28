require(MASS)
require(geoR)

grid_w = 5
grid_h = 4

# Lab 3 - Assignment 1

data = read.table("data/rainfall.txt", header=FALSE)[,1]

n = length(data)

# (a)

data_mean = mean(data)

# Prior parameters for sigma (variance)
v0 = 1
sigma0 = 1 / v0

n_draws = 4000

# Initial value for sigma
sigma = rinvchisq(n=1, v0, sigma0)

gibbs_draws = matrix(0,n_draws,2)
for(i in 1:n_draws){
  mu = rnorm(n=1, mean=data_mean, sd=(sigma/n))
  sigma = rinvchisq(n=1, n, (v0*sigma0 + sum((data - mu)^2))/(n + v0)) 
  gibbs_draws[i,] = c(mu, sigma)
}


mean_draws = gibbs_draws[,1]
var_draws = gibbs_draws[,2]


# Calculate mean of batches of 5 draws to visualize the
# auto correlation between sequential draws
mean_means = c()
mean_vars = c()
for (i in 1:n_draws){
  if(i%%5 == 0){
    mean_means = c(mean_means, mean(mean_draws[i-4:i]))
    mean_vars = c(mean_vars, mean(var_draws[i-4:i]))
  }
}

pdf("plots/3_1_gibbs_conv_mu.pdf", width=grid_w, height=grid_h)

# Plot the auto correlation (convergence) between draws of mu
min_mean = min(mean_means)
max_mean = max(mean_means)
plot(mean_means, 
     type="l", 
     ylim=c(min_mean, max_mean), 
     cex=.1,
     main=expression(paste("Convergence of Gibbs Sampling ", "(", mu, ")", sep=" ")),
     xlab="Batches of draws",
     ylab=expression(paste("Mean of batches of sequential draws of ", mu, sep=" ")))

dev.off()

pdf("plots/3_1_gibbs_conv_sigma.pdf", width=grid_w, height=grid_h)

# Plot the auto correlation (convergence) between draws of sigma
min_var = min(mean_vars)
max_var = max(mean_vars)
plot(mean_vars, 
     type="l", 
     ylim=c(min_var, max_var), 
     cex=.1,
     main=expression(paste("Convergence of Gibbs Sampling ", "(", sigma^2, ")", sep=" ")),
     xlab="Batches of draws",
     ylab=expression(paste("Mean of batches of sequential draws of ", sigma^2, sep=" ")))

dev.off()
