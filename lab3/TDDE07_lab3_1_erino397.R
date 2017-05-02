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
sigma_sq0 = 1 / v0

n_draws = 4000

# Initial value for sigma
sigma_sq = rinvchisq(n=1, v0, sigma0)


gibbs_draws = matrix(0,n_draws,2)
for(i in 1:n_draws){
  mu = rnorm(n=1, mean=data_mean, sd=sqrt(sigma_sq/n))
  sigma_sq = rinvchisq(n=1, n, (v0*sigma_sq0 + sum((data - mu)^2))/(n + v0)) 
  gibbs_draws[i,] = c(mu, sigma_sq)
}


mean_draws = gibbs_draws[,1]
var_draws = gibbs_draws[,2]


# Calculate mean of batches of 2 draws to visualize the
# auto correlation between sequential draws
mean_means = c()
mean_vars = c()
for (i in 1:n_draws){
  if(i%%2 == 0){
    mean_means = c(mean_means, mean(mean_draws[i-1:i]))
    mean_vars = c(mean_vars, mean(var_draws[i-1:i]))
  }
}

# Plots displaying convergence of the Normal hyper 
# parameters during sampling

pdf("plots/3_1_1_conv_mu.pdf", width=grid_w, height=grid_h)

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


pdf("plots/3_1_1_conv_sigma.pdf", width=grid_w, height=grid_h)

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

# (c)

# Kernel density estimate
pdf("plots/3_1_3_kernel_dens_est.pdf", width=grid_w, height=grid_h)

kernel_density = density(data)

plot(kernel_density$x, 
     kernel_density$y, 
     type="l",
     cex=.1,
     ylab="Density",
     xlab="Precipitation",
     main="Rainfall: Kernel density estimate")

dev.off()


# Normal density from (a)
pdf("plots/3_1_3_normal_density.pdf", width=grid_w, height=grid_h)

mean_mean = mean(mean_draws)
mean_var = mean(var_draws)

x_grid = seq(mean_mean - 2, mean_mean + 2, 0.0001)

normal_density = dnorm(x_grid, mean=mean_mean, sd=sqrt(mean_var/n))

plot(x_grid, 
     normal_density, 
     type="l",
     cex=.1,
     ylab="Density",
     xlab="Precipitation",
     main="Rainfall: Normal density")

dev.off()


# - Mixture of normals in './template/gaussian_mixture.R'
