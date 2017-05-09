require(MASS)
require(geoR)

grid_w = 5
grid_h = 4

# Lab 3 - Assignment 1

data = read.table("data/rainfall.txt", header=FALSE)[,1]

n = length(data)

# (a)

pdf("plots/3_1_1_precipitation.pdf")

par(mfrow=c(2,1))

# Plot the precipitation
plot(data, 
     type="h", 
     lwd=2,
     ylab="Precipitation",
     xlab="Time",
     xaxt="n",
     col="black",
     main="Precipitation 1948-1983")

axis(1, 
     at=seq(0, n, n/(1983-1948)), 
     labels=seq(1948, 1983))

# Plot the precipitation density
prec_density = density(data)
plot(prec_density, 
     type="l", 
     lwd=2,
     xlab="Precipitation",
     ylab="Density",
     main="The Daily Precipitation")

dev.off()

data_mean = mean(data)

# Prior parameters for sigma (variance)
v0 = 1
sigma_sq0 = 1 / v0

n_draws = 4000

# Initial value for sigma
sigma_sq = rinvchisq(n=1, v0, sigma_sq0)


gibbs_draws = matrix(0,n_draws,2)
for(i in 1:n_draws){
  mu = rnorm(n=1, mean=data_mean, sd=sqrt(sigma_sq/n))
  sigma_sq = rinvchisq(n=1, v0 + n, (v0*sigma_sq0 + sum((data - mu)^2))/(n + v0)) 
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

pdf("plots/3_1_1_gibbs_sampl_conv.pdf")

par(mfrow=c(2,1))

# Plot the auto correlation (convergence) between draws of mu
min_mean = min(mean_means)
max_mean = max(mean_means)
plot(mean_means, 
     type="l", 
     ylim=c(min_mean, max_mean), 
     cex=.1,
     lwd=2,
     main=expression(paste("Convergence of Gibbs Sampling ", "(", mu, ")", sep=" ")),
     xlab="Batches of sequential draws",
     ylab=expression(paste("Mean of seq. draws of ", mu, sep=" ")))



# Plot the auto correlation (convergence) between draws of sigma
min_var = min(mean_vars)
max_var = max(mean_vars)
plot(mean_vars, 
     type="l", 
     ylim=c(min_var, max_var), 
     cex=.1,
     lwd=2,
     main=expression(paste("Convergence of Gibbs Sampling ", "(", sigma^2, ")", sep=" ")),
     xlab="Batches of sequential draws",
     ylab=expression(paste("Mean of seq. draws of ", sigma^2, sep=" ")))

dev.off()

# (c)

# Kernel density estimate
pdf("plots/3_1_3_dens_comp.pdf")

par(mfrow=c(3,1))

kernel_density = density(data)

plot(kernel_density$x, 
     kernel_density$y, 
     type="l",
     cex=.1,
     lwd=2,
     ylab="Density",
     xlab="Precipitation",
     main="Rainfall: Kernel density estimate")

# Normal density from (a)

mean = mean(mean_draws)
std_dev = sqrt(mean(var_draws)/n)

x_grid = seq(mean - 2, mean + 2, 0.0001)

normal_density = dnorm(x_grid, mean=mean, sd=std_dev)

plot(x_grid, 
     normal_density, 
     type="l",
     cex=.1,
     lwd=2,
     ylab="Density",
     xlab="Precipitation",
     main="Rainfall: Normal density")

# Mixture of normals from (b)
# (run './template/gaussian_mixture.R' first)

hist(x, breaks = 20, cex=.1, border="lightgray", freq = FALSE, xlim = c(xGridMin,xGridMax), xlab="Precipitation", ylab="Density", main = "Rainfall: Mixture of Normals")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "black")
lines(xGrid, dnorm(xGrid, mean = mean(x), sd = apply(x,2,sd)), type = "l", lwd = 2, col = "gray")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("lightgray","black","gray"), lwd = 2)


dev.off()

