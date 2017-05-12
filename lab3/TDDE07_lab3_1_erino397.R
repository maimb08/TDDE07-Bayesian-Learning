require(MASS)
require(geoR)

grid_w = 5
grid_h = 4

# ----------------------
#  Lab 3 - Assignment 1
# ----------------------

data = read.table("data/rainfall.txt", header=FALSE)[,1]

n = length(data)

# -----
#  (a)
# -----

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

# -----
#  (b)
# -----

pdf("plots/3_1_3_dens_comp.pdf", height=grid_h)

kernel_density = density(data)

plot(kernel_density$x, 
     kernel_density$y, 
     type="l",
     cex=.1,
     lwd=2,
     ylab="Density",
     xlab="Precipitation",
     main="Rainfall: Kernel density estimate")

col1 = rgb(240, 240, 240, maxColorValue=255)
col1_b = "lightgray"
col2 = "black"
col3 = "black"

# Kernel density
polygon(kernel_density$x, 
        kernel_density$y,
        col=col1,
        border=col1_b,
        lwd=2)


# Mixture of normals from (b)
# (run './template/gaussian_mixture.R' first)
lines(xGrid, 
      mixDensMean, 
      type="l", 
      lwd=2,
      col=col2)


# Normal density
mean = mean(mean_draws)
std_dev = sqrt(mean(var_draws))
normal_density = dnorm(xGrid, mean=mean, sd=std_dev)

lines(xGrid, 
      normal_density, 
      col=col3,
      lwd=2,
      lty=2)


legend("topright", 
       box.lty = 0, 
       legend = c("Kernel Density","Mixture density","Normal density"), 
       col=c(col1_b, col2, col3), 
       lty=c(1,1,2),
       lwd=2)


dev.off()