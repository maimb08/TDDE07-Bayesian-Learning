require("geoR")

grid_w = 6
grid_h = 5


# 2.1

y = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
mean = 3.5
n = length(y)
n_draws = 10000

# Sample variance
tau_sq = sum((log(y) - mean)^2) / n
# X ~ chi(n)
X_draws = rchisq(n_draws, n)
# This is a draw from inv_chi(n, tausq)
sigma_sq = n * tau_sq / X_draws
# sigma_sq = rinvchisq(n_draws, n, tau_sq)
interval = seq(min(sigma_sq), max(sigma_sq), 0.001)
invchisq = dinvchisq(interval, n, tau_sq)

pdf("plots/2_1_chi_draws.pdf", width=grid_w, height=grid_h)
# Plot draws against density
hist(sigma_sq, 
     500, 
     border="gray", 
     prob=TRUE, 
     xlim=c(0,1.5), 
     xlab="Posterior of the variance", 
     ylab="", 
     main = '2.1 Posterior draws against posterior distribution')
lines(interval, invchisq, lwd=2)
dev.off()


# 2.2

z = sqrt(sigma_sq / 2)
G = 2*pnorm(z)-1

pdf("plots/2_2_gini_draws.pdf", width=grid_w, height=grid_h)
hist(G, 
     100, 
     xlab="Gini Coefficient", 
     ylab="", 
     main = '2.2 Posterior distribution of the Gini coefficient')
dev.off()

# 2.3

# Equal tail interval
equal_tail_interval = quantile(G, probs=c(0.025, 0.975))

# Highest Posterior Density
gd = density(G)
# Order x and y by y from high to low
ordered_x = gd$x[order(-gd$y)]
ordered_y = gd$y[order(-gd$y)]

# Iterate until 95% of prob. is captured
prob_mass = 0
total_mass = sum(gd$y)
for(i in 1:length(gd$y)){
  prob_mass = prob_mass + ordered_y[i]
  if(prob_mass / total_mass >= 0.95){
    break
  }
}

# Calculate the interval
a = min(ordered_x[1:i])
b = max(ordered_x[1:i])
highest_posterior_density = c(a, b)

pdf('plots/2_3_posterior_intervals.pdf', width=grid_w, height=grid_h)
plot(gd, col="red", lwd=2, main="2.3 Credibility intervals")
lines(equal_tail_interval, rep(0.12, 2), col="black", lwd=3)
lines(highest_posterior_density, rep(0.02, 2), col="gray", lwd=3)
legend("topright", 
       legend = c("ETI","HPD"),
       fill = c("black", "gray"),
       inset = 0.02)
dev.off()

