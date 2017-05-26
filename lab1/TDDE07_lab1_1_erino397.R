grid_w = 6
grid_h = 5

# Experiment Setup
a0 = 2
b0 = 2
s = 14
f = 6

# Posterior alpha and beta
a_n = a0 + s
b_n = b0 + f

# Theoretical mean and standard deviations
mean = a_n / (a_n + b_n)
std_dev = sqrt(a_n*b_n / (((a_n + b_n)**2)*(a_n + b_n + 1)))

# 1.1

# Simulate from posterior beta(a0+s,b0+f)
draw_post_mean_st <- function(n_draw){
  post_draws = rbeta(n_draw, a_n, b_n)
  m = mean(post_draws)
  s = sd(post_draws)
  
  return(c(m, s))
}

# Calculate mean and variation for 2 to 10000 draws from posterior
data = sapply(2:10000, draw_post_mean_st)

# Save plot for the convergence of mean and variation

pdf("plots/1_1_mean.pdf", width=grid_w, height=grid_h)
  mean_values = data[1,]
  plot(mean_values, 
       type="l", 
       col="gray", 
       xlab="Iterations", 
       ylab="Mean", 
       main = '1.1 Mean convergence')
  abline(h=mean)
dev.off()

pdf("plots/1_1_std_dev.pdf", width=grid_w, height=grid_h)
  standard_deviations = data[2,]
  plot(standard_deviations, 
       type="l", 
       col="gray", 
       xlab="Iterations", 
       ylab="Standard Deviation", 
       main = '1.1 Standard deviation convergence ')
  abline(h=std_dev)
dev.off()


# 1.2

n_draws = 10000
post_draws = rbeta(n_draws, a_n, b_n)
p_draws = 100 * sum(post_draws <= 0.4) / length(post_draws)
p_true = 100 * pbeta(.4, a_n, b_n)

# 1.3

n_draws = 10000
post_draws = rbeta(n_draws, a_n, b_n)
lodds = log(post_draws/(1-post_draws))

pdf("plots/1_3_log_odds.pdf", width=grid_w, height=grid_h)
hist(lodds, 100, 
     prob=TRUE, 
     col="gray", 
     xlab="", 
     ylab="", 
     main = '1.3 Log Odds Simulation')
lines(density(lodds), lwd=2)
dev.off()
