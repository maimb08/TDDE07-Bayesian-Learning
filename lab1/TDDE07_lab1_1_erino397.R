
# Experiment Setup
a = 2
b = 2
n = 20
s = 14
p = s / n

# Posterior alpha and beta
a_n = a + n * p
b_n = b + n * (1 - p)

# Theoretical mean and standard deviations
mean = a_n / (a_n + b_n)
std_dev = sqrt(a_n*b_n / (((a_n + b_n)**2)*(a_n + b_n + 1)))

# 1.1

calc_mean_stddev <- function(n_draw){
  post_draws = rbeta(n_draw, a_n, b_n)
  m = mean(post_draws)
  s = sd(post_draws)
  
  return(c(m, s))
}

# Calculate mean and variation for 2 to 10000 draws from posterior
data = sapply(2:10000, calc_mean_stddev)

# Save plot for the convergence of mean and variation

pdf("1_mean")
  mean_values = data[1,]
  plot(mean_values, type="l", main = '1.1 Mean convergence')
  abline(h=mean)
dev.off()

pdf("1_std_dev")
  standard_deviations = data[2,]
  plot(standard_deviations, type="l", main = '1.1 Standard deviation convergence ')
  abline(h=std_dev)
dev.off()


# 1.2

p_less_point_four = pbeta(.4, a_n, b_n)

# 1.3

n_draws = 10000
post_draws = rbeta(n_draws, a_n, b_n)
lodds = log(post_draws/(1-post_draws))
hist(lodds, 100, prob=TRUE, col="purple")
lines(density(lodds), lwd=2)
