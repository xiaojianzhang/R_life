y <- c(16,9,10,13,19,20,18,17,35,55)
n <- c(16+58,9+90,10+48,13+57,19+103,20+57,18+86,17+112,35+273,55+64)
raw_proportions <- y/n

# (b) compute the marginal posterior density of 
# the hyperparameters and draw simulations from the 
# joint posterior distribution of the parameters and
# hyperparameter

#natural transformed scale:
# u <- log(alpha/beta)
# v <- log(alpha + beta)

u <- seq(-3,0,0.02)
v <- seq(1,5,0.02)

beta_matrix <- matrix(NA, length(u), length(v))
for(i in 1:length(u)){
  for(j in 1:length(v)){
    beta_matrix[i,j] <- exp(v[j])/(1+exp(u[i])) 
  }
} 
alpha_matrix <- matrix(NA, length(u), length(v))
for(i in 1:length(u)){
  for(j in 1:length(v)){
    alpha_matrix[i,j] <- beta_matrix[i,j] * exp(u[i]) 
  }
} 

# p(y|alpha, beta)
likelihood_density <- function(a, b, n, y){
  tem <- 1
  for(i in 1:10){
    tem <- tem * beta(a+y[i],b+n[i]-y[i])/beta(a,b)
  }
  return(tem)
}


prior_density <- function(alpha, beta){
  return(alpha * beta * (alpha + beta)^(-5/2))
}  

posterior_matrix <- matrix(NA, length(u), length(v))
for(i in 1:length(u)){
  for(j in 1:length(v)){
    posterior_matrix[i,j] <- prior_density(alpha_matrix[i,j], beta_matrix[i,j]) *
      likelihood_density(alpha_matrix[i,j], beta_matrix[i,j], n, y)
  }
}

max_post <- max(log(posterior_matrix))
log_post <- log(posterior_matrix) - max_post 
marginal_posterior <- exp(log_post)
contour(u,v,marginal_posterior, nlevels = 10, level = ,
        xlim = c(-3,0), ylim = c(1,5), main="Contour Plot",
        xlab="log(alpha/beta)",ylab="log(alpha+beta)")

#draw 1000 samples from marginal posterior distributions
#compute the marginal posterior distribution of u
mar_pos_u <- rowSums(marginal_posterior)

#draw 1000 samples of u from the marginal posterior distribution of u
sampleDist <- function(n, x, prob){ 
  return(sample(x, n, replace = T, prob))
}
sample_u <- sampleDist(1000, u, mar_pos_u)

#draw 1000 samples of v from the marginal posterior distribution of v given u
for(i in 1:length(u)){
  marginal_posterior[i, ] <- marginal_posterior[i, ]/mar_pos_u[i]
}
sample_v <- numeric()
for(i in 1:length(sample_u)){
  k <- (sample_u[i]-(-3))/0.02 + 1
  sample_v <- c(sample_v, sampleDist(1, v, marginal_posterior[k,]))
}

#For each of the sampled alpha and beta, add a uniform random jitter centered
#at zero with a width equal to the spacing of the sampling grid.
sample_u <- sample_u + runif(1000,-0.02/2,0.02/2)
sample_v <- sample_v + runif(1000,-0.02/2,0.02/2)

plot(sample_u, sample_v,
     main = "scatter plot",
     type= "p",
     xlab= "log(alpha/beta)",
     ylab= "log(alpha+beta)",
     xlim = c(-3,0),
     ylim = c(1,5),
     pch = ".")

# draw 1000 samples of theta 
sample_tehta <- matrix(NA, 1000, 10)
for(i in 1:1000){
  beta <- exp(sample_v[i])/(1+exp(sample_u[i])) 
  alpha <- beta * exp(sample_u[i])
  for(j in 1:10){
    sample_tehta[i,j] <- rbeta(1,alpha+y[j], beta+n[j]-y[j])
  }
}

# (c)compare the posterior distribution of the parameters
# theta_j to the raw proportions
plot(raw_proportions, raw_proportions,
     main = "compare the posterior distribution of the parameters
     theta_j to the raw proportions",
     type= "p",
     xlab= "raw proportions",
     ylab= "95% posterior interval for theta(i)",
     xlim = c(0,0.5),
     ylim = c(0,0.6),
     pch = 16)

for(i in 1:10){
  lines(c(raw_proportions[i], raw_proportions[i]),
       c(sort(sample_tehta[,i])[25],sort(sample_tehta[,i])[975]))
}

# (d) 95% posterior interval for the average underlying
# proportion of traffic that is bicycles

average_theta <- rowSums(sample_tehta)/10
interval <- sort(average_theta)[c(25,975)]
sprintf("95%% posterior interval for the number of fatal accidents in 1986 is [%f, %f]",
        interval[1], interval[2])

# (e) 95% posterior interval for the number of those vehicles
# that are bicycles
n_bicycles <- numeric()
for(i in 1:1000){
  n_bicycles <- c(n_bicycles, rbinom(1,100,average_theta[i]))
}
interval <- sort(n_bicycles)[c(25,975)]
sprintf("95%% posterior interval for the number of fatal accidents in 1986 is [%f, %f]",
        interval[1], interval[2])
