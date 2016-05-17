# generate prior function:
# Note: given info about witch's hat distribution, the distribution
# can be written as a piecewise function as follows:
witch_hat_prior <- function(theta){
  if(theta<0 | theta > 1){
    return (0)
  }else if(theta < 0.385 | theta > 0.585){
    return (0.5)
  }else if(theta < 0.485){
    return (50 * theta - 18.75)
  }else{
    return (29.75 - 50 * theta)
  }
}

n = 980
y = 437

posterior_density <- function(prior_den, n, y){
  
  theta <- seq(0, 1, 0.001)
  p_y_theta <- choose(n, y) * theta^y * (1-theta)^(n-y)  # likelihood
  p_theta <- sapply(theta, prior_den)   # prior
  
  prob_y <- function(y){      # prior predictive
    integrand <- function(x){
      return (choose(n, y) * x^y * (1-x)^(n-y) * witch_hat_prior(x))
    }
    p <- (integrate(integrand, 0, 1))
    return (p[[1]])
  }
  
  return ((p_y_theta * p_theta) / prob_y(y))   # posterior
  
}
theta <- seq(0, 1, 0.001)
post_den <- posterior_density(witch_hat_prior, n, y)
# (a) a plot of the posterior density
plot(theta, post_den, type = "l")

# (b) the posterior mean and variance
# evaluation by numerical integration(Trapezoidal rule)
integration <- function(integrand){
  sum <- 0.001 * integrand[1] / 2
  for(i in 2:(length(integrand)-1)){
    sum <- sum + 0.001 * integrand[i]
  }
  sum <- sum + 0.001 * tail(integrand,1) / 2
  return (sum)
}

integrand1 <- theta * post_den
mean <- integration(integrand1)

integrand2 <- (theta - mean)^2 * post_den 
variance <- integration(integrand2) 

# (c) the posterior median, and
# (d) the 95% central posterior interval

sum <- 0.001 * post_den[1] / 2
for(i in 2:1000){
  sum <- sum + 0.001 * post_den[i]
  if(sum > 0.5){
    median <- theta[i-1]
    break
  }
}

sum <- 0.001 * post_den[1] / 2
for(i in 2:1000){
  sum <- sum + 0.001 * post_den[i]
  if(sum > 0.025){
    left <- theta[i]
    break
  }
}

sum <- 0.001 * post_den[1] / 2
for(i in 2:1000){
  sum <- sum + 0.001 * post_den[i]
  if(sum > 1-0.025){
    right <- theta[i]
    break
  }
}

sprintf("mean = %f,  variance = %f", mean, variance)
sprintf("median = %f", median)
sprintf("95%% central posterior interval = [%f,%f]", left, right)
