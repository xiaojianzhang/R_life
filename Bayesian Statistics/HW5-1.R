#Chapter 5 Exercise 5
#setup
library(geoR)
ydata <- c(10, 10, 12, 11, 9)
n = 5
y_bar = 10.4
s_square = 1.3

#(c) How do the incorrect and correct posterior
# distributions differ?

#(1)Consider incorrect posterior
#draw sigma_square from inverse chi square(n-1, s^2)
sample_sigma_square <- rinvchisq(1000, n-1, s_square)

#draw miu from normal(y_bar, sigma_square/n)
sample_mu <- rnorm(1000, mean=y_bar, sd=sqrt(sample_sigma_square/n))

#sample mean
mean_mu <- mean(sample_mu)
mean_sigma_square <- mean(sample_sigma_square)
print(mean_mu)
print(mean_sigma_square)
#sample variance
var_mu <- var(sample_mu)
var_sigma_square <- var(sample_sigma_square)
print(var_mu)
print(var_sigma_square)

#contour plot
mu <- seq(0,20,0.05)
log_sigma <- seq(-2,4,0.02)
sigma <- exp(log_sigma)
log_post_mu_log_sigma <- function(mu, sigma, y){
  z <- 0
  for(i in 1:length(y)){
    z <- z + log(dnorm(y[i], mean=mu, sd=sigma))
  }
  return(z)
}
log_post <- outer(mu, sigma, log_post_mu_log_sigma, ydata)
post <- exp(log_post - max(log_post))
contours <- c(.0001,.001,.01,seq(.05,.95,.05))
contour (main="contour plot of incorrect posterior dist",
         mu, log_sigma, post, levels=contours,
         xlab="mu", ylab="log sigma", cex=2)

#(2)Consider correct posterior
log_post_mu_log_sigma <- function(mu, sigma, y){
  z <- 0
  for(i in 1:length(y)){
    z <- z + log(pnorm(y[i] + 0.5, mean=mu, sd=sigma) -
                   pnorm(y[i] - 0.5, mean=mu, sd=sigma))
  }
  return(z)
}
log_post <- outer(mu, sigma, log_post_mu_log_sigma, ydata)
post <- exp(log_post - max(log_post))
contours <- c(.0001,.001,.01,seq(.05,.95,.05))
contour (main="contour plot of correct posterior dist",
         mu, log_sigma, post, levels=contours,
         xlab="mu", ylab="log sigma", cex=2)

normalized_post <- post / sum(post)
post_mu <- rowSums(post)
mu_index <- sample (1:length(mu), 500, replace=T,
                   prob=post_mu)
mu_sample <- mu[mu_index]
for(i in 1:length(post_mu)){
  normalized_post[i, ] <- normalized_post[i, ] / post_mu[i] 
}
sigma_square_sample <- rep(NA, 500)
for(i in 1:500){
  sigma_square_sample[i] <- exp(sample(log_sigma, 1, 
                    prob=normalized_post[mu_index[i], ]))^2
}

#sample mean
mean_mu <- mean(mu_sample)
mean_sigma_square <- mean(sigma_square_sample)
print(mean_mu)
print(mean_sigma_square)

#sample variance
var_mu <- var(mu_sample)
var_sigma_square <- var(sigma_square_sample)
print(var_mu)
print(var_sigma_square)

#(d) draw simulatons from the posterior dist of z.
# compute the posterior mean of (z_1 - z_2)^2
z <- matrix(0, 500, 5)
for(i in 1:500){
  for(j in 1:5){
    rn = rnorm(1,mean=mu_sample[i],sd=sqrt(sigma_square_sample[i]))
    while(rn >= ydata[j] + 0.5 || rn <= ydata[j] - 0.5){
      rn = rnorm(1,mean=mu_sample[i],sd=sqrt(sigma_square_sample[i]))
    }
    z[i,j] <- rn
  }
}
#posterior mean of (z[1]-z[2])^2
print(mean((z[,1]-z[,2])^2))
