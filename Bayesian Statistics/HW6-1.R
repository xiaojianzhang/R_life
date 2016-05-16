#HW6-1:Chapter 11 Exercise 3.
library('geoR')
library('rstan')
ydata <- cbind(c(83,92,92,46,67),c(117,109,114,104,87),
                    c(101,93,92,86,67),c(105,119,116,102,116),
                    c(79,97,103,79,92),c(57,92,104,77,100))
J <- 6
n_j <- 5
n <- 30

theta_update <- function(){
  theta_hat <- (mu/tau^2 + colSums(ydata)/sigma^2)/(1/tau^2 + n_j/sigma^2)
  V_theta <- 1/(1/tau^2 + 1/sigma^2)
  rnorm(J, theta_hat, sqrt(V_theta))
}

mu_update <- function(){
  rnorm(1, mean(theta), tau/sqrt(J))
}

sigma_update <- function(){
  sigma_hat <- 0
  for(j in 1:J){
    sigma_hat <- sigma_hat + sum((ydata[,j]-theta[j])^2)
  }
  sqrt(rinvchisq(1, n, sigma_hat/n))
}

tau_update <- function(){
  sqrt(rinvchisq(1, J-1, sum((theta-mu)^2)/(J-1)))  
}

chains <- 5
iter <- 4000
sims <- array(NA, c(iter, chains, J+3))
dimnames(sims) <- list(NULL,NULL,
                       c(paste('theta[', 1:6, ']', sep=''), 'mu', 'tau', 'sigma'))
for(m in 1:chains){
  theta <- ydata[m,]
  mu <- mean(theta)
  for(t in 1:iter){
    tau <- tau_update()
    sigma <- sigma_update()
    theta <- theta_update()
    mu <- mu_update()
    sims[t,m,] <- c(theta, mu, tau, sigma)
  }
}
monitor(sims)

#(i) the posterior distribution
#of the mean of the quality measurements
#of the sixth machine
theta_sixth <- rowSums(sims[2001:4000, , 6])/5
hist(theta_sixth, breaks=50, 
    main='the posterior distribution
          of the mean of the quality measurements
          of the sixth machine',
    xlab = 'theta_6', freq = FALSE)

#(ii) the predictive
#distribution for another quality measurement
#of the sixth machine
pred <- rep(NA, 2000)
mu <- rowSums(sims[2001:4000, , 7])/5
tau <- rowSums(sims[2001:4000, , 8])/5
sigma <- rowSums(sims[2001:4000, , 9])/5
for(i in 1:2000){
  theta_hat <- (mu[i]/tau[i]^2 + sum(ydata[,6])/sigma[i]^2)/(1/tau[i]^2 + n_j/sigma[i]^2)
  V_theta <- 1/(1/tau[i]^2 + 1/sigma[i]^2)
  pred[i] <- rnorm(1, theta_hat, sqrt(V_theta))
}
hist(pred, breaks=50, 
     main='the predictive
          distribution for another quality measurement
          of the sixth machine',
     xlab = 'y_6', freq = FALSE)

#the posterior
#distribution of the mean of the 
#quality measurements of the seventh machine
theta_seventh <- rowSums(cbind(rowSums(sims[2001:4000, , 1])/5,
                       rowSums(sims[2001:4000, , 2])/5,
                       rowSums(sims[2001:4000, , 3])/5,
                       rowSums(sims[2001:4000, , 4])/5,
                       rowSums(sims[2001:4000, , 5])/5,
                       rowSums(sims[2001:4000, , 6])/5))/6
hist(theta_seventh, breaks=50, 
     main='the posterior
          distribution of the mean of the 
          quality measurements of the seventh machine',
     xlab = 'theta_7', freq = FALSE)