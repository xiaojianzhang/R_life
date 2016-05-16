ydata <- c(0.028, -0.741, -0.541, -0.246, 0.069, -0.584, -0.512,
           -0.079, -0.424, -0.335, -0.213, -0.039, -0.593, 0.282,
           -0.321, -0.135, 0.141, 0.322, 0.444, -0.218, -0.591, -0.608)

sd_data <- c(0.850, 0.483, 0.565, 0.138, 0.281, 0.676, 0.139, 0.204, 
             0.274, 0.117, 0.195, 0.229, 0.425, 0.205, 0.298, 0.261, 
             0.364, 0.553, 0.717, 0.260, 0.257, 0.272)

J <- 22

#(a)plot the posterior density of tau over an appropriate range
#   that includes essentially all of the posterior density
#using the formula (5.21) on page 117

tau <- seq(-5,5,0.01)
posterior_density_of_tau <- function(tau){
  z <- 1.0
  sigma_square <- sd_data^2
  V_miu <- 1.0/sum(1/(tau^2 + sigma_square))
  z<- z*sqrt(V_miu)
  miu_hat <- sum(ydata/(sigma_square+tau^2)) / sum(1.0/(tau^2 + sigma_square))
  tem <- (sigma_square+tau^2)^(-0.5) * exp(-(ydata-miu_hat)^2/(2*(sigma_square+tau^2)))
  z <- z * prod(tem)
  return(z)
}
posterior_density_tau <- sapply(tau, posterior_density_of_tau)
plot (tau, posterior_density_tau,xlim = c(0,1),
      main = "marginal posterior density of tau",
      type="l", xlab="tau", ylab="", xaxs="i",
      yaxs="i", yaxt="n", bty="n", cex=2)

#draw 1000 samples of tau from the posterior distribution of tau
sampleDist <- function(n, x, prob){ 
  return(sample(x, n, replace = T, prob))
}
sample_tau <- sampleDist(1000, tau, posterior_density_tau)
sample_tau <- sample_tau + runif(1000,-0.01/2,0.01/2)
sample_tau <- sort(sample_tau)
#(b)display how the posterior means and standard deiations
#of the theta_j's depend on tau
#simulating miu using formula (5.20)
sample_miu <- c()
for(i in 1:1000){
  sigma_square <- sd_data^2
  miu_hat <- sum(ydata/(sigma_square+sample_tau[i]^2)) / sum(1.0/(sample_tau[i]^2 + sigma_square))
  V_miu <- 1.0/sum(1.0/(sample_tau[i]^2 + sigma_square))
  sample_miu <- c(sample_miu, rnorm(1,mean=miu_hat,sd=sqrt(V_miu)))
}

#simulating theta_j using formula (5.17)
sample_theta <- matrix(0,1000,J)
for(i in 1:1000){
  for(j in 1:J){
    theta_hat <- ((ydata[j]/sd_data[j]^2) + 
                    (sample_miu[i]/sample_tau[i]^2))/((1.0/sd_data[j]^2)+(1.0/sample_tau[i]^2))
    V <- 1.0/((1.0/sd_data[j]^2)+(1.0/sample_tau[i]^2))
    sample_theta[i,j] <- rnorm(1,mean=theta_hat,sd=sqrt(V)) 
  }
}  

conditional_posterior_means <- matrix(0,1000,J)
conditional_posterior_variance <- matrix(0,1000,J)
for(i in 1:1000){
  miu_hat <- sum(ydata/(sigma_square+sample_tau[i]^2)) / sum(1.0/(sample_tau[i]^2 + sigma_square))
  V_miu <- 1.0/sum(1.0/(sample_tau[i]^2 + sigma_square))
  for(j in 1:J){
    conditional_posterior_means[i,j] <- ((ydata[j]/sd_data[j]^2) + (miu_hat/sample_tau[i]^2))/
    ((1.0/sd_data[j]^2) + (1.0/sample_tau[i]^2))
    
    conditional_posterior_variance[i,j] <- 1.0/((1.0/sd_data[j]^2) + (1.0/sample_tau[i]^2)) +
    ((1.0/sample_tau[i]^2)/((1.0/sd_data[j]^2) + (1.0/sample_tau[i]^2)))^2*V_miu
  }
}
matplot(sample_tau, conditional_posterior_means,
        type ="l",lty =1:J ,
        main ="conditional posterior means",
        xlab ="tau",ylab ="posterior mean", xlim = c(0,1))

matplot(sample_tau, sqrt(conditional_posterior_variance),
        type ="l",lty =1: J ,
        main ="conditional posterior standard deviation",
        xlab ="tau",ylab ="standard deviation", xlim = c(0,1))

#(c) scatterplot of the crude effect estimates vs
#   the posterior median effect estimates of the 22 studies
median_effect <- numeric(J)
for(i in 1:J){
  median_effect[i] <- median(sample_theta[,i])
}

plot(median_effect, ydata,
     main = "crude effect estimates vs posterior median effect",
     type= "p",
     xlab= "posterior median effectcrude",
     ylab= "crude effect estimates",
     pch = "o")

#(d)plot a histogram of the simulations of new treament effect
theta_new <- numeric(1000)
for(i in 1:1000){
  theta_new[i] <- rnorm(1,mean=sample_miu[i],sd=sqrt(sample_tau^2))
}

hist(theta_new, xlab="new treatment effect", yaxt="n",
      breaks=seq(-1.6, 1.1,.1), cex=2)