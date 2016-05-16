# data
ydata <- c(16+58, 9+90, 10+48, 13+57, 19+103, 20+57, 18+86, 17+112, 35+273, 55+64)
# gammda parameter Gammda(phi, phi)
phi <- 1e-3

u <- seq(3,7,0.02)
v <- seq(-8,-1,0.02)

beta_matrix <- matrix(0, length(u), length(v))
for(i in 1:length(u)){
  for(j in 1:length(v)){
    beta_matrix[i,j] <- exp(v[j]) 
  }
} 
alpha_matrix <- matrix(0, length(u), length(v))
for(i in 1:length(u)){
  for(j in 1:length(v)){
    alpha_matrix[i,j] <- exp(v[j]+u[i]) 
  }
}

#
# (b) compute the joint posterior density of (alpha, beta)
# p(alpha, beta | y)
gamma_approx <- function(alpha,y){
  z <- sum(log(seq(alpha,alpha+y-1,1)))
  return(z)
}
posterior <- function(alpha,beta){
  z <- phi*log(alpha*beta)-phi*(alpha + beta) + 
    10*alpha*log(beta) - (10*alpha+sum(ydata))*log(beta+1) + 
    sum(mapply(gamma_approx,alpha,ydata))
  return(z)
}

posterior_matrix <- matrix(0, length(u), length(v))
for(i in 1:length(u)){
  for(j in 1:length(v)){
    posterior_matrix[i,j] <- posterior(alpha_matrix[i,j],beta_matrix[i,j])
  }
}

max_post <- max(posterior_matrix)
log_post <- posterior_matrix - max_post
posterior_matrix <- exp(log_post)
normalized_posterior <- posterior_matrix / sum(posterior_matrix)
contour(u, v, normalized_posterior, nlevels = 10, level = ,
        xlim = c(3,7), ylim = c(-8,-1), main="Contour Plot",
        xlab="log(alpha/beta)",ylab="log(beta)")

#draw 1000 samples from marginal posterior distributions
#compute the marginal posterior distribution of u
mar_pos_u <- rowSums(normalized_posterior)
#draw 1000 samples of u from the marginal posterior distribution of u
sampleDist <- function(n, x, prob){ 
  return(sample(x, n, replace = T, prob))
}
sample_u <- sampleDist(1000, u, mar_pos_u)

#draw 1000 samples of v from the marginal posterior distribution of v given u
for(i in 1:length(u)){
  normalized_posterior[i, ] <- normalized_posterior[i, ]/mar_pos_u[i]
}
sample_v <- numeric()
for(i in 1:length(sample_u)){
  k <- (sample_u[i]-3)/0.02 + 1
  sample_v <- c(sample_v, sampleDist(1, v, normalized_posterior[k,]))
}

sample_u <- sample_u + runif(1000,-0.02/2,0.02/2)
sample_v <- sample_v + runif(1000,-0.02/2,0.02/2)

plot(sample_u, sample_v,
     main = "scatter plot",
     type= "p",
     xlab= "log(alpha/beta)",
     ylab= "log(beta)",
     xlim = c(3,7),
     ylim = c(-8,-1),
     pch = ".")

#(e) draw 1000 samples of theta_j
sample_alpha <- exp(sample_u+sample_v)
sample_beta <- exp(sample_v)
theta <- matrix(0,1000,10)
for(i in 1:1000){
  for(j in 1:10){
    theta[i,j] <- rgamma(1,sample_alpha[i] + ydata[j], sample_beta[i]+1)
  }
}