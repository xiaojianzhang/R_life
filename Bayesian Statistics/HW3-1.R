# (b)sketch contours of informative
# prior distribution

alpha <- seq(23.9,32.9,0.05)
beta <- seq(-1.02,-0.82,0.001)
#informative prior density function
info_prior_density <- function(x1,x2,mu1,sigma1,mu2,sigma2,corr){
  z <- (x1-mu1)^2/(sigma1^2) - (2*corr*(x1-mu1)*(x2-mu2))/(sigma1*sigma2) + (x2-mu2)^2/(sigma2^2)
  return(1/(2*pi*sigma1*sigma2*sqrt(1-corr^2)) * exp(-z/(2*(1-corr^2))))
}
#informative prior density value on the grid
info_prior_matrix <- matrix(NA, length(alpha), length(beta))
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    info_prior_matrix[i,j] <- info_prior_density(alpha[i],beta[j], 28.9,2,-0.92,0.05,0)
  }
}

contour(alpha,beta, info_prior_matrix, nlevels = 10, level = ,
        xlim = c(24,33), ylim = c(-1.1,-0.8), main="Contour Plot",
        xlab="alpha",ylab="beta")

# (f) plot the contours and take 1000 draws from
# the joint posterior density of (alpha, beta)
t <- seq(1,10,1)
y_t <- c(24,25,31,31,22,21,26,20,16,22)
alpha <- seq(20,40,0.05)
beta <- seq(-3,2,0.05)

prior_matrix <- matrix(NA, length(alpha), length(beta))
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    if(alpha[i] >= 0 & alpha[i] + 12 * beta[j] >= 0){
      prior_matrix[i,j] <- 1
    }else{
      prior_matrix[i,j] <- 0 
    }
  }
}

likelihood_density <- function(alpha,beta,t,y_t){
  tem <- 1
  for(i in 1:10){
    tem <- tem * exp(-(alpha + beta*t[i])) * (alpha + beta*t[i])^y_t[i]/factorial(y_t[i])
  }
  return(tem)
}

posterior_matrix <- matrix(NA, length(alpha), length(beta))
for(i in 1:length(alpha)){
  for(j in 1:length(beta)){
    posterior_matrix[i,j] <- prior_matrix[i,j] * likelihood_density(alpha[i],beta[j],t,y_t)
  }
}
posterior_matrix <- posterior_matrix/sum(posterior_matrix)
contour(alpha,beta,posterior_matrix, nlevels = 10, level = ,
        xlim = c(20,40), ylim = c(-3,2), main="Contour Plot",
        xlab="alpha",ylab="beta")

#draw 1000 samples from posterior distributions
#compute the marginal posterior distribution of mu
mar_pos_alpha <- rowSums(posterior_matrix)

#draw 1000 samples of alpha from the marginal posterior distribution of alpha
sampleDist <- function(n, x, prob){ 
  return(sample(x, n, replace = T, prob))
}
sample_alpha <- sampleDist(1000, alpha, mar_pos_alpha)

#draw 1000 samples of beta from the marginal posterior distribution of beta given alpha
for(i in 1:length(alpha)){
  posterior_matrix[i, ] <- posterior_matrix[i, ]/mar_pos_alpha[i]
}
sample_beta <- numeric()
for(i in 1:length(sample_alpha)){
  k <- (sample_alpha[i]-20)/0.05 + 1
  sample_beta <- c(sample_beta, sampleDist(1, beta, posterior_matrix[k,]))
}

#For each of the sampled alpha and beta, add a uniform random jitter centered
#at zero with a width equal to the spacing of the sampling grid.
sample_alpha <- sample_alpha + runif(1000,-0.05/2,0.05/2)
sample_beta <- sample_beta + runif(1000,-0.05/2,0.05/2)

# (g) plot a histogram of the posterior density
# for expected number of fatal accidents in 1986
hist(sample_alpha + sample_beta*11, 
     main="expected number of fatal accidents in 1986",
     xlab = "number of fatal accidents", yaxt = "n", 
     breaks = seq(5,35,0.5), cex = 2)

# (h)95% predective interval for the 
# number of fatal accidents in 1986
predictive <- numeric()
for(i in 1:1000){
  predictive <- c(predictive, rpois(1,sample_alpha[i]+11*sample_beta[i]))
}
interval <- sort(predictive)[c(25,975)]
sprintf("95%% posterior interval for the number of fatal accidents in 1986 is [%f, %f]",
        interval[1], interval[2])
