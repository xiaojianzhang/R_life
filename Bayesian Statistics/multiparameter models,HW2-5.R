library("boot")
x <- c(-0.86, -0.30, -0.05, 0.73)
n <- c(5, 5, 5, 5)
y <- c(0, 1, 3, 5)

x1 <- seq(-5,10,0.05)
x2 <- seq(-10,40,0.05)
#define bivariate gaussian density function
bivariate_gaussian_density <- function(x1,x2,mu1,sigma1,mu2,sigma2,corr){
  z <- (x1-mu1)^2/(sigma1^2) - (2*corr*(x1-mu1)*(x2-mu2))/(sigma1*sigma2) + (x2-mu2)^2/(sigma2^2)
  return(1/(2*pi*sigma1*sigma2*sqrt(1-corr^2)) * exp(-z/(2*(1-corr^2))))
}
# prior distribution
prior_matrix <- matrix(NA, length(x1), length(x2))
for(i in 1:length(x1)){
  for(j in 1:length(x2)){
    prior_matrix[i,j] <- bivariate_gaussian_density(x1[i],x2[j], 0,2,10,10,0.5)
  }
}

tem <- numeric()

likelihood <- function(a,b,n,x){
  for(i in 1:4){
    tem <- c(tem, inv.logit(a+b*x[i])^y[i] * (1-inv.logit(a+b*x[i]))^(n[i]-y[i]))
  }
  c <- prod(tem)
  return(c)
}

posterior_matrix <- matrix(NA, length(x1), length(x2))
for(i in 1:length(x1)){
  for(j in 1:length(x2)){
    posterior_matrix[i,j] <- prior_matrix[i,j] * likelihood(x1[i],x2[j],n,x)
  }
}
#contour plot
par(mfrow=c(1,2))
contour(x1,x2,posterior_matrix, nlevels = 10, level = ,
        xlim = c(-5,10), ylim = c(-10,40), main="Contour Plot",
        xlab="alpha",ylab="beta")

# normalize posterior
posterior_matrix <- posterior_matrix/sum(posterior_matrix)

#draw 1000 samples from posterior distributions
#compute the marginal posterior distribution of mu
mar_pos_alpha <- rowSums(posterior_matrix)

#draw 1000 samples of alpha from the marginal posterior distribution of alpha
sampleDist <- function(n, x, prob){ 
  return(sample(x, n, replace = T, prob))
}
sample_alpha <- sampleDist(1000, x1, mar_pos_alpha)

#draw 1000 samples of beta from the marginal posterior distribution of beta given alpha
for(i in 1:length(x1)){
  posterior_matrix[i, ] <- posterior_matrix[i, ]/mar_pos_alpha[i]
}
sample_beta <- numeric()
for(i in 1:length(sample_alpha)){
  k <- (sample_alpha[i]-(-5))/0.05 + 1
  sample_beta <- c(sample_beta, sampleDist(1, x2, posterior_matrix[k,]))
}

#For each of the sampled alpha and beta, add a uniform random jitter centered
#at zero with a width equal to the spacing of the sampling grid.
sample_alpha <- sample_alpha + runif(1000,-0.05/2,0.05/2)
sample_beta <- sample_beta + runif(1000,-0.05/2,0.05/2)
plot(sample_alpha, sample_beta,
     main = "scatter plot",
     type= "p",
     xlab= "alpha",
     ylab= "beta",
     xlim = c(-5,10),
     ylim = c(-10,40),
     pch = ".")

#the posterior distribution of the LD50
par(mfrow=c(1,1))
hist(-sample_alpha/sample_beta, main="LD50 histogram", xlab = "LD50", yaxt = "n", breaks = seq(-0.7,0.7,0.01), cex = 2)

