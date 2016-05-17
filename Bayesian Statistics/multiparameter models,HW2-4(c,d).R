library("foreign")
library("nleqslv")
y <- c(16/(58+16), 9/(9+90), 10/(10+48), 13/(13+57), 19/(19+103), 20/(20+57), 18/(18+86), 17/(17+112), 35/(35+273), 55/(55+64))
z <- c(12/(12+113), 1/(1+18), 2/(2+14), 4/(4+44), 9/(9+208), 7/(7+67), 9/(9+29), 8/(154+8))
# grid from 0.01 to 0.99 with gap 0.01
x1 <- seq(0.01, 0.99, 0.01)
x2 <- seq(0.01, 0.25, 0.01)
grid_matrix <- matrix(NA,99,25)
for(i in 1:length(x1)){
  mu_sigma <- x1[i]*(1-x1[i])
  for(j in 1:length(x2)){
    if(x2[j] < mu_sigma){
      grid_matrix[i,j] <- 1/mu_sigma
    }else{
      grid_matrix[i,j] <- 0
    }
  }
}
# solve for alpha, beta using grid points produced above
# x[1]:alpha, x[2]:beta
sys_eqs <- function(x){
  y <- numeric(2)
  
  y[1] <- x[1] / (x[1] + x[2]) - mu
  y[2] <- (x[1]*x[2])/((x[1]+x[2]) * (x[1]+x[2]+1)) - sigma_square
  
  return(y)
}
prior_matrix <- grid_matrix
alpha_matrix <- matrix(NA,99,25)
beta_matrix <- matrix(NA,99,25)
for(i in 1:length(x1)){
  for(j in 1:length(x2)){
    if(grid_matrix[i,j] > 0){
      mu <- x1[i]
      sigma_square <- x2[j]
      xstart <- c(4.5,20)
      x_end <- (nleqslv(xstart, sys_eqs, control = list(btol=1e-6)))[[1]]
      alpha_matrix[i,j] <- x_end[1]
      beta_matrix[i,j] <- x_end[2]
    }else{
      alpha_matrix[i,j] <- 0
      beta_matrix[i,j] <- 0
    }
  }
}

# produce posterior density of theta_y and theta_z, respectively
dbeta_y <- numeric()
dbeta_z <- numeric()

likelihood_y <- function(a,b){
  for(i in 1:10){
    dbeta_y <- c(dbeta_y, dbeta(y[i], a, b))
  }
  c <- prod(dbeta_y)
  return(c)
}

likelihood_z <- function(a,b){
  for(i in 1:8){
    dbeta_z <- c(dbeta_z, dbeta(z[i], a, b))
  }
  c <- prod(dbeta_z)
  return(c)
}

posterior_y_matrix <- matrix(NA,99,25)
posterior_z_matrix <- matrix(NA,99,25)

for(i in 1:99){
  for(j in 1:25){
    if(prior_matrix[i,j] > 0){
      posterior_z_matrix[i,j] <- prior_matrix[i,j]*likelihood_z(alpha_matrix[i,j],beta_matrix[i,j])
    }else{
      posterior_z_matrix[i,j] <- 0
    }
  }
}

for(i in 1:99){
  for(j in 1:25){
    if(prior_matrix[i,j] > 0){
      posterior_y_matrix[i,j] <- prior_matrix[i,j]*likelihood_y(alpha_matrix[i,j],beta_matrix[i,j])
    }else{
      posterior_y_matrix[i,j] <- 0
    }
  }
}

# normalize posterior
posterior_y_matrix <- posterior_y_matrix/sum(posterior_y_matrix)
posterior_z_matrix <- posterior_z_matrix/sum(posterior_z_matrix)

#draw 1000 samples from posterior distributions
#compute the marginal posterior distribution of mu
mar_pos_mu_y <- rowSums(posterior_y_matrix)
mar_pos_mu_z <- rowSums(posterior_z_matrix)

#draw 1000 samples of mu from the marginal posterior distribution of mu
sampleDist <- function(n, x, prob){ 
  return(sample(x, n, replace = T, prob))
}
sample_mu_y <- sampleDist(1000, x1, mar_pos_mu_y)
sample_mu_z <- sampleDist(1000, x1, mar_pos_mu_z)
diff_mu_y_mu_z <- sample_mu_y - sample_mu_z
hist(diff_mu_y_mu_z, xlab = "mu_y - mu_z", yaxt = "n", breaks = seq(-0.2,0.4,0.01), cex = 2)
