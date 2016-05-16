library("nleqslv")
y <- c(43, 44, 45, 46.5, 47.5)
# find the posterior mode of theta
# x:theta
sys_eqs <- function(x){
  z <- numeric()
  
  z <- (x-y[1])/(1+(y[1]-x)^2) + (x-y[2])/(1+(y[2]-x)^2) + (x-y[3])/(1+(y[3]-x)^2) + (x-y[4])/(1+(y[4]-x)^2) + (x-y[5])/(1+(y[5]-x)^2)
  
  return(z)
}
xstart <- c(44)
x_end <- (nleqslv(xstart, sys_eqs, control = list(btol=1e-7)))[[1]]

sec_der <- 0
for(i in 1:5){
  sec_der <- sec_der + 2*((y[i]-x_end)^2-1)/(1+(y[i]-x_end)^2)^2
}

I_theta <- 1/(-sec_der)

# plot the approximate normal density 
par(mfrow=c(1,2))
x <- seq(-40, 50, 0.01)
z <- dnorm(x, 44.862, sqrt(0.727))
plot(x, z,type="l", 
     main= "normal approximation",
     xlab= "theta", ylab = "",
     xaxs = "i", yaxs ="i",
     yaxt="n",bty="n",
     col= "blue", xlim = c(40,50), cex=2)

# compute the density in Exercise 2.11
theta <- seq(0,100,0.01)
prior_theta <- vector(length = length(theta))

for(i in 1:length(theta)){
  prior_theta[i] <- 1/100
}

tem <- numeric()
likelihood <- function(a){
  for(i in 1:5){
    tem <- c(tem, 1/(1+(y[i]-a)^2))
  }
  c <- prod(tem)
  return(c)
}

posterior_theta <- numeric()
for(i in 1:length(prior_theta)){
  posterior_theta <- c(posterior_theta, prior_theta[i]*likelihood(theta[i]))
}
posterior_theta <- posterior_theta/sum(posterior_theta)
plot(theta, posterior_theta, type="l", 
     main= "exact posterior density",
     xlab= "theta", ylab = "",
     xaxs = "i", yaxs ="i",
     yaxt="n",bty="n",
     col= "blue", xlim = c(40,50), cex=2)