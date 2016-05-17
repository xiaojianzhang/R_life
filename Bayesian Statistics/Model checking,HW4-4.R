A <- 1e+5
y_bar <- 5.1
T <- numeric(1000)
for(i in 1:1000){
  y_rep <- rnorm(100,mean = y_bar, sd = sqrt(1+0.01))
  T[i] <- max(abs(y_rep))
}

hist(T, xlab="T(y_rep)", yaxt="n",
     breaks=seq(6, 10, 0.1), cex=2)

p_value <- sum(T>=8.1)/1000
