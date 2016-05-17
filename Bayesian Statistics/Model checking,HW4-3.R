ydata<-c(1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0,0)
sample_theta <- rbeta(1000,8,14)
T <- numeric(1000)
for(i in 1:1000){
  num <- 0
  y <- c()
  while(num < 13){
    y_rep <- rbinom(1,1,sample_theta[i])
    y <- c(y,y_rep)
    if(y_rep == 0){
      num <- num + 1
    }
  }
  T[i] <- sum(abs(y[2:length(y)]-y[1:length(y)-1]))
  
}
hist(T, xlab="T(y^rep)", yaxt="n",
     breaks=seq(0, 25,1), cex=2)