#HW6-3:Chapter 14 Exercise 1.
#question(a)
y1 <- c (5.0, 13.0, 7.2, 6.8, 12.8, 5.8, 9.5, 6.0, 3.8, 14.3, 1.8, 6.9, 4.7, 9.5)
y2 <- c (0.9, 12.9, 2.6, 3.5, 26.6, 1.5, 13.0, 8.8, 19.5, 2.5, 9.0, 13.1, 3.6, 6.9)
y3 <- c (14.3, 6.9, 7.6, 9.8, 2.6, 43.5, 4.9, 3.5, 4.8, 5.6, 3.5, 3.9, 6.7)
basement1 <- c(1,1,1,1,1,0,1,1,1,0,1,1,1,1)
basement2 <- c(0,1,1,0,1,1,1,1,1,0,1,1,1,0)
basement3 <- c(1,0,1,0,1,1,1,1,1,1,1,1,1)

make_indicators <- function (x){
  ux <- unique(x)
  matrix1 <- matrix(x, nrow=length(x), ncol=length(ux))
  matrix2 <- matrix(ux, nrow=length(x), ncol=length(ux), byrow=TRUE)
  (matrix1==matrix2)*1}

counties <- rep(1:3,c(length(y1),length(y2),length(y3)))
y <- c(y1,y2,y3)
x <- cbind(c(basement1,basement2,basement3), make_indicators(counties))
ls_out <- lsfit (x, log(y), intercept=F)
lsd <- ls.diag(ls_out)

nsims <- 10000
n <- nrow(x)
k <- ncol(x)
sigma <- rep(NA, nsims)
beta <- array(NA, c(nsims, k))
for (i in 1:nsims){
  sigma[i] <- lsd$std.dev*sqrt((n-k)/rchisq(1,n-k))
  beta[i,] <- ls_out$coef + (sigma[i]/lsd$std.dev)*lsd$std.err*t(chol(lsd$corr))%*%rnorm(k)}

output <- exp (cbind (beta[,2], beta[,1]+beta[,2], beta[,3],
                      beta[,1] + beta[,3], beta[,4], beta[,1] + beta[,4], beta[,1], sigma))

for(i in 1:ncol(output)){
  print (round(quantile(output[,i],c(.25,.5,.75)),1))
}
  
#question(b)
theta <- rbeta(nsims, 13, 3)
b <- rbinom(nsims, 1, theta)
logy_rep <- rnorm (nsims, beta[,3] + b*beta[,1], sigma)
y_rep <- exp(logy_rep)
print (round(quantile(y_rep,c(.025,.25,.5,.75,.975)),1))
hist (y_rep[y_rep<40], breaks=0:40,
      main = 'the posterior predictive distribution for its radon measurement',
      xlab="radon measurement (new house)", freq = FALSE,cex=2)

