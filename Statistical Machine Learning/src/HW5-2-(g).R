#NAME: Sheng Zhang, UNI:sz2553, Class:STAT W4400, Assignment Number: HW #5
#Problem 2(g)

#Generate 256 exp(theta=1) data:
x = rexp(256, rate=1)

#hyperparameter:
alpha0 = 2
beta0 = 0.2

theta = seq(0,4,0.01)

posterior <- function(theta, n){
  post <- dgamma(theta, shape= alpha0 + n,rate= beta0+sum(x[1:n]))
  return(post)
}

post <- mapply(posterior,theta,4)
plot(theta, post, 'l', ylim = c(0,7), col=1, ylab='p(theta|x_1:n)')
k<-2
for(n in c(8,16,256)){
  post <- mapply(posterior,theta,n)
  lines(theta, post, col=k)
  k <- k+1
}
legend('topright',c('n=4', 'n=8', 'n=16', 'n=256'), lty=1, col=c(1,2,3,4),bty='n',cex=1)
