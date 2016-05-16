#Perceptron algorithm implementation
#NAME: Sheng Zhang UNI:sz2553

#Inputs
#w:  w[1:d] is the normal vector of a hyperplane, 
#    w[d+1] = -c is the negative offset parameter. 
#n: sample size

#Outputs
#S: n by (d+1) sample matrix with last col 1
#y: vector of the associated class labels

fakedata <- function(w, n){

if(! require(MASS))
{
	install.packages("MASS")
}
if(! require(mvtnorm))
{
	install.packages("mvtnorm")
}

require(MASS)
require(mvtnorm)

# obtain dimension
d <- length(w)-1

# compute the offset vector and a Basis consisting of w and its nullspace
offset <- -w[length(w)] * w[1:d] / sum(w[1:d]^2)
Basis <- cbind(Null(w[1:d]), w[1:d])	 

# Create samples, correct for offset, and extend
# rmvnorm(n,mean,sigme) ~ generate n samples from N(0,I) distribution
S <- rmvnorm(n, mean=rep(0,d),sigma = diag(1,d)) %*%  t(Basis) 
S <- S + matrix(rep(offset,n),n,d,byrow=T)
S <- cbind(S,1)

# compute the class assignments
y <- as.vector(sign(S %*% w))

# add corrective factors to points that lie on the hyperplane.
S[y==0,1:d] <- S[y==0,1:d] + runif(1,-0.5,0.5)*10^(-4)
y = as.vector(sign(S %*% w))
return(list(S=S, y=y))

} # end function fakedata

# To evaluate a Perceptron solution
classify <- function(S,z){
  y <- vector()
  for(i in 1:dim(S)[1]){
    if((S[i,] %*% z) > 0){
      y[i] <- 1
    }else{
      y[i] <- -1
    }
  }
  return(y)
}

#Implement the Perceptron training algorithm with alpha(k)=1/k
perceptron <- function(S,y){
  #compute the value of perceptron cost function
  cost_fun <- function(S,y,z){
    tem <- S %*% z
    indicator <- sign(tem) - y
    cost <- 0
    for(i in 1:length(y)){
      if(indicator[i] != 0){
        cost <- cost + abs(tem[i])
      }
    }
    return(cost)
  }
  
  #compute the value of gradient of the cost function
  gradient_cost <- function(S,y,z){ 
    tem <- S %*% z
    indicator <- sign(tem) - y
    gradient <- rep(0, dim(S)[2])
    for(i in 1:length(y)){
      if(indicator[i] != 0){
        gradient <- gradient + (-y[i]) * S[i,]
      }
    }
    return(gradient)
  }
  
  # do interations
  z <- rep(1,dim(S)[2])
  Z_history <- vector()
  for(k in 1:1000){
    if(cost_fun(S,y,z) > 1e-5){
      Z_history <- rbind(Z_history, z)
      z <- z - (1/k) * gradient_cost(S,y,z)
    }else{
      break
    }
  }
  return(list(z=z, Z_history=Z_history))
}

#generate training data by fakedata
w <- c(1,-1, 5)
n <- 100
train_data <- fakedata(w,n)
S_train <- train_data$S
y_train <- train_data$y

train <- perceptron(S_train,y_train)
z <- train$z
Z_history <- train$Z_history

#generate test data by fakedata
test_data <- fakedata(w,n)
S <- test_data$S
y <- test_data$y
#check whether the test data is correctly classified
y_test <- classify(S,z)
error_number <- sum(0.5*abs(y_test - y))
test_correct_rate <- 100 * (1 - error_number/n)

#plot1  
xdata <- seq(-6,2,0.01)
ydata <- (-z[1]/z[2])*xdata - (z[3]/z[2])
plot(xdata, ydata, type="l", 
     main= "test data set and the classifier hyperplane",
     xlab= "x", ylab = "y",
     col= "blue", cex=2)
grid(5,5,col = "gray")
for(i in 1:length(y)){
  if(y[i] == 1){
    points(S[i,1],S[i,2],pch="+")
  }else{
    points(S[i,1],S[i,2],pch="-")
  }
  
}

#plot2
xdata <- seq(-6,2,0.01)
ydata <- (-z[1]/z[2])*xdata - (z[3]/z[2])
plot(xdata, ydata, type="l", 
     main= "training data and visualizing Z_history",
     xlab= "x", ylab = "y",
     col= "blue", cex=2)
grid(5,5,col = "gray")
for(i in 1:length(y_train)){
  if(y_train[i] == 1){
    points(S_train[i,1],S_train[i,2],pch="+")
  }else{
    points(S_train[i,1],S_train[i,2],pch="-")
  }
}

for(i in 1:dim(Z_history)[1]){
  z_his <- Z_history[i,]
  z_his <- z_his/norm(z_his,"2")
  xdata <- c(-4,0)
  ydata <- (-z_his[1]/z_his[2])*xdata - (z_his[3]/z_his[2])
  abline(xdata,ydata)
}