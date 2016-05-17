#NAME: Sheng Zhang, UNI:sz2553, Class:STAT W4400, Assignment Number: HW #3
#Question 1.Boosting(70 points)

#read data into R
X <- read.table("Desktop/Statistical Machine Learning/uspsdata/uspsdata.txt")
label <- read.table("Desktop/Statistical Machine Learning/uspsdata/uspscl.txt")
y <- label[,1]

# function <train > for AdaBoost implementation
train <- function(X, w, y){
  n <- dim(X)[1]
  l <- dim(X)[2]
  m <- rep(0, l)
  theta <- rep(0, l)
  cost <- rep(0, l)
  # find optimal theta for every dimension j
  for (j in 1:l){
    # sort the samples along dimension j
    indx <- order(X[,j])
    x_j <- X[indx ,j]
    
    # using a cumulative sum , count the weight when progressively
    # shifting the threshold to the right
    w_cum <- cumsum(w[indx] * y[indx])
    
    # handle multiple occurrences of same x_j value : threshold
    # point must not lie between elements of same value
    
    w_cum[duplicated(x_j)==1] <- NA
    # find the optimum threshold and classify accordingly
    maxIndx <- min(which(abs(w_cum)==max(abs(w_cum), na.rm= TRUE)))
    m[j] <- (w_cum[maxIndx] < 0)*2 - 1
    theta[j] <- x_j[maxIndx]
    c <- ((x_j > theta [j])*2 - 1)* m[j]
    cost[j] <- w %*% (c != y)
  }
  
  # determine optimum dimension , theta and comparison mode
  j_star <- min(which(cost==min(cost)))
  pars <- list(j = j_star, theta = theta[j_star], m = m[j_star])
  return(pars)
}

classify <- function(X, pars){
  labels <- (2*(X[, pars$j] > pars$theta) - 1) * pars$m
  return(labels)
}

agg_class <- function(X, alpha, allPars){
  n <- dim(X)[1]
  B <- length(alpha)
  labels <- matrix(0, n, B)
  # determine labeling for each base classifier
  for(b in 1:B){
    labels[,b] <- classify(X, allPars[[b]])
  }
  # weight classifier response with respective alpha coefficient
  labels <- labels %*% alpha
  c_hat <- sign(labels)
  return(c_hat)
}

AdaBoost <- function(X,y,B){
  n <- dim(X)[1]
  w <- rep(1/n, n)
  alpha <- rep(0, B)
  allPars <- rep(list(list()) ,B)
  
  # boost base classifiers
  for(b in 1:B){
    
    #train weak learner
    allPars[[b]] <- train(X,w,y)
    
    #compute error
    misslabeled <- (y != classify (X, allPars[[b]]))
    error <- (w %*% misslabeled/sum(w))[1]
    
    #compute voting weight
    alpha [b] <- log((1-error)/error)
    
    #recompute weights
    w <- w*exp(alpha[b]* misslabeled)
  }
  return (list(allPars = allPars , alpha = alpha))
}

B_max <- 30
n_cross_validation <- 5

testErrorRate <- matrix(0, B_max , n_cross_validation)
trainErrorRate <- matrix(0, B_max , n_cross_validation)
n <- dim(X)[1]

testErrorRate <- matrix (0, nrow =B_max , ncol =nCV)
trainErrorRate <- matrix (0, nrow =B_max , ncol =nCV )

for (i in 1:n_cross_validation){
  # randomly split data in training and test half
  p <- sample.int(n)
  trainIndx <- p[1:round(n/2)]
  testIndx <- p[-(1:round(n/ 2))]
  ada <- AdaBoost(X[trainIndx,], y[trainIndx], B_max)
  allPars <- ada$allPars
  alpha <- ada$alpha
  # determine error rate , depending on the number of base classifier
  for(B in 1:B_max){
    c_hat_test <- agg_class(X[testIndx,], alpha[1:B], allPars[1:B])
    testErrorRate[B,i] <- mean(y[testIndx] != c_hat_test)
    c_hat_train<- agg_class(X[trainIndx,], alpha[1:B], allPars[1:B])
    trainErrorRate[B,i] <- mean(y[trainIndx] != c_hat_train)
  }
}

B_max <- 50
n_cross_validation <- 5

testErrorRate <- matrix(0, B_max , n_cross_validation)
trainErrorRate <- matrix(0, B_max , n_cross_validation)
n <- dim(X)[1]

#split our data into 5 parts for cross-validation
p <- sample.int(n)
Index <- matrix(0,200/n_cross_validation,n_cross_validation)
Index[,1] <- p[1:40] 
Index[,2] <- p[41:80]
Index[,3] <- p[81:120] 
Index[,4] <- p[121:160] 
Index[,5] <- p[161:200] 

#training and testing
for(B in 1:B_max){
  for(i in 1:n_cross_validation){
    trainIndex <- vector('numeric',0)
    for(j in 1:n_cross_validation){
      if(j != i){
        trainIndex <- c(trainIndex, Index[,j])
      }
    }
    testIndex <- Index[,i]
    Ada <- AdaBoost(X[trainIndex,], y[trainIndex], B_max)
    allPars <- Ada$allPars
    alpha <- Ada$alpha
    # determine error rate , depending on the number of base classifier
    c_hat_test <- agg_class(X[testIndex,], alpha[1:B], allPars[1:B])
    testErrorRate[B,i] <- mean(y[testIndex] != c_hat_test)
    c_hat_train <- agg_class(X[trainIndex,], alpha[1:B], allPars[1:B])
    trainErrorRate[B,i] <- mean(y[trainIndex] != c_hat_train)
  }  
}

# plot results
matplot(trainErrorRate , type ="l",lty =1: n_cross_validation , main =" training error ",
        xlab ="number of weak learners",ylab ="error rate",ylim =c(0,0.5))

matplot(testErrorRate , type ="l",lty =1: n_cross_validation , main =" test error ",
        xlab ="number of weak learners",ylab ="error rate",ylim =c(0,0.5))
