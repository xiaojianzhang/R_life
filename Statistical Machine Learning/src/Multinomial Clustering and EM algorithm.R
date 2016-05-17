#NAME: Sheng Zhang, UNI:sz2553, Class:STAT W4400, Assignment Number: HW #4
#Problem 3.Multinomial Clustering(50 points)

#Load data in R:
H <- matrix(readBin(file("/Users/zhangsheng/Desktop/Statistical Machine Learning/histograms.bin", "rb"), "double", 640000), 40000, 16)
dim_H <- dim(H)

# add a samll constant(1e-3) to the input histograms in case of numerical problems
H[which(H==0)] <- H[which(H==0)] + 1e-3

MultinomialEM <- function(H,K,tau){
  #Choose K clusters at random and normalize each:
  t <- H[sample(1:dim_H[1], size = K, replace = FALSE), ]
  t <- normalize(t)
  c <- rep(1,K)/K
  # do iterations:
  #(a)E-step:
  phi <- exp(H %*% t(t))
  A <- Generate_A(c,phi)
  
  #(b)M-step:
  c <- colSums(A)/dim(H)[1]
  b <- t(A) %*% H
  t <- normalize(b)
  A_old <- matrix(0, dim(A)[1], dim(A)[2])
  while(norm(A-A_old, type="o") >= tau){
    A_old <- A
    #(a)E-step:
    phi <- exp(H %*% t(t))
    A <- Generate_A(c,phi)
    
    #(b)M-step:
    c <- colSums(A)/dim(H)[1]
    b <- t(A) %*% H
    t <- normalize(b)
    
  }
  m <- numeric(dim(H)[1])
  for(i in 1:dim(H)[1]){
    m[i] <- which.max(A[i, ])
  }
  return(m)
}

normalize <- function(x){
  for(i in 1:dim(x)[1]){
    x[i, ] <- x[i, ] / sum(x[i, ])
  }
  return(x)
}

Generate_A <- function(c,phi){
  denominator <- phi %*% c
  for(i in 1:dim(phi)[2]){
    phi[, i] <- phi[, i] * c[i]
  }
  for(i in 1:dim(phi)[1]){
    phi[i, ] <- phi[i, ] / denominator[i] 
  }
  return(phi)
}

image_plot <- function(m){
  M <- matrix(m, nrow = 200, byrow = TRUE)
  im <- t(M[nrow(M):1, ])
  image(x=1:200, y=1:200, im, axes=FALSE, col = grey((0:256)/256), xlab = "", ylab = "")
}

for(i in 3:5){
  m <- MultinomialEM(H,i,1e-4)
  image_plot(m)
}
