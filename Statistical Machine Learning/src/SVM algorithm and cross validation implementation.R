#SVM algorithm and cross validation implementation
#NAME: Sheng Zhang UNI:sz2553

library("e1071")
library("scatterplot3d")
data <- read.table("Desktop/Statistical Machine Learning/uspsdata/uspsdata.txt")
label <- read.table("Desktop/Statistical Machine Learning/uspsdata/uspscl.txt")
image_label <- label[,1]
image_data <- vector()
for(i in 1:200){
  image_data <- rbind(image_data, unlist(data[i,]))
}

#randomly select 20% of the data
sample <- sort(sample(seq(1,200,1),40,replace = F,rep(1/200,200)))
test_set <- vector()
test_label <- vector()
train_set <- vector()
train_label <- vector()
for(i in 1:200){
  if(is.element(i,sample)){
    test_set <- rbind(test_set, image_data[i,])
    test_label <- c(test_label, image_label[i])
  }else{
    train_set <- rbind(train_set, image_data[i,])
    train_label <- c(train_label, image_label[i])
  }
}

## classification mode
#model <- svm(train_set, train_label)
tune.linear <- tune.svm(train_set, as.factor(train_label), kernel="linear",cost=2^seq(-10,0,0.5))
print(tune.linear$performances)
model_linear <- tune.linear$best.model
pred <- predict(model_linear,train_set)
table(pred,train_label)
pred <- predict(model_linear,test_set)
table(pred,as.factor(test_label))

tune.non <- tune(svm, train_set, as.factor(train_label), kernel="radial",
                 ranges=list(cost=seq(0.1,3,0.2), gamma=seq(1e-3,1e-2,1e-3)))
print(tune.non)
model_non <- tune.non$best.model
pred <- predict(model_non,train_set)
table(pred,train_label)
pred <- predict(model_non,test_set)
table(pred,as.factor(test_label))

#plot1
xdata <- 2^seq(-10,0,0.5)
ydata <- tune.linear$performances$error
plot(xdata, ydata, 
     main= "margin parameter v.s. misclassification rate(linear case)",
     xlab= "margin parameter cost", ylab = "incorrect rate",
     col= "blue", cex=0.5)

#plot2
xdata <- tune.non$performances$cost
ydata <- tune.non$performances$gamma
zdata <- tune.non$performances$error
scatterplot3d(xdata, ydata, zdata, highlight.3d=TRUE, col.axis="blue",
              main = "margin parameter + the kernel bandwidth v.s. misclassification rate(nonlinear case)",
              xlab = "margin parameter cost",
              ylab = "the kernel bandwidth gamma", 
              zlab = "incorrect rate",
              col.grid="lightblue", pch=20)

