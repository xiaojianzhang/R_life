library("foreign")
library("nleqslv")

pew_research_center <- read.dta("Desktop/Bayesian Statistics/pew_research_center_june_elect_wknd_data.dta")
ElectionResult <- read.csv("Desktop/Bayesian Statistics/2008ElectionResult.csv")
ideo <- pew_research_center[['ideo']]
state_pew <- pew_research_center[["state"]]
state_elec <- ElectionResult[['state']]
vote_Obama_pct <- ElectionResult[['vote_Obama_pct']]
state_pew_table <- table(state_pew)
idex_ideo <- which(ideo %in% "very liberal")
ideo_table <- table(state_pew[idex_ideo])

xdata <- numeric()
ydata <- numeric()
zdata <- vector()
for(i in 1:51){
  if(i != 2 & i != 9 & i != 12){
    ydata <- c(ydata, ideo_table[i] / state_pew_table[i])
    xdata <- c(xdata, vote_Obama_pct[i])
  }
}
for(i in 1:50){
  if(i != 2 & i != 11){
    zdata <- c(zdata, state.abb[i])
  }
}

num_res = vector()
for(i in 1:51){
  if(i != 2 & i != 9 & i != 12){
    num_res <- c(num_res, state_pew_table[i])
  }
}

# Estimating the parameters of the prior distribution, i.e. alpha and beta
# y_j: the number of people who label themselves "very liberal" in the state j
# n_j: the number of respondents in the state j
# From the appendix A we know the mean and variance of the beta-binomial distribution
# E[y_j/n_j] = alpha / (alpha + beta)
# Var[y_j/n_j] = alpha*beta*(alpha+beta+n_j)(alpha+beta)^2*(alpha+beta+1)
# Note: we substitute n_j in variance formula with E[n_j] i.e. the sample average of n_j
y_n_sample <- ydata
E_y_n <- mean(y_n_sample)
V_y_n <- var(y_n_sample)
E_n <- mean(num_res)

# x[1]:alpha, x[2]:beta
sys_eqs <- function(x){
  y <- numeric(2)
  
  y[1] <- x[1] / (x[1] + x[2]) - E_y_n
  y[2] <- (x[1]*x[2]*(x[1] + x[2] + E_n)) / (E_n * (x[1] + x[2])^2 * (x[1] + x[2] + 1)) - V_y_n
  
  return(y)
}

xstart <- c(4.5,100)
x_end <- (nleqslv(xstart, sys_eqs, control = list(btol=1e-6)))[[1]]
alpha <- x_end[1]
beta <- x_end[2]

# posterior density p(theta_j|y_j) ~ Beta(alpha + y_j, beta + n_j - y_j)
# thus, posterior mean E(theta_j|y_j) = (alpha + y_j) / (alpha + beta + n_j)
post_mean <- numeric()
for(i in 1:51){
  if(i != 2 & i != 9 & i != 12){
    post_mean <- c(post_mean, (alpha + ideo_table[i]) / (alpha + beta + state_pew_table[i]))
  }
}

plot(xdata, ydata, 
     main= "proportion liberal in each state vs. Obama vote share",
     xlab= "Obama vote share (in %)",
     ylab= "proportion liberal in each state",
     col= "blue", pch = 20, cex = 0.7, lty = "solid", lwd = 1)

text(xdata, ydata, labels=zdata, cex= 0.5, pos = 4)

plot(xdata, post_mean, 
     main= "posterior mean in each state vs. Obama vote share",
     xlab= "Obama vote share (in %)",
     ylab= "posterior mean in each state",
     col= "blue", pch = 20, cex = 0.7, lty = "solid", lwd = 1)

text(xdata, post_mean, labels=zdata, cex= 0.5, pos = 4)

plot(num_res, ydata, 
     main= "proportion liberal in each state vs. number of respondents",
     xlab= "number of respondents",
     ylab= "proportion liberal in each state",
     col= "blue", pch = 20, cex = 0.7, lty = "solid", lwd = 1)

text(num_res, ydata, labels=zdata, cex= 0.5, pos = 4)

plot(num_res, post_mean, 
     main= "posterior mean in each state vs. number of respondents",
     xlab= "number of respondents",
     ylab= "posterior mean in each state",
     col= "blue", pch = 20, cex = 0.7, lty = "solid", lwd = 1)

text(num_res, post_mean, labels=zdata, cex= 0.5, pos = 4)
