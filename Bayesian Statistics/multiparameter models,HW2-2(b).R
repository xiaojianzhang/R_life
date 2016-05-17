# draw 1000 samples from two t distributions, respectively
mu_c <- rt(1000,31) * (0.24/sqrt(32)) + 1.013
mu_t <- rt(1000,35) * (0.20/sqrt(36)) + 1.173

# use samples to approximate
diff <- mu_t - mu_c

# plot the histogram of the diff
hist(diff, xlab = "mu_t - mu_c", yaxt = "n", breaks = seq(-0.1,0.4,0.01), cex = 2)

# 95% posterior interval for diff
interval <- sort(diff)[c(25, 976)]
sprintf("95%% posterior interval for mu_t - mu_c is [%f, %f]", interval[1], interval[2])