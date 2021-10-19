Rcpp::sourceCpp("Functions_Merge.cpp")

###############################################################
## Requires functions from Functions_Merge.cpp
library(mvtnorm)

N <- 5000
J <- 7
set.seed(2021)
x0 <- rep(1, N)
x1 <- rnorm(N)
x2 <- rnorm(N, sd = 0.4)
X <- matrix(c(x0, x1, x2), ncol = 3)
beta_true <- c(1.8, -0.34, -0.65)
z_true <- as.vector(X %*% beta_true) + rnorm(N, 0, 1)
cutpoints_true <- c(-1000000, 0, 0.75, 1.25, 1.8, 2.5, 2.9, 1000000)
hist(z_true, breaks = 50)
abline(v = cutpoints_true, lty = "dashed", col = "blue", lwd = 2)
y <- rep(0, N)
for(j in 1:length(cutpoints_true)){
  for(i in 1:N){
    if(z_true[i] >= cutpoints_true[j]){ 
      y[i] <- j
    }
  }
}
hist(y)

# Initiate parameters for sampler
cutpoints0 <- c(0, 0.8, 1.3, 2.0, 2.5, 3.0, 10000)
cutpoints0[J] <- 10 ^ 9
beta0 <- c(1.65, -.29, -.6)
B0 <- matrix(c(c(0.4,0,0), 
               c(0,0.3,0),
               c(0,0,0.2)), nrow = 3)

# function arguments - temporary
beta <- beta0
cutpoints <- cutpoints0
cp_vect <- cutpoints[2:(J-1)]
beta <- as.vector(beta)
mus <- as.vector(X %*% beta)

trace <- list(cutpoints, rep(0, N), beta0, log_post_cpp(cutpoints[2:(J-1)], y, mus), 0)

M <- 2000

#trace <- list(cutpoints, rep(0, N), beta0, log_post_cpp(cutpoints[2:(J-1)], y, mus), 0)
#for(m in 1:M){
#  beta <- trace[[5*m - 2]]
#  mus <- as.vector(X %*% as.vector(beta))
#  cutpoints <- trace[[5*m-4]]
#  trace <- c(trace, Algorithm_2_cpp(y, mus, cutpoints, X, B0, beta0))
#  ## trace <- c(trace, Alg_2_cpp(y, X, trace[[5*m - 2]], trace[[5*m-4]]))
#}

start_time <- Sys.time()
trace <- ordinal_Sampler_cpp(trace, M, y, X, B0, beta0)
end_time <- Sys.time()
end_time - start_time

########################################################################
## Plot trace vectors: ################

# MH Step Acceptance Points
plot(sapply(trace[c(FALSE, FALSE, FALSE, FALSE, TRUE)], "[[", 1), 
     main ="MH_Accept")

# Cutpoints
plot(sapply(trace[c(TRUE, FALSE, FALSE, FALSE, FALSE)], "[[", 2), type = "n", 
     main ="First Free Cutpoint")
lines(sapply(trace[c(TRUE, FALSE, FALSE, FALSE, FALSE)], "[[", 2))
abline(h = cutpoints_true[3], lty = "dashed", lwd = 3, col = "blue")
plot(sapply(trace[c(TRUE, FALSE, FALSE, FALSE, FALSE)], "[[", 3), type = "n", 
     main ="Second Free Cutpoint")
lines(sapply(trace[c(TRUE, FALSE, FALSE, FALSE, FALSE)], "[[", 3))
abline(h = cutpoints_true[4], lty = "dashed", lwd = 3, col = "blue")
plot(sapply(trace[c(TRUE, FALSE, FALSE, FALSE, FALSE)], "[[", 6), type = "n", 
     main ="Last Free Cutpoint")
lines(sapply(trace[c(TRUE, FALSE, FALSE, FALSE, FALSE)], "[[", 6))
abline(h = cutpoints_true[7], lty = "dashed", lwd = 3, col = "blue")

plot(sapply(trace[c(FALSE, FALSE, TRUE, FALSE, FALSE)], "[[", 1), type = "n", 
     main ="beta0")
lines(sapply(trace[c(FALSE, FALSE, TRUE, FALSE, FALSE)], "[[", 1))
abline(h = beta_true[1], lty = "dashed", lwd = 3, col = "blue")
abline(h = beta0[1], lty = "dashed", lwd = 3, col = "red")

plot(sapply(trace[c(FALSE, FALSE, TRUE, FALSE, FALSE)], "[[", 2), type = "n", 
     main ="beta1")
lines(sapply(trace[c(FALSE, FALSE, TRUE, FALSE, FALSE)], "[[", 2))
abline(h = beta0[2], lty = "dashed", lwd = 3, col = "red")
abline(h = beta_true[2], lty = "dashed", lwd = 3, col = "blue")

plot(sapply(trace[c(FALSE, FALSE, TRUE, FALSE, FALSE)], "[[", 3), type = "n", 
     main ="beta2")
lines(sapply(trace[c(FALSE, FALSE, TRUE, FALSE, FALSE)], "[[", 3))
abline(h = beta_true[3], lty = "dashed", lwd = 3, col = "blue")
abline(h = beta0[3], lty = "dashed", lwd = 3, col = "red")

plot(sapply(trace[c(FALSE, FALSE, FALSE, TRUE, FALSE)], "[[", 1), 
     main ="Log Posterior")
abline(h = log_post_cpp(cutpoints_true[3:(J)], y, as.vector(X %*% beta_true)), 
       lty = "dashed", lwd = 3, col = "blue")


library(ggplot2)
ggplot(data = data.frame(i = seq(1, length(y)), 
                         z = z_true, 
                         y = as.factor(y))) +
  geom_point(aes(x = i, y = z, col = y)) +
  geom_hline(yintercept = cutpoints_true[2:J], lty = "dashed", col = "blue") +
  ggtitle("True Latent Zs and Cutpoints")

ggplot(data = data.frame(i = seq(1, length(y)), 
                         z = trace[[5*M - 3]], 
                         y = as.factor(y))) +
  geom_point(aes(x = i, y = z, col = y)) +
  geom_hline(yintercept = cutpoints_true[2:J], lty = "dashed", col = "blue") +
  geom_hline(yintercept = trace[[5*M-4]][1:(J-1)]) +
  ggtitle("Latent Zs and Sampled/True Cutpoints, Mth iteration")