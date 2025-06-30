library(Matrix)
library(mvtnorm)
library(sn)
library(mgcv)
library(ldr)

### Matrix Distance Function ###
mat_dist <- function(B1, B2) { 
                P1 <- B1 %*% solve(t(B1) %*% B1) %*% t(B1)
                P2 <- B2 %*% solve(t(B2) %*% B2) %*% t(B2)
                P_diff <- P1 - P2
                sqrt(sum(diag(P_diff %*% P_diff)))
}
################################

beta1 <- c(2, 1, 1, 0, 0, 0, 0)
beta2 <- c(1, 0, 0, 0, 1, 1, 2)
p <- length(beta1)
B <- cbind(beta1, beta2)
Bs <- qr.Q(qr(B))

n <- 1000
mu_X <- rep(0, p)
Sigma_X <- 0.3^abs(outer(1:p, 1:p, "-"))
eta1 <- c(2, 2, 2, 2, 2, 2, 2)
eta2 <- c(-2.5, -2, 0)
X <- rmvnorm(n, mean = mu_X, sigma = Sigma_X)
eps <- rnorm(n)
y <- -2 + (X %*% beta1 + 10)^2/36 * atan(X %*% beta2) + eps * 2/3
s_y <- cbind(1, y, y^2) %*% eta2
s <- as.numeric(X %*% eta1 + s_y)
prbs <- pnorm(s)
idx0 <- rbinom(n, size = 1, prob = prbs)
idx <- which(idx0 > 0)
length(idx)/n
index <- rep("Mis", n)
index[idx] <- "Obs"
index <- factor(index, levels = c("Mis", "Obs"))
d <- 2
X_obs <- X[idx, ]; y_obs <- y[idx]
X_mis <- X[-idx, ]; y_mis <- y[-idx]
dir1_obs <- X_obs %*% Bs[, 1]; dir1_mis <- X_mis %*% Bs[, 1]
dir2_obs <- X_obs %*% Bs[, 2]; dir2_mis <- X_mis %*% Bs[, 2]
dir1 <- X %*% Bs[, 1]; dir2 <- X %*% Bs[, 2]
dev.new(height = 4.5, width = 5.4, nuit = "in")
op <- par(mfrow = c(1, 1))
plot(y ~ dir1, type = "n")
points(y_obs~dir1_obs, pch = 20, col = 4, cex = 0.8)
points(y_mis~dir1_mis, pch = 20, col = 2, cex = 0.8)
par(op)
dev.new(height = 4.5, width = 5.4, nuit = "in")
op <- par(mfrow = c(1, 1))
plot(y ~ dir2, type = "n")
points(y_obs~dir2_obs, pch = 20, col = 4, cex = 0.8)
points(y_mis~dir2_mis, pch = 20, col = 2, cex = 0.8)
par(op)

z <- y
d <- 2
Fy <- bF_o(y, case = "pdisc", degree = 3, nslices = 3, center = TRUE, scale = FALSE)
Fz <- bF_s(z, case = "pdisc", degree = 0, nslices = 12, center = FALSE, scale = FALSE)

res_try <- LSM_IR(X_obs, X_mis, y, Fy, Fz, index, d = d, 
                  psiA = NULL, csA = NULL, sigmasA = NULL, method = "MaxSkew", 
                  maxit = 2e04, eps = 1e-06, graph = TRUE)

res0x4[[1]]
mat_dist(B[, 1:d], res0x4[[1]])


