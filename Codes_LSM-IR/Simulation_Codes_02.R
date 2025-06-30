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

beta1 <- c(2, rep(0, 7), c(1, 1, 0, 0, 0, 0))
beta2 <- c(1, rep(0, 7), c(0, 0, 0, 1, 1, 2))
p <- length(beta1)
B <- cbind(beta1, beta2)
N <- 500
n <- 1000
link_index <- 1:2
intcpts <- c(-0.4, 1.8, 
             -3.0, -0.7, 
             -5.4, -3, 
             -7.8, -5.3, 
             -10.4, -7.8)

pair <- cbind(intcpts, rep(1:2, 5))
pair_list <- lapply(apply(pair, 1, "list"), "unlist")

require(mvtnorm)
mu_X <- rep(0, p)
Sigma_X <- 0.3^abs(outer(1:p, 1:p, "-"))
eta1 <- rep(2, p)
dp1 <- list(mu_X = mu_X, Sigma_X = Sigma_X, eta1 = eta1)

links02 <- as.list(1:6)
links02[[1]] <- function(n, intcpt, dp1) { 
                   epsilon0 <- rnorm(n)
                   mu_X <- dp1$mu_X; Sigma_X <- dp1$Sigma_X; eta1 <- dp1$eta1
                   X0 <- rmvnorm(n, mean = mu_X, sigma = Sigma_X)
                   y0 <- -2 + 2 * sign(X0 %*% beta1) * abs(X0 %*% beta1)^(1/3) + epsilon0
                   eta2 <- c(intcpt, 0, 0.5)
                   s_xy <- as.numeric(X0 %*% eta1 + cbind(1, y0, y0^2) %*% eta2)
                   prbs <- pnorm(s_xy)
                   idx0 <- rbinom(n, size = 1, prob = prbs)
                   idx <- which(idx0 > 0)
                   index <- rep("Mis", n)
                   index[idx] <- "Obs"
                   index <- factor(index, levels = c("Mis", "Obs"))
                   d <- 1
                   data <- list(X0 = X0, y0 = y0, index = index, idx = idx, 
                                d = d, psiA = eta1, csA = eta2)
                   return(data)
}

links02[[2]] <- function(n, intcpt, dp1) { 
                   epsilon0 <- 2/3 * rnorm(n)
                   mu_X <- dp1$mu_X; Sigma_X <- dp1$Sigma_X; eta1 <- dp1$eta1
                   X0 <- rmvnorm(n, mean = mu_X, sigma = Sigma_X)
                   y0 <- -2 + (X0 %*% beta1 + 10)^2/36 * atan(X0 %*% beta2) + epsilon0
                   eta2 <- c(intcpt, -2)
                   s_xy <- as.numeric(X0 %*% eta1 + cbind(1, y0) %*% eta2)
                   prbs <- pnorm(s_xy)
                   idx0 <- rbinom(n, size = 1, prob = prbs)
                   idx <- which(idx0 > 0)
                   index <- rep("Mis", n)
                   index[idx] <- "Obs"
                   index <- factor(index, levels = c("Mis", "Obs"))
                   d <- 2
                   data <- list(X0 = X0, y0 = y0, index = index, idx = idx, 
                                d = d, psiA = eta1, csA = eta2)
                   return(data)
}

library(parallel)
cl <- makeCluster(detectCores() - 1)
parLSMIR02 <- function(pair, N) { 

    clusterExport(cl, c("beta1", "beta2", "B", "p"), envir = globalenv())
    clusterExport(cl, c("mu_X", "Sigma_X", "eta1", "dp1", "links02"), envir = globalenv())
    clusterExport(cl, c("bF_o", "bF_s", "pdm_power", "zeta1", "zeta2", "mode_HW", "DomDir_Skew", 
                        "JADE_Skew", "init_outeq", "init_slceq", "init_est", "LSMIR_init_estep", 
                         "LSMIR_init_cmstep", "LSMIR_init_ecm", "LSMIR_main_estep", "LSMIR_main_cmstep", 
                         "LSMIR_main_ecm", "LSM_IR"), 
                  envir = globalenv())
    
    clusterEvalQ(cl, lapply(c("Matrix", "MASS", "GrassmannOptim", "sn", 
                              "mvtnorm", "mgcv", "splines", "nlme"), 
                            require, character.only = TRUE))

    res02 <- parSapply(cl, 1:N, function(x, pair) { 
                           intcpt <- pair[1]
                           link <- pair[2]
                           n <- 1000
                           d <- ifelse(link > 2, 2, 1)
                           data <- links02[[link]](n, intcpt, dp1)
                           y0 <- data$y0; X0 <- data$X0
                           index <- data$index; idx <- data$idx
                           X_obs <- X0[idx, ]; X_mis <- X0[-idx, ]
                           y_obs <- y0[idx]; y_mis <- y0[-idx]

                           Fy <- bF_o(y0, case = "pdisc", degree = 3, nslices = 3, center = TRUE)
                           Fz1 <- bF_s(y0, case = "pdisc", degree = 0, nslices = 12, center = FALSE)
                           
                           data$Fy <- Fy; data$Fz1 <- Fz1

                         proc21 <- system.time(
                           tmp21 <- try(LSM_IR(X_obs, X_mis, as.numeric(y0), Fy, Fz1, index, 
                                                     d = d, psiA = NULL, csA = NULL, sigmasA = NULL, 
                                                     method = "MaxSkew", maxit = 2e4, eps = 1e-5, 
                                                     graph = FALSE), 
                                        silent = FALSE)
                        )
                          tmp21[[length(tmp21) + 1]] <- proc21

                           tmp2 <- list()
                           tmp2$PDS <- tmp21
                           tmp2$Data_list <- data

                           return(tmp2)
                       }, pair)
    return(res02)
}
stopCluster(cl)

N <- 500
res_Sim02 <- vector(mode = "list", nrow(pair))
(T1 <- Sys.time())
cl <- makeCluster(detectCores() - 2)
for (i in 1:nrow(pair)) { 
    res9070I16[[i]] <- parLSMIR02(pair[i, ], N)
    print(pair[i, ])
    print(Sys.time())
}
stopCluster(cl)
(T2 <- Sys.time())
T2 - T1

