library(Matrix)
library(sn)
library(mgcv)

#######################################################################################################
############################## ECM Function for Initialisation Procedure ##############################
#######################################################################################################
LSMIR_init_ecm <- function(X_obs, X_mis = NULL, Fy, Fz, index, d, ini_param = list(), 
                       maxit = 1e04, eps = 1e-04, graph = FALSE) { 
    ### Preparation ###
    if (!is.matrix(X_obs)) 
        stop("The sample data of the observed covariates must be of a matrix form.\n", 
             call. = FALSE)
    if (!is.matrix(Fy)) 
        stop("The data of the outcome response bases must be of a matrix form.\n", 
             call. = FALSE)
    if (!is.matrix(Fz)) 
        stop("The data of the selection covariates bases must be of a matrix form.\n", 
             call. = FALSE)
    if (nrow(Fy) != nrow(Fz)) 
        stop("Dimensions are not matched between outcome and selection equations.\n", 
             call. = FALSE)
    if (is.factor(index)) { labels <- levels(index) } else { labels <- unique(index) }
    if (length(labels) != 2) 
        stop("Incorrect structure of the selection/missing index.\n", call. = FALSE)
    Fy_m <- as.matrix(Fy[which(index == labels[1]), ])
    Fy_o <- as.matrix(Fy[which(index == labels[2]), ])
    Fz_m <- as.matrix(Fz[which(index == labels[1]), ])
    Fz_o <- as.matrix(Fz[which(index == labels[2]), ])
    Ff_m <- cbind(Fy_m, Fz_m); Ff_o <- cbind(Fy_o, Fz_o)
    fy_bar <- colMeans(Fy); fz_bar <- colMeans(Fz)
    if (length(index) != nrow(Fy)) { 
        stop("Dimensions are not matched between missingness index and response sample.\n", 
             call. = FALSE)
    } else { n <- nrow(Fy) }
    n1x <- nrow(X_obs)
    if (n1x != nrow(Fy_o)) 
        stop("Dimensions are not matched between observed covariates and response.\n", 
             call. = FALSE)
    n1 <- nrow(Fy_o); n2 <- n - n1
    p <- ncol(X_obs); r1 <- ncol(Fy); r2 <- ncol(Fz)
    if (p <= 1 | r1 < 1 | is.null(r2)) 
        stop("The dimensions of observed covariates and response bases are too small.\n", 
             call. = FALSE)
    if (qr(X_obs)$rank != p) 
        stop("The observed covariates matrix must be of full column rank.\n", 
             call. = FALSE)
    if (d < 1 | d > min(p, r1)) 
        stop("Incorrect choice of the structural dimension or 'r1'.\n", call. = FALSE)
    ### Initialising ###
    log_Liks <- rep(c(-Inf), 1, maxit + 1)
    if (!is.list(ini_param)) 
        stop("The initial values must bu organised in a list\n", call. = FALSE)
    if (length(ini_param) < 6) 
        stop("Incorrect settings of initial values\n", call. = FALSE)
    ini_param_list <- ini_param
    Es1_list <- LSMIR_init_estep(ini_param_list, X_obs, Fy, Fz, index)
    param1_list <- LSMIR_init_cmstep(Es1_list, ini_param_list, X_obs, Fy, Fz, index, d)
    log_Liks[1] <- param1_list$log_Lik
    lLinf_new <- log_Liks[1]
    pre_param_list <- param1_list
    for (it in 1:maxit) { 
        Es_list <- LSMIR_init_estep(pre_param_list, X_obs, Fy, Fz, index)
        param_list <- LSMIR_init_cmstep(Es_list, pre_param_list, X_obs, Fy, Fz, index, d)
        log_Liks[1 + it] <- param_list$log_Lik
        pre_param_list <- param_list
        if (it >= 2) { 
            acc <- (log_Liks[it + 1] - log_Liks[it]) / 
                       (log_Liks[it] - log_Liks[it - 1])
            lLinf_old <- lLinf_new
            lLinf_new <- try(log_Liks[it] + 1 / (1 - acc) * 
                                 (log_Liks[it + 1] - log_Liks[it]))
            if (abs(lLinf_new - lLinf_old) < eps | is.na(lLinf_new)) break
        }
    }
    xi_hat <- param_list$xi
    Gamma_hat <- param_list$Gamma
    D_hat <- param_list$D
    C_hat <- param_list$C
    h_hat <- param_list$h
    cs_hat <- param_list$cs
    Delta_hat <- D_hat + tcrossprod(h_hat)
    B_hat <- qr.Q(qr(solve(Delta_hat) %*% Gamma_hat))
    colnames(B_hat) <- paste("Dir", 1:d, sep = "")
    omega_hat <- sqrt(diag(Delta_hat))
    alpha_hat <- omega_hat * as.vector(solve(D_hat) %*% h_hat) / 
                     sqrt(1 + as.numeric(t(h_hat) %*% solve(D_hat) %*% h_hat))
    log_Lik_hat <- param_list$log_Lik
    nbprms <- p * (p + 5) / 2 + d * (p + r1 - d) + r2
    AIC <- log_Lik_hat - nbprms
    BIC <- log_Lik_hat - log(n) / 2 * nbprms
    if (graph) { 
        dirs_o <- matrix(X_obs %*% B_hat, nrow = n1)
        y_obs <- y[which(index == labels[2])]
        if (!is.null(X_mis)) { 
            dirs_m <- matrix(X_mis %*% B_hat, nrow = n2)
            y_mis <- y[which(index == labels[1])]
        }
        dev.new(height = 400, width = 400, unit = "in")
        if (d == 1) { 
          if (is.null(X_mis)) { 
            op <- par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1.1, 1.1), 
                      mgp = c(1.2, 0.4, 0), cex = 0.8)
            plot(y_obs ~ dirs_o[, 1], xlab = "1st Direction: Observed (Blue)", 
                 ylab = "Response: Observed (Blue)", type = "p", main = "", 
                 pch = 20, col = "blue")
            model01 <- lowess(dirs_o[, 1], y_obs, f = 2/3, iter = 5)
            lines(model01, col = 1, lty = 2, lwd = 2)
            plot(log_Liks[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
                 xlab = "Iterations", ylab = "Log-Likelihood Space Observations")
            par(op)
          } else { 
            op <- par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1.1, 1.1), 
                      mgp = c(1.2, 0.4, 0), cex = 0.8)
            yf <- c(y_obs, y_mis)
            dirs <- rbind(dirs_o, dirs_m)
            plot(yf ~ dirs[, 1], type = "n", xlab = "1st Direction: Obsd.(Blue) & Miss.(Red)", 
                 ylab = "Response: Observed (Blue) & Missing (Red)")
            points(y_obs ~ dirs_o[, 1], pch = 20, col = "blue")
            points(y_mis ~ dirs_m[, 1], pch = 20, col = "red")
            model01 <- lowess(dirs[, 1], yf, f = 2/3, iter = 5)
            lines(model01, col = 1, lty = 2, lwd = 2)
            plot(log_Liks[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
                 xlab = "Iterations", ylab = "Log-Likelihood Space Observations")
            par(op)
          }
        } else if (d > 1) { 
          if (is.null(X_mis)) { 
            op <- par(mar = c(2.2, 2.2, 1.1, 1.1), mgp = c(1.2, 0.4, 0), cex = 0.8)
            layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE))
            plot(y_obs ~ dirs_o[, 1], xlab = "1st Direction: Obsd.(Blue)", 
                 ylab = "Response: Observed (Blue)", type = "p", main = "", 
                 pch = 20, col = "blue")
            model01 <- lowess(dirs_o[, 1], y_obs, f = 2/3, iter = 5)
            lines(model01, col = 1, lty = 2, lwd = 2)
            plot(y_obs ~ dirs_o[, 2], xlab = "2nd Direction: Obsd.(Blue)", 
                 ylab = "Response: Observed (Blue)", type = "p", main = "", 
                 pch = 20, col = "blue")
            model02 <- lowess(dirs_o[, 2], y_obs, f = 2/3, iter = 5)
            lines(model02, col = 1, lty = 2, lwd = 2)
            plot(log_Liks[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
                 xlab = "Iterations", ylab = "Log-Likelihood Space Observations")
            par(op)
          } else { 
            op <- par(mar = c(2.2, 2.2, 1.1, 1.1), mgp = c(1.2, 0.4, 0), cex = 0.8)
            layout(matrix(c(1, 2, 3, 3), nrow = 2, ncol = 2, byrow = TRUE))
            yf <- c(y_obs, y_mis)
            dirs <- rbind(dirs_o, dirs_m)
            plot(yf ~ dirs[, 1], type = "n", xlab = "1st Direction: Obsd.(Blue) & Miss.(Red)", 
                 ylab = "Response: Observed (Blue) & Missing (Red)")
            points(y_obs ~ dirs_o[, 1], pch = 20, col = "blue")
            points(y_mis ~ dirs_m[, 1], pch = 20, col = "red")
            model01 <- lowess(dirs[, 1], yf, f = 2/3, iter = 5)
            lines(model01, col = 1, lty = 2, lwd = 2)
            plot(yf ~ dirs[, 2], type = "n", xlab = "2nd Direction: Obsd.(Blue) & Miss.(Red)", 
                 ylab = "Response: Observed (Blue) & Missing (Red)")
            points(y_obs ~ dirs_o[, 2], pch = 20, col = "blue")
            points(y_mis ~ dirs_m[, 2], pch = 20, col = "red")
            model02 <- lowess(dirs[, 2], yf, f = 2/3, iter = 5)
            lines(model02, col = 1, lty = 2, lwd = 2)
            plot(log_Liks[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
                 xlab = "Iterations", ylab = "Log-Likelihood Space Observations")
            par(op)
          }
        }
    }
    res_list <- list(B_hat = B_hat, Gamma_hat = Gamma_hat, xi_hat = xi_hat, 
                     Delta_hat = Delta_hat, h_hat = h_hat, C_hat = C_hat, 
                     cs_hat = cs_hat, D_hat = D_hat, alpha_hat = alpha_hat, 
                     iteration = it + 1, 
                     minimum_change = min(abs(diff(log_Liks[1:(1 + it)]))), 
                     latest_increase = log_Liks[1 + it] - log_Liks[it], 
                     log_Lik = log_Lik_hat, AIC = AIC, BIC = BIC)
    return(res_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
