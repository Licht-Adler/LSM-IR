library(Matrix)
library(sn)
library(mgcv)

#######################################################################################################
################################ ECM-PX Algorithm for the LSM-IR Model ################################
#######################################################################################################
LSMIR_main_ecm <- function(X_obs, X_mis = NULL, y, Foy, Fsy, index, d, ini_param = list(), psiA = NULL, 
                         csA = NULL, sigmasA = NULL, maxit = 1e04, eps = 1e-04, graph = FALSE) { 
    ### Function vech ###
    vech <- function(x) { 
        if (is.vector(x)) {
            if (length(x) == 1) 
                return(x)
            else stop("vech undefined for vectors")
        }
        else if (is.matrix(x)) {
            d <- ncol(x)
            if (d != nrow(x)) 
                stop("vech only defined for square matrices")
            vechx <- vector()
            for (j in 1:d) vechx <- c(vechx, x[j:d, j])
            return(vechx)
        }
    }
    ### Preparation ###
    Fy <- Foy; Fz <- Fsy
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
    if (length(psiA) < 0 | length(psiA) > p) 
        stop("Incorrect length of the vector 'psiA'.\n", call. = FALSE)
    if (!is.null(psiA) & !is.vector(psiA)) 
        stop("Incorrect format of the vector 'psiA'.\n", call. = FALSE)
    if (!is.null(psiA) & any(is.character(psiA))) 
        stop("Incorrect setup of the vector 'psiA'.\n", call. = FALSE)
    if (!is.null(csA) & length(csA) != ncol(Fz)) 
        stop("Incorrect length of the vector 'csA'.\n", call. = FALSE)
    if (!is.null(csA) & !is.numeric(csA)) 
        stop("Incorrect format of the vector 'csA'.\n", call. = FALSE)
    if (!is.null(csA) & any(is.na(csA))) 
        stop("Incorrect setup of the vector 'csA'.\n", call. = FALSE)
    K_par <- p * (p + 5) / 2 + (p + r1) * d + r2 + 1
    if (K_par >= n) 
        stop("Sample size too small or number of parameters too large.\n", call. = FALSE)
    ### Initialising ###
    log_Liks <- log_Lik_ac <- log_Lik_cd <- rep(c(-Inf), 1, maxit + 1)
    pars_cnvg <- rep(c(Inf), 1, maxit)
    par_vec0 <- rep(0, K_par)
    if (!is.list(ini_param)) 
        stop("The initial values must bu organised in a list\n", call. = FALSE)
    if (length(ini_param) < 6) 
        stop("Incorrect settings of initial values\n", call. = FALSE)
    ini_param_list <- ini_param
    Es1_list <- LSMIR_main_estep(ini_param_list, X_obs, Fy, Fz, index)
    param1_list <- LSMIR_main_cmstep(Es1_list, ini_param_list, X_obs, Fy, Fz, index, d, 
                                   psiA, csA, sigmasA)
    log_Liks[1] <- log_Lik_ac[1] <- param1_list$log_Lik_ac
    log_Lik_cd[1] <- param1_list$log_Lik_cd
    par_vec0[1:p] <- as.vector(param1_list$xi)
    par_vec0[(1 + p):(p + p * d)] <- as.vector(param1_list$Gamma)
    par_vec0[(1 + p + p * d):(p + (p + r1) * d)] <- as.vector(param1_list$Co)
    par_vec0[(1 + p + (p + r1) * d):(K_par - r2 - p - 1)] <- vech(param1_list$Delta)
    par_vec0[(K_par - r2 - p):(K_par - r2 - 1)] <- as.vector(param1_list$psi)
    par_vec0[(K_par - r2):(K_par - 1)] <- as.vector(param1_list$cs)
    par_vec0[K_par] <- as.numeric(param1_list$sigma_s)
    lLinf_new <- log_Liks[1]
    par_vec_new <- par_vec0
    pre_param_list <- param1_list
    for (it in 1:maxit) { 
        Es_list <- LSMIR_main_estep(pre_param_list, X_obs, Fy, Fz, index)
        param_list <- LSMIR_main_cmstep(Es_list, pre_param_list, X_obs, Fy, Fz, index, d, 
                                      psiA, csA, sigmasA)
        log_Liks[1 + it] <- log_Lik_ac[1 + it] <- param_list$log_Lik_ac
        log_Lik_cd[1 + it] <- param_list$log_Lik_cd
        par_vec0[1:p] <- as.vector(param_list$xi)
        par_vec0[(1 + p):(p + p * d)] <- as.vector(param_list$Gamma)
        par_vec0[(1 + p + p * d):(p + (p + r1) * d)] <- as.vector(param_list$Co)
        par_vec0[(1 + p + (p + r1) * d):(K_par - r2 - p - 1)] <- vech(param_list$Delta)
        par_vec0[(K_par - r2 - p):(K_par - r2 - 1)] <- as.vector(param_list$psi)
        par_vec0[(K_par - r2):(K_par - 1)] <- as.vector(param_list$cs)
        par_vec0[K_par] <- as.numeric(param_list$sigma_s)
        par_vec_old <- par_vec_new
        par_vec_new <- par_vec0
        pars_cnvg[it] <- sqrt(sum((par_vec_new - par_vec_old)^2))
        pre_param_list <- param_list
        if (it >= 2) { 
            acc <- (log_Liks[it + 1] - log_Liks[it]) / 
                       (log_Liks[it] - log_Liks[it - 1])
            lLinf_old <- lLinf_new
            lLinf_new <- try(log_Liks[it] + 1 / (1 - acc) * 
                                 (log_Liks[it + 1] - log_Liks[it]))
            change_par <- pars_cnvg[it]
            if (abs(change_par) < eps | is.na(change_par)) break
        }
    }
    xi_hat <- param_list$xi
    Gamma_hat <- param_list$Gamma
    Delta_hat <- param_list$Delta
    Co_hat <- param_list$Co
    psi_hat <- param_list$psi
    cs_hat <- param_list$cs
    sigma_s_hat <- param_list$sigma_s
    sigma0 <- param_list$sigma0
    B_hat <- qr.Q(qr(solve(Delta_hat) %*% Gamma_hat))
    colnames(B_hat) <- paste("Dir", 1:d, sep = "")
    comp_aux <- param_list$comp_aux
    tau_o <- param_list$tau_o
    tau_m <- param_list$tau_m
    if (is.null(psiA)) { 
        df_psi <- p
    } else { 
        df_psi <- p - length(which(is.na(psiA)))
    }
    if (is.null(csA)) { df_cs <- r2 } else { df_cs <- 0 }
    log_Lik_ac_hat <- param_list$log_Lik_ac
    log_Lik_cd_hat <- param_list$log_Lik_cd
    nbprms <- p * (p + 3) / 2 + d * (p + r1 - d) + df_psi + df_cs + 1
    AIC <- (log_Lik_ac_hat - nbprms) * (-2)
    BIC <- (log_Lik_ac_hat - log(n) / 2 * nbprms) * (-2)
    muX_m <- Es_list$muX_m
    VX_m <- Es_list$VX_m
    muS_m <- Es_list$muS_m
    vS_m <- Es_list$vS_m
    muS_o <- Es_list$muS_o
    vS_o <- Es_list$vS_o
    kappaXS_m <- Es_list$kappaXS_m
    X_imp <- rbind(X_obs, t(muX_m))
    X_hat0 <- t(xi_hat + Gamma_hat %*% t(Co_hat) %*% t(Fy))
    X_hat <- X_hat0[c(which(index == labels[2]), which(index == labels[1])), ]
    s_imp <- c(muS_o, muS_m)
    s_hat0 <- as.numeric(Fz %*% cs_hat + X_hat0 %*% psi_hat)
    s_hat <- s_hat0[c(which(index == labels[2]), which(index == labels[1]))]
    ### Plots drawing part ###
    if (graph) { 
        dirs_o <- matrix(X_obs %*% B_hat, nrow = n1)
        y_obs <- y[which(index == labels[2])]
        if (!is.null(X_mis)) { 
            if (any(is.na(X_mis))) { 
                X_mis_t <- t(X_mis)
                pos_na <- which(is.na(X_mis_t))
                X_mis_t[pos_na] <- muX_m[pos_na]
                X_mis <- t(X_mis_t)
            }
            dirs_m <- matrix(X_mis %*% B_hat, nrow = n2)
            y_mis <- y[which(index == labels[1])]
        } else { 
            X_mis <- t(muX_m)
            dirs_m <- matrix(X_mis %*% B_hat, nrow = n2)
            y_mis <- y[which(index == labels[1])]
        }
        dev.new(height = 400, width = 400, unit = "in")
        if (d == 1) { 
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
            plot(log_Lik_ac[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
                 xlab = "Iterations", ylab = "Actual-Constrained Log-Likelihood Space Observations")
            par(op)
        } else if (d > 1) { 
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
            plot(log_Lik_ac[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
                 xlab = "Iterations", ylab = "Actual-Constrained Log-Likelihood Space Observations")
            par(op)
        }
        dev.new(height = 400, width = 400, unit = "in")
        op <- par(mfrow = c(2, 1), mar = c(2.2, 2.2, 1.1, 1.1), 
                      mgp = c(1.2, 0.4, 0), cex = 0.8)
        plot(log_Lik_cd[1:(1 + it)], type = 'p', pch = 20, col = "red", cex = 0.8, 
             xlab = "Iterations", ylab = "Complete_Data Log-Likelihood Space Observations")
        plot(log10(pars_cnvg[1:it]), type = 'p', pch = 20, col = "red", cex = 0.8, xlab = "Iterations", 
             ylab = "Logarithm of L-1 Distance between New/Old Parameters Estimated")
        par(op)
    }
    res_list <- list(B_hat = B_hat, xi_hat = xi_hat, Gamma_hat = Gamma_hat, 
                     Co_hat = Co_hat, Delta_hat = Delta_hat, psi_hat = psi_hat, 
                     cs_hat = cs_hat, sigma_s_hat = sigma_s_hat, 
                     iteration = it + 1, 
                     minimum_change = min(abs(diff(log_Lik_ac[1:(1 + it)]))), 
                     latest_increase = log_Lik_ac[1 + it] - log_Lik_ac[it], 
                     log_Lik_ac = log_Lik_ac_hat, AIC = AIC, BIC = BIC, 
                     log_Lik_cd =  log_Lik_cd_hat, pars_cnvg_fin = pars_cnvg[it], 
                     muX_m = muX_m, VX_m = VX_m, muS_m = muS_m, vS_m = vS_m, muS_o = muS_o, 
                     vS_o = vS_o, kappaXS_m = kappaXS_m, comp_aux = comp_aux, 
                     tau_o = tau_o, tau_m = tau_m, sigma0 = sigma0)

    return(res_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
