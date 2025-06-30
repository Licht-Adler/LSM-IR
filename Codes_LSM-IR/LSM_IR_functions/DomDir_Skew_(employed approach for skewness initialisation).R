library(Matrix)
library(sn)
library(mgcv)

#######################################################################################################
################ Skewness Maximisation method for dominant direction (Loperfido, 2018) ################
#######################################################################################################
DomDir_Skew <- function(data) { 
    if (!is.matrix(data)) stop("The object must be of matrix type.\n")
    n <- nrow(data); p <- ncol(data)
    if (n <= p) stop("The sample size must be larger than the number of random variables.\n")
    ### Obtain the directional skewness vector for bivariate case ###
    ### Best rank-1 approximation of symmetric binary tensors (LCMV algorithm) ###
    ### Givens rotation for the case of (2 * 1)-vector ###
    ### ref. Loperfido (2018) ###
    ### ref. de Lathauwer, Comon, de Moor & Vandewalle (1996) ###
    ### ref. de Lathauwer, de Moor & Vandewalle (2000b) ###
    MaxSkewBiv <- function(x1, x2) { 
        data2 <- cbind(x1, x2)
        n <- length(x1)
        skew_tmp <- mean((scale(x2))^3) / ((1 - 1/n)^1.5)
        if (abs(cor(x1, x2)) > 0.99 & skew_tmp >= 0) { 
            v1_Biv <- c(0, 1)
            sgn_adj <- 1
            gamma1_Biv <- skew_tmp
        } else if (abs(cor(x1, x2)) > 0.99 & skew_tmp < 0) { 
            v1_Biv <- c(0, -1)
            sgn_adj <- -1
            gamma1_Biv <- -skew_tmp
        } else if (abs(cor(x1, x2)) <= 0.99) { 
            Sigma_Biv <- cov(data2) * (1 - 1/n)
            G0_Biv <- pdm_power(Sigma_Biv, -0.5)
            mu_Biv <- colMeans(data2)
            Z_Biv <- t(G0_Biv %*% (t(data2) - mu_Biv))
            z1 <- Z_Biv[, 1]
            z2 <- Z_Biv[, 2]
            a111 <- mean(z1^3)
            a112 <- mean(z1^2 * z2)
            a122 <- mean(z1 * z2^2)
            a222 <- mean(z2^3)
            poly_coef <- c(-a122, a222 - 2 * a112, 2 * a122 - a111, a112)
            if (round(sqrt(sum(poly_coef^2)), 8) == 0) poly_coef <- c(1, 0, 0, 0)
            slv <- polyroot(poly_coef)
            slv_real <- Re(slv[round(Im(slv), 8) == 0])
            if (length(slv_real) == 1) { 
                v1_Biv <- c(slv_real, 1) / sqrt(1 + slv_real^2)
                gamma1_Biv <- mean((scale(Z_Biv %*% v1_Biv))^3) / ((1 - 1/n)^1.5)
            } else if (length(slv_real) > 1) { 
                v1_Biv_cdd <- sapply(slv_real, function(s) { c(s, 1) / sqrt(1 + s^2) }, 
                                     simplify = TRUE)
                skew_Biv_cdd <- apply(v1_Biv_cdd, 2, function(v) { mean((scale(Z_Biv %*% v))^3) })
                v1_Biv <- v1_Biv_cdd[, which.max(abs(skew_Biv_cdd))]
                gamma1_Biv <- skew_Biv_cdd[which.max(abs(skew_Biv_cdd))] / ((1 - 1/n)^1.5)
            }
            sgn_adj <- 1
            if (gamma1_Biv < 0) { 
                v1_Biv <- -v1_Biv
                gamma1_Biv <- -gamma1_Biv
                sgn_adj <- -1
            }
        }
        res <- list(v1_Biv = v1_Biv, sgn_adj = sgn_adj, gamma1_Biv = gamma1_Biv)
        return(res)
    }
    ### Whitening data ###
    Sigma <- cov(data) * (1 - 1/n)
    G0 <- pdm_power(Sigma, -0.5)
    mu <- colMeans(data)
    Z <- t(G0 %*% (t(data) - mu))
    ### Obtain directional skewness vector when p == 2 ###
    if (p == 2) { 
        z1 <- Z[, 1]
        z2 <- Z[, 2]
        res_2 <- MaxSkewBiv(z1, z2)
        dir1_MS <- as.vector(G0 %*% res_2$v1_Biv)
        return(dir1_MS)
    }
    ### Initialisation based on HOSVD when p > 2 ###
    ### ref. de Lathauwer, de Moor & Vandewalle (2000a) ###
    M3_z <- matrix(rowMeans(apply(Z, 1, function(z) { z %x% t(z) %x% z })), ncol = p)
    v1_ini <- as.vector(svd(M3_z)$v[, 1])
    gamma1_ini <- mean((scale(Z %*% v1_ini))^3) / ((1 - 1/n)^1.5)
    if (gamma1_ini < 0) { 
        v1_ini <- -v1_ini
        gamma1_ini <- -gamma1_ini
    }
    ### Generalise the Bivariate Skewness Maximisation to the p-variate case ... ###
    ### by Using the Jacobi-type iteration algorithm ###
    ### ref. Loperfido (2018) ###
    ### ref. Golub & van Loan (1989) ###
    iter_max <- 50; iter <- 0
    gamma1_iter <- NULL; gamma1_prv <- 0
    tol <- 1e-08; diff <- Inf
    v1_iter <- v1_ini
    while (iter <= iter_max & abs(diff) >= tol) { 
        iter <- iter + 1
        gamma1_iter <- cbind(gamma1_iter, rep(0, p + 1))
        gamma1_iter[1, iter] <- gamma1_prv
        for (j in 1:p) { 
            y1 <- Z[, j]
            v1_iter <- v1_iter / sqrt(sum(v1_iter^2))
            v1_iter_old <- v1_iter[j]
            v1_iter[j] <- 0
            y2 <- as.vector(Z %*% v1_iter)
            res_iter <- MaxSkewBiv(y1, y2)
            if (res_iter$gamma1_Biv >= gamma1_iter[j, iter]) { 
                gamma1_iter[1 + j, iter] <- res_iter$gamma1_Biv
                v1_iter[j] <- res_iter$v1_Biv[1] / res_iter$v1_Biv[2] * res_iter$sgn_adj
            } else { 
                gamma1_iter[1 + j, iter] <- gamma1_iter[j, iter]
                v1_iter[j] <- v1_iter_old
            }
        }
        v1_iter <- v1_iter / sqrt(sum(v1_iter^2))
        if (iter == 1) gamma1_iter[1, 1] <- gamma1_iter[2, 1]
        gamma1_prv <- gamma1_iter[1 + p, iter]
        diff <- gamma1_prv - gamma1_iter[1, iter]
    }
    ### Obtain the final directional skewness vector ###
    v1_MaxSkew <- v1_iter
    gamma1_MaxSkew <- mean((scale(Z %*% v1_MaxSkew))^3) / ((1 - 1/n)^1.5)
    if (gamma1_MaxSkew < 0) { 
        v1_MaxSkew <- -v1_MaxSkew
        gamma1_MaxSkew <- -gamma1_MaxSkew
    }
    dir1_MaxSkew <- as.vector(G0 %*% v1_MaxSkew)
    return(dir1_MaxSkew)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
