library(Matrix)
library(sn)
library(mgcv)

#######################################################################################################
################### JADE method for skewness direction (Cardoso & Souloumiac, 2013) ###################
#######################################################################################################
JADE_Skew <- function(data, eps = 1e-06, maxiter = 100) { 
    X <- na.fail(data)
    if (!all(sapply(X, is.numeric))) 
        stop("The data matrix must be of numeric type .\n", call. = FALSE)
    X <- as.matrix(X)
    n <- nrow(X); p <- ncol(X)
    ### Mean values removal ###
    mu_X <- colMeans(X)
    X_c <- t(t(X) - mu_X)
    ### Whitening and projection onto signal subspace ###
    eig0 <- eigen(crossprod(X_c) / (n - 1), symmetric = TRUE)
    lambda0 <- eig0$values
    U0 <- eig0$vectors
    W0 <- U0 %*% diag(sqrt(1/lambda0))
    Z <- X_c %*% W0
    ### Estimation of the cumulant matrices ####
    pos0 <- cbind(rep(1:p, times = p), rep(1:p, each = p))
    pos <- pos0[which(pos0[, 1] - pos0[, 2] >= 0), ]
    pos_list <- lapply(apply(pos, 1, "list"), "unlist")
    CMs <- lapply(pos_list, function(stell, data) { 
                      n <- nrow(data); m <- ncol(data)
                      constant1 <- ifelse(diff(stell) == 0, 1, 0)
                      constant2 <- ifelse(diff(stell) == 0, 1, sqrt(2))
                      Eij <- matrix(0, nrow = m, ncol = m)
                      Eij[stell[1], stell[2]] <- 1
                      Kij <- t(((Z[, stell[1]] * Z[, stell[2]]) %*% t(rep(1/n, p))) * Z) %*% Z
                      Qij <- constant2 * (Kij - Eij - t(Eij) - constant1 * diag(p))
                      Qij
                  }, data = Z)
    CMs <- do.call("rbind", CMs)
    ### Joint diagonalization of the cumulant matrices ###
        ### starting values of the Jacobi-type iteration ###
    nrow_CMs <- nrow(CMs)
    num_CM <- nrow_CMs / p
    V_rot <- diag(p)
    wieder <- TRUE
    iter <- 0
        ### joint diagonalization proper ###
    while (wieder == TRUE) { 
        iter <- iter + 1
        wieder <- FALSE
        for (i in 1:(p - 1)) { 
            for (j in (i + 1):p) { 
                I_i <- seq(i, nrow_CMs, p)
                I_j <- seq(j, nrow_CMs, p)
            ### computation of Givens angle ###
                g <- cbind(CMs[I_i, i] - CMs[I_j, j], CMs[I_j, i] + CMs[I_i, j])
                gg <- crossprod(g)
                t_diag <- gg[1, 1] - gg[2, 2]
                t_off <- gg[1, 2] + gg[2, 1]
                theta <- 1/2 * atan2(t_off, t_diag + sqrt(t_diag^2 + t_off^2))
                cos_theta <- cos(theta)
                sin_theta <- sin(theta)
            ### Jacobi update ###
                if (abs(theta) > eps) { 
                    wieder <- TRUE
                    row_I_i <- CMs[I_i, ]
                    row_I_j <- CMs[I_j, ]
                    CMs[I_i, ] <- cos_theta * row_I_i + sin_theta * row_I_j
                    CMs[I_j, ] <- cos_theta * row_I_j - sin_theta * row_I_i
                    col_i <- CMs[, i]
                    col_j <- CMs[, j]
                    CMs[, i] <- cos_theta * col_i + sin_theta * col_j
                    CMs[, j] <- cos_theta * col_j - sin_theta * col_i
                    temp_V_rot <- V_rot[i, ]
                    V_rot[i, ] <- cos_theta * V_rot[i, ] + sin_theta * V_rot[j, ]
                    V_rot[j, ] <- cos_theta * V_rot[j, ] - sin_theta * temp_V_rot
                }
            }
        }
        if (iter > maxiter) wieder <- FALSE
    }
    V_z <- W0 %*% t(V_rot)
    IC_z <- X_c %*% V_z
    ### Choose the vector maximizing the skewness of the indpendent components respectively ###
    ### Select the two largest absolute skewness values ###
    ### If their difference is very small, choose the one with larger kurtosis ###
    skew_z <- colMeans(IC_z^3, na.rm = TRUE)
    kurt_z <- colMeans(IC_z^4, na.rm = TRUE) - 3
    sgn_sk <- sign(skew_z)
    idx_biv <- order(abs(skew_z), decreasing = TRUE)[1:2]
    if (abs(diff(abs(skew_z[idx_biv]))) < 5e-03 ) { 
        idx <- idx_biv[which.max(kurt_z[idx_biv])]
    } else { idx <- idx_biv[1] }
    ### Obtain the desired direction related to the skewness vector ###
    dir1 <- as.vector(V_z[, idx]) * sgn_sk[idx]
    return(dir1)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
