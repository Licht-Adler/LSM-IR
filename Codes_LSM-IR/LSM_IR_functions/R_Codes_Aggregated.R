library(Matrix)
library(sn)
library(mgcv)

#######################################################################################################
################ Generate the basis functions $F_o(y)$ in the inverse outcome equation ################
#######################################################################################################
bF_o <- function(y, case = c("poly", "categ", "fourier", "pdisc", "pcont", "nps", "BS"), 
                 degree = 1, nslices = 3, center = FALSE, scale = FALSE) { 
            require(splines)
            require(mgcv)
            Indicator <- function(H, y) { as.numeric(y %in% H) }
            rightends <- function(n, nslices) { 
                             if ((n - round(n) != 0) | (nslices - round(nslices) != 0))
                                 stop("'n' and 'nslices' must be both integers.\n", call. = FALSE)
                             quotient <- floor(n / nslices)
                             multiples <- (1:nslices) * quotient
                             remainder <- n %% nslices
                             ends <- rep(0, nslices)
                             if (remainder > 0) { 
                                 rest <- c(1:remainder, rep(remainder, nslices - remainder)) 
                                 ends <- multiples + rest
                             } else { 
                                 ends <- multiples
                             }
                             if (prod(ends) <= 0) 
                                 stop("Error happens in choosing endpoints.\n", call. = FALSE)
                             return(ends)
            }
            slicing <- function(y, nslices = 3) { 
                           n <- length(y)
                           y_asce <- sort(y)
                           ends_r <- rightends(n, nslices)
                           ends_l <- c(0, ends_r[-nslices]) + 1
                           ends <- cbind(ends_l, ends_r)
                           if (nslices == 1) { 
                               bins_y <- list(y_asce)
                           } else if (nslices == n) { 
                               bins_y <- as.list(y_asce)
                           } else {
                               bins_y <- apply(ends, 1, function(x, y) { list(y[x[1]:x[2]]) }, y_asce)
                               bins_y <- lapply(bins_y, unlist)
                           }
                           return(bins_y)
            }
            case <- match.arg(case)
            basis_method <- c("poly", "categ", "fourier", "pdisc", "pcont", "nps", "BS")
            if (!(case %in% basis_method)) stop("Incorrect choice of basis functions.\n", call. = FALSE)
            n <- length(y)
            if (nslices > n) 
                stop("'nslices' must not be larger than 'n'.\n", call. = FALSE)
            if (nslices <= 0) 
                stop("'nslices' must be positive.\n", call. = FALSE)
            if (case == "categ") { 
                if (is.factor(y)) y <- as.character.factor(y)
                bins_y <- unique(sort(y))
                r <- length(bins_y) - 1
                fy0 <- sapply(bins_y, Indicator, y)
                fy <- fy0[, 1:r]
            }
            else if (case == "fourier") { 
                fy <- array(0, c(n, 2 * degree))
                fy[, 2 * (1:degree) - 1] <- sapply(1:degree, function(dg, y) { cos(2 * pi * y * dg) }, y)
                fy[, 2 * (1:degree)] <- sapply(1:degree, function(dg, y) { sin(2 * pi * y * dg) }, y)
            }
            else if (case == "poly") { 
                if (degree == 0) stop("This case is not defined.\n", call. = FALSE)
                fy <- sapply(1:degree, function(dg, y) { y^dg }, y)
            }
            else if (case == "pdisc") { 
                if ((nslices == 0) | (nslices == 1)) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                r <- (degree + 1) * nslices - 1
                bins_y <- slicing(y, nslices)
                F0 <- sapply(bins_y, Indicator, y)
                Y0 <- matrix(y, nrow = n, ncol = nslices) - 
                          t(matrix(sapply(bins_y, function(x) { x[1] }), nrow = nslices, ncol = n))
                fy0 <- array(0, c(n, r + 1))
                if (degree == 0) { 
                    fy0 <- F0
                    fy <- fy0[, 1:r]
                } else if (degree == 1) { 
                    fy0[, 2 * (1:nslices) - 1] <- F0
                    fy0[, 2 * (1:nslices)] <- F0 * Y0
                    fy <- fy0[, -(r + 1 - degree)]
                } else if (degree == 2) { 
                    fy0[, 3 * (1:nslices) - 2] <- F0
                    fy0[, 3 * (1:nslices) - 1] <- F0 * Y0
                    fy0[, 3 * (1:nslices)] <- F0 * Y0^2
                    fy <- fy0[, -(r + 1 - degree)]
                } else if (degree == 3) { 
                    fy0[, 4 * (1:nslices) - 3] <- F0
                    fy0[, 4 * (1:nslices) - 2] <- F0 * Y0
                    fy0[, 4 * (1:nslices) - 1] <- F0 * Y0^2
                    fy0[, 4 * (1:nslices)] <- F0 * Y0^3
                    fy <- fy0[, -(r + 1 - degree)]
                }
            }
            else if (case == "pcont") { 
                if ((nslices == 0) | (nslices == 1)) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                if (degree == 0) 
                    stop("Piecewise constant continuous is not defined.\n", call. = FALSE)
                r <- nslices * degree
                fy <- array(0, c(n, r))
                bins_y <- slicing(y, nslices)
                F0 <- sapply(1:nslices, function(s, L) { rowSums(cbind(L[, s:nslices], 0)) }, 
                                 sapply(bins_y, Indicator, y, simplify = "matrix"), 
                                 simplify = "matrix")
                Y0 <- matrix(y, nrow = n, ncol = nslices) - 
                          t(matrix(sapply(bins_y, function(x) { x[1] }), nrow = nslices, ncol = n))
                Y0[, 1] <- y
                if (degree == 1) { 
                    fy <- F0 * Y0
                } else if (degree == 2) { 
                    fy[, 2 * (1:nslices) - 1] <- F0 * Y0
                    fy[, 2 * (1:nslices)] <- F0 * Y0^2
                } else if (degree == 3) { 
                    fy[, 3 * (1:nslices) - 2] <- F0 * Y0
                    fy[, 3 * (1:nslices) - 1] <- F0 * Y0^2
                    fy[, 3 * (1:nslices)] <- F0 * Y0^3
                }
            }
            if (case == "nps") { 
                if ((nslices == 0) | (nslices == 1)) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                if (degree < 3) 
                    stop("The minimum degree of power is 3.\n", call. = FALSE)
                r <- nslices + degree - 1
                fy <- array(0, c(n, r))
                bins_y <- slicing(y, nslices)
                F0 <- sapply(1:nslices, function(s, L) { rowSums(cbind(L[, s:nslices], 0)) }, 
                                 sapply(bins_y, Indicator, y, simplify = "matrix"), 
                                 simplify = "matrix")
                Y0 <- matrix(y, nrow = n, ncol = nslices) - 
                          t(matrix(sapply(bins_y, function(x) { x[1] }), nrow = nslices, ncol = n))
                fy[, 1:degree] <- sapply(1:degree, function(dg, y) { y^dg }, y)
                fy[, (1 + degree):(nslices + degree - 1)] <- F0[, -1] * Y0[, -1]^degree
            }
            else if (case == "BS") { 
                if (nslices <= 1) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                if (degree <= 0) 
                    stop("The polynomials must have positive powers.\n", call. = FALSE)
                df <- nslices + degree
                Knots <- smoothCon(s(y, bs = "bs", k = df, m = c(degree, 2)), data = data.frame(y), 
                                   knots = NULL, absorb.cons = TRUE)[[1]]$knots
                Knots[(degree + 1):(df + 1)] <- quantile(y, seq(from = 0, to = 1, 
                                                                length.out = nslices + 1))
                Knots <- data.frame(y = Knots)
                fy <- smoothCon(s(y, bs = "bs", k = df, m = c(degree, 2)), data = data.frame(y), 
                                knots = Knots, absorb.cons = TRUE)[[1]]$X
            }

            Fy0 <- Re(fy)
            Foy <- scale(Fy0, center = center, scale = scale)
            return(Foy)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
################### Generate the basis functions $F_s(y)$ in the selection equation ###################
#######################################################################################################
bF_s <- function(z, case = c("const", "poly", "categ", "fourier", "pcont", "pdisc","nps", "BS"), 
                 degree = 1, nslices = 1, center = FALSE, scale = FALSE) { 
            require(splines)
            require(mgcv)
            basis_method <- c("const", "poly", "categ", "fourier", "pdisc", "pcont", "nps", "BS")
            if (!(case %in% basis_method)) 
                stop("Incorrect choice of basis functions.\n", call. = FALSE)
            if (case != "const") {
                Fz0 <- bF_o(z, case = case, degree = degree, nslices = nslices, 
                            center = center, scale = scale)
                Fz1 <- scale(Fz0, center = center, scale = scale)
                Fz <- cbind(1, Fz1)
                return(Fz)
            } else {
                Fz <- matrix(1, nrow = length(z), ncol = 1)
                return(Fz)
            }
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
################################# Power of A Positive Definite Matrix #################################
#######################################################################################################
pdm_power <- function(A, power) { 
                 if (!is.numeric(power) | length(power) > 1 | power^2 < 0) 
                     stop("the power must be a real number\n", call. = FALSE)
                 if (!is.matrix(A) | nrow(A) != ncol(A)) 
                     stop("The object must be a square matrix\n", call. = FALSE)
                 if (abs(max(A - t(A))) > 1e-04)
                     warning("The object is not a symmetric matrix\n", call. = FALSE)
                 A <- Re(A + t(A)) / 2
                 eigens <- eigen(A)
                 eig_val <- eigens[[1]]
                 if (any(Re(eig_val) <= 0) | any(Im(eig_val) != 0)) 
                     stop("The object must be a positive definite matrix\n", call. = FALSE)
                 eig_vec <- eigens[[2]]
                 eig_val_new <- sign(eig_val) * (abs(eig_val))^power
                 A_power <- eig_vec %*% diag(eig_val_new) %*% t(eig_vec)
                 return(A_power)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
############################### Inverse Mills Ratio with its Derivative ###############################
#######################################################################################################
zeta1 <- function(x) { 
             z1 <- ifelse(x > -50, 
                          exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE)), 
                          -x / (1 - 1 / (x^2 + 2) + 1 / ((x^2 + 2) * (x^2 + 4)) - 5 / ((x^2 + 2) * 
                              (x^2 + 4) * (x^2 + 6)) + 9 / ((x^2 + 2) * (x^2 + 4) * (x^2 + 6) * 
                              (x^2 + 8)) - 129 / ((x^2 + 2) * (x^2 + 4) * (x^2 + 6) * (x^2 + 8) * 
                              (x^2 + 10)))
                          )
              return(z1)
}

zeta2 <- function(x) { 
             z2 <- -zeta1(x) * (x + zeta1(x))
             return(z2)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
############################## Joint Mode Approximation (Hsu & Wu, 2013) ##############################
#######################################################################################################
mode_HW <- function(X) { 
    if (!is.matrix(X)) stop("The object must be of matrix type.\n")
    n <- nrow(X); p <- ncol(X)
    if (n <= p) stop("The sample size must be larger than the number of random variables.\n")
    if (min(X) <= 0) { c0 <- 2 } else { c0 <- 1 }
    X0 <- X + c0 * abs(min(X))
    BoxCox <- function(A, lambda) { 
        A <- as.matrix(A)
        if (min(A) <= 0) stop("All entries of data matrix must be positive.\n")
        if (!is.vector(lambda)) stop("The parameter lambda must be of vector type.\n")
        if (ncol(A) != length(lambda)) 
            stop("The column number of data matrix must match with the length of lambda.\n")
        A_bc <- A
        for (j in 1:ncol(A)) { 
            if (abs(lambda[j]) < 1e-06) { 
                A_bc[, j] <- log(A[, j])
            } else { 
                A_bc[, j] <- (A[, j]^lambda[j] - 1) / lambda[j]
            }
        }
        return(A_bc)
    }
    f1_obj <- function(lambda, Z) { 
        if (ncol(Z) != length(lambda)) 
            stop("The column number of data matrix must match with the length of lambda.\n")
        n <- nrow(Z); p <- ncol(Z)
        Q1 <- sum(colSums(log(Z)) * lambda)
        Z_bc <- BoxCox(Z, lambda)
        mu_bc <- colMeans(Z_bc)
        Sigma_bc <- tcrossprod(t(Z_bc) - mu_bc) / n
        obj <- Q1 - n / 2 * log(det(Sigma_bc)) - n * p / 2
        return(obj)
    }
    lambda_ini <- rep(0, p)
    opt_l <- optim(lambda_ini, fn = f1_obj, Z = X0, method = "BFGS", 
                   control = list(fnscale = -1, maxit = 1e04))
    lambda_hat <- opt_l$par
    X1 <- BoxCox(X0, lambda_hat)
    mu1_hat <- colMeans(X1)
    Sigma1_hat <- tcrossprod(t(X1) - mu1_hat) / n
    require(nleqslv)
    f2_obj <- function(m, lambda, mu1, Sigma1) { 
        m1 <- as.vector(BoxCox(matrix(m, nrow = 1), lambda))
        diag(m^lambda) %*% solve(Sigma1) %*% (m1 - mu1) - (lambda - 1)
    }
    m_ini <- apply(X0, 2, function(z) { 
                        bw0 <- bw.nrd(z)
                        ds0 <- density(z, bw = bw0)
                        ds0$x[which.max(ds0$y)]
                   })
    slv_m <- nleqslv(m_ini, fn = f2_obj, method = "Broyden", global = "dbldog", 
                     lambda = lambda_hat, mu1 = mu1_hat, Sigma1 = Sigma1_hat, 
                     control = list(allowSingular = TRUE))
    m_hat <- slv_m$x
    mode_hat <- m_hat - c0 * abs(min(X))
    return(mode_hat)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

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

#######################################################################################################
############################### Initialisation for the Outcome Equation ###############################
#######################################################################################################
init_outeq <- function(X, Fy, rank, method) { 
    ### It's full-rank if rank == p, otherwise it's reduced-rank ###
    ### Usually it's set p == d for reduced-rank case if d is already known ###
    ### method = c("MM", "ML", "FOBI", "JADE", "MODE", "MaxSkew") ###
    if (!is.matrix(X)) 
        stop("The sample data of covariates must be used in a matrix form.\n", call. = FALSE)
    if (!is.matrix(Fy))
        stop("he data of 'Fy' must be used in a matrix form.\n", call. = FALSE)
    n1 <- nrow(X); n2 <- nrow(Fy)
    if (n1 != n2) {
        stop("The dimensions between two data matrices X and Fy are not conformable.\n", 
             call. = FALSE)
    } else { n <- n1 }
    p <- ncol(X); r <- ncol(Fy)
    if (p <= 1 | r < 1) 
        stop("The dimension of data matrix is too small.\n", call. = FALSE)
    if (qr(X)$rank < p)
        stop("Data matrix 'X' must have of full column rank.\n", call. = FALSE)
    if (rank <= 0 | rank > p) 
        stop("Inappropriate choice of rank for the inverse regression.\n", call. = FALSE)
    ### Initialisation based on the conditional inverse mean regression function ###
    if (rank == p) { 
    ### Full-Rank inver regression ###
    ### ref. Azzalini & Capitanio (2014) ###
        fit01 <- lm.fit(x = cbind(1, Fy), y = X, method = "qr")
        Coef01 <- as.matrix(coef(fit01))
        Beta_st <- Coef01[-1, ]
        mu0 <- as.vector(Coef01[1, ])
        rsd0 <- resid(fit01)
        Beta_st_P1 <- diag(p)
        Beta_st_P2 <- Beta_st %*% solve(Beta_st_P1)
    } else if (rank < p) { 
    ### Reduced-Rank inverse regression based on canonical variates ###
    ### ref. Izenman (2008) ###
        Cov_XX <- cov(X)
        Cov_FF <- cov(Fy)
        Cov_FX <- cov(Fy, X)
        G_nh <- pdm_power(Cov_XX, -0.5)
        W0 <- G_nh %*% t(Cov_FX) %*% solve(Cov_FF) %*% Cov_FX %*% G_nh
        V0 <- Re(eigen(W0)[[2]][, 1:rank])
        Beta_st_P1 <- solve(Cov_FF) %*% Cov_FX %*% G_nh %*% V0
        Beta_st_P2 <- solve(G_nh) %*% V0
        mu_X <- colMeans(X)
        mu_F <- colMeans(Fy)
        mu0 <- as.vector(mu_X - Beta_st_P2 %*% t(Beta_st_P1) %*% mu_F)
        rsd0 <- t(t(X) - mu0) - Fy %*% Beta_st_P1 %*% t(Beta_st_P2)
    }
    ### Initialize the parameters related to the multivariate skew-normal distribution ###
    if (method == "MM") { 
    ### Moment-based Method ###
    ### ref. Arellano-Valle & Azzalini (2008) ###
    ### ref. Azzalini & Capitanio (2014) ###
        rsd <- as.matrix(rsd0)
        mu_rsd <- apply(rsd, 2, mean)
        sigma_rsd <- apply(rsd, 2, sd)
        rsd_std <- t(t(rsd) - mu_rsd) / sigma_rsd
        gamma1 <- apply(rsd_std^3, 2, mean)
        out <- (abs(gamma1) > 0.99527)
        gamma1[out] <- sign(gamma1[out]) * 0.995
        c0 <- sign(gamma1) * (2 * abs(gamma1) / (4 - pi))^(1/3)
        mu_z <- c0 / sqrt(1 + c0^2)
        eta <- sqrt(pi/2) * mu_z
        omega <- sigma_rsd / sqrt(1 - mu_z^2)
        Sigma_st <- var(rsd)
        h_st <- omega * eta
        tmp01 <- (pi - 2) * sum((chol(solve(Sigma_st)) %*% h_st)^2) / pi
        if (tmp01 >= 1) { h_st <- h_st * 0.995 / sqrt(tmp01) } else { h_st <- h_st }
        Delta_st <- Sigma_st + 2/pi * tcrossprod(h_st)
        omega_st <- sqrt(diag(Delta_st))
        D_st <- Delta_st - tcrossprod(h_st)
        alpha_st <- as.vector(solve(D_st) %*% h_st) * omega_st / 
                        sqrt(1 + sum((chol(solve(D_st)) %*% h_st)^2))
        xi_st <- mu0 + mu_rsd - sqrt(2 / pi) * h_st
        xi_ini <- xi_st
        D_ini <- D_st
        Delta_ini <- Delta_st
        h_ini <- h_st
        alpha_ini <- alpha_st
    } else if (method == "ML") { 
    ### Maximum-Likelihood Method ###
    ### ref. Azzalini & Capitanio (1999, 2014) ###
    ### Finding moment-based starting values for MLE ###
        rsd <- as.matrix(rsd0)
        mu_rsd <- apply(rsd, 2, mean)
        sigma_rsd <- apply(rsd, 2, sd)
        rsd_std <- t(t(rsd) - mu_rsd) / sigma_rsd
        gamma1 <- apply(rsd_std^3, 2, mean)
        out <- (abs(gamma1) > 0.99527)
        gamma1[out] <- sign(gamma1[out]) * 0.995
        c0 <- sign(gamma1) * (2 * abs(gamma1) / (4 - pi))^(1/3)
        mu_z <- c0 / sqrt(1 + c0^2)
        eta <- sqrt(pi/2) * mu_z
        omega <- sigma_rsd / sqrt(1 - mu_z^2)
        Sigma_st <- var(rsd)
        h_st <- omega * eta
        tmp01 <- (pi - 2) * sum((chol(solve(Sigma_st)) %*% h_st)^2) / pi
        if (tmp01 >= 1) { h_st <- h_st * 0.995 / sqrt(tmp01) } else { h_st <- h_st }
        Delta_st <- Sigma_st + 2/pi * tcrossprod(h_st)
        omega_st <- sqrt(diag(Delta_st))
        D_st <- Delta_st - tcrossprod(h_st)
        alpha_st <- as.vector(solve(D_st) %*% h_st) * omega_st / 
                        sqrt(1 + sum((chol(solve(D_st)) %*% h_st)^2))
        beta_st <- mu_rsd - sqrt(2 / pi) * h_st
        start_ML <- list(beta = beta_st, Omega = Delta_st, alpha = alpha_st, omega = omega_st)
        ### SELM Fitting from R Package: sn ###
        require(sn)
        fit_ML <- msn.mle(y = rsd, start = start_ML, opt.method = "BFGS", 
                          control = list(maxit = 1e04))
        dp_ML <- fit_ML$dp
        xi_ini <- mu0 + as.vector(dp_ML$beta)
        Delta_ini <- dp_ML$Omega
        alpha_ini <- as.vector(dp_ML$alpha)
        h_ini <- as.vector(cov2cor(Delta_ini) %*% alpha_ini) * sqrt(diag(Delta_ini)) / 
                     sqrt(1 + fit_ML$aux$alpha.star^2)
        D_ini <- Delta_ini - tcrossprod(h_ini)
    } else if (method == "FOBI") { 
    ### FOBI: Fourth-Order Blind Identification ###
    ### It's an ICS (Invariant Coordinate Selection) based, and ###
    ### also an ICA (Independent Component Analysis) based method ###
    ### This type of method is originated from the canonical form of MSN distribution ###
    ### ref. Capitanio (2020), Azzalini & Capitanio (2014), Loperfido (2010) ###
    ### ref. Tyler, Critchley, D黰bgen & Oja (2009), Cardoso (1989), Nordhausen & Oja (2018) ###
    ### ref. Oja, Sirki?& Eriksson (2006), Norhausen, Oja & Ollila (2010) ###
    ### ref. Ilmonen, Nevalainen & Oja (2010), Ilmonen, Oja & Serfling (2012) ###
    ### ref. Miettinen, Taskinen, Nordhausen & Oja (2015), Arevalillo & Navarro (2021) ###
        rsd <- as.matrix(rsd0)
    ### Whitening residuals ###
        mu_rsd <- apply(rsd, 2, mean)
    ### The sample covariance matrix is also the first scatter matrix ###
        Sigma_rsd <- cov(rsd) * (1 - 1/n)
        G0_rsd <- pdm_power(Sigma_rsd, -0.5)
        Z_rsd <- t(t(rsd) - mu_rsd) %*% G0_rsd
    ### The kurtosis matrix is a one-step M-functional and also the second scatter matrix ###
        K_z <- matrix(rowMeans(apply(Z_rsd, 1, function(z) { 
                                         as.numeric(crossprod(z)) * tcrossprod(z)
                                     })), ncol = p) / (p + 2)
        V_z <- eigen(K_z)[[2]]
        IC_z <- Z_rsd %*% V_z
    ### Choose the vector maximizing the skewness of the indpendent components respectively ###
    ### Select the two largest absolute skewness values ###
    ### If their difference is very small, choose the one with larger kurtosis ###
        skew_z <- colMeans(IC_z^3)
        kurt_z <- colMeans(IC_z^4) - 3
        sgn_sk <- sign(skew_z)
        idx_biv <- order(abs(skew_z), decreasing = TRUE)[1:2]
        if (abs(diff(abs(skew_z[idx_biv]))) < 5e-03 ) { 
            idx <- idx_biv[which.max(kurt_z[idx_biv])]
        } else { idx <- idx_biv[1] }
    ### Obtain the desired direction related to the skewness vector ###
        dir1 <- as.vector((G0_rsd %*% V_z)[, idx]) * sgn_sk[idx]
        dir1_st <- dir1 / sqrt(sum(dir1^2))
    ### Obtain the initial estimates of MSN based on Loperfido (2010) ###
    ### Compute the multivariate skewness based on Mardia (1970) ###
        b1_ast <- min(0.99527^2 * 0.995, mean((tcrossprod(Z_rsd))^3))
        quad01 <- (b1_ast / (2 * (4 - pi)^2))^(1/3) * 2
        alpha_ast_st <- quad01 * pi / (2 - quad01 * (pi - 2))
        Sigma_MM <- var(rsd)
        quad02 <- sum((chol(Sigma_MM) %*% dir1_st)^2)
        size_dir1_st <- alpha_ast_st * (pi + (pi - 2) * alpha_ast_st) / 
                            (pi * quad02 * (1 + alpha_ast_st))
        psi_st <- dir1_st * sqrt(size_dir1_st)
        Delta_inv <- solve(Sigma_MM) - 2 / (pi + (pi - 2) * alpha_ast_st) * tcrossprod(psi_st)
        Delta_ini <- solve(Delta_inv)
        omega_st <- sqrt(diag(Delta_ini))
        alpha_ini <- psi_st * omega_st
        h_ini <- as.vector(Delta_ini %*% psi_st) / sqrt(1 + alpha_ast_st)
        D_ini <- Delta_ini - tcrossprod(h_ini)
        xi_ini <- mu0 + mu_rsd - sqrt(2 / pi) * h_ini
    } else if (method == "JADE") { 
    ### JADE: Joint Approximate Diagonalisation of Eigen-matrices ###
    ### JADE is also an IC (invariant coordinate) functional if ###
    ### cyclic Jacobi iteration procedure is employed in the algorithm ###
    ### ref. Tyler et al. (2009), Cardoso & Souloumiac (1993), Nordhausen & Oja (2018) ###
    ### ref. Oja et al. (2006), Norhausen et al. (2010), Ilmonen et al. (2010, 2012) ###
    ### ref. Miettinen et al. (2015), Miettinen, Nordhausen & Taskinen (2017) ###
        rsd <- as.matrix(rsd0)
        mu_rsd <- apply(rsd, 2, mean)
        dir1 <- JADE_Skew(rsd)
        dir1_st <- dir1 / sqrt(sum(dir1^2))
        Sigma_rsd <- cov(rsd)
        G0_rsd <- pdm_power(Sigma_rsd, -0.5)
        Z_rsd <- t(t(rsd) - mu_rsd) %*% G0_rsd
    ### Obtain the initial estimates of MSN based on Loperfido (2010) ###
    ### Compute the multivariate skewness based on Mardia (1970) ###
        b1_ast <- min(0.99527^2 * 0.995, mean((tcrossprod(Z_rsd))^3))
        quad01 <- (b1_ast / (2 * (4 - pi)^2))^(1/3) * 2
        alpha_ast_st <- quad01 * pi / (2 - quad01 * (pi - 2))
        Sigma_MM <- var(rsd)
        quad02 <- sum((chol(Sigma_MM) %*% dir1_st)^2)
        size_dir1_st <- alpha_ast_st * (pi + (pi - 2) * alpha_ast_st) / 
                            (pi * quad02 * (1 + alpha_ast_st))
        psi_st <- dir1_st * sqrt(size_dir1_st)
        Delta_inv <- solve(Sigma_MM) - 2 / (pi + (pi - 2) * alpha_ast_st) * tcrossprod(psi_st)
        Delta_ini <- solve(Delta_inv)
        omega_st <- sqrt(diag(Delta_ini))
        alpha_ini <- psi_st * omega_st
        h_ini <- as.vector(Delta_ini %*% psi_st) / sqrt(1 + alpha_ast_st)
        D_ini <- Delta_ini - tcrossprod(h_ini)
        xi_ini <- mu0 + mu_rsd - sqrt(2 / pi) * h_ini
    } else if (method == "MODE") { 
    ### Based on the mode equation from Capitanio (2020) ###
    ### Use the method proposed by Hsu & Wu (2013) to solve the joint mode ###
    ### ref. Capitanio (2020), Azzalini & Capitanio (2014), Loperfido (2010) ###
    ### ref. Hsu & Wu (2013), Bickel (2003), Chac髇 (2020) ###
    ### ref. Kirschstein, Liebscher, Porzio, Ragozini (2016), Sager (1978, 1979) ###
    ### ref. Klemel?2005, 2009), Devroye (1977) ###
        rsd <- as.matrix(rsd0)
        mode_rsd <- mode_HW(rsd)
        mu_rsd <- apply(rsd, 2, mean)
        dir1 <- mode_rsd - mu_rsd
        dir1_h <- dir1 / sqrt(sum(dir1^2))
        Sigma_rsd <- cov(rsd)
        G0_rsd <- pdm_power(Sigma_rsd, -0.5)
        Z_rsd <- t(t(rsd) - mu_rsd) %*% G0_rsd
        skew_MO <- mean(scale(Z_rsd %*% dir1_h)^3)
        if (skew_MO < 0) { 
            dir1_h <- -dir1_h
            skew_MO <- -skew_MO
        }
    ### Following the principle of Loperfido (2010) ###
        b1_ast <- min(0.99527^2 * 0.995, mean((tcrossprod(Z_rsd))^3))
        quad01 <- (b1_ast / (2 * (4 - pi)^2))^(1/3) * 2
        alpha_ast_st <- quad01 * pi / (2 - quad01 * (pi - 2))
        Sigma_MM <- var(rsd)
        quad03 <- sum((chol(solve(Sigma_MM)) %*% dir1_h)^2)
        size_dir1_h <- pi * (alpha_ast_st) / (pi + (pi - 2) * alpha_ast_st) / quad03
        h_ini <- dir1_h * sqrt(size_dir1_h)
        Delta_ini <- Sigma_MM + 2/pi * tcrossprod(h_ini)
        D_ini <- Delta_ini - tcrossprod(h_ini)
        omega_st <- sqrt(diag(Delta_ini))
        alpha_ini <- as.vector(solve(D_ini) %*% h_ini) * omega_st / sqrt(1 + alpha_ast_st)
        xi_ini <- mu0 + mu_rsd - sqrt(2 / pi) * h_ini
    } else if (method == "MaxSkew") { 
    ### Based on Loperfido (2010, 2018), Malkovich & Afifi (1973) ###
    ### ref. Golub & van Loan (1989), de Lathauwer, de Moor & Vandewalle (2000a, b) ###
    ### ref. de Lathauwer, Comon, de Moor & Vandewalle (1996), Arevalillo & Navarro (2021) ###
        rsd <- as.matrix(rsd0)
        mu_rsd <- apply(rsd, 2, mean)
        dir1 <- DomDir_Skew(rsd)
        dir1_st <- dir1 / sqrt(sum(dir1^2))
        Sigma_rsd <- cov(rsd)
        G0_rsd <- pdm_power(Sigma_rsd, -0.5)
        Z_rsd <- t(t(rsd) - mu_rsd) %*% G0_rsd
    ### Obtain the initial estimates of MSN based on Loperfido (2010) ###
    ### Compute the multivariate skewness based on Mardia (1970) ###
        b1_ast <- min(0.99527^2 * 0.995, mean((tcrossprod(Z_rsd))^3))
        quad01 <- (b1_ast / (2 * (4 - pi)^2))^(1/3) * 2
        alpha_ast_st <- quad01 * pi / (2 - quad01 * (pi - 2))
        Sigma_MM <- var(rsd)
        quad02 <- sum((chol(Sigma_MM) %*% dir1_st)^2)
        size_dir1_st <- alpha_ast_st * (pi + (pi - 2) * alpha_ast_st) / 
                            (pi * quad02 * (1 + alpha_ast_st))
        psi_st <- dir1_st * sqrt(size_dir1_st)
        Delta_inv <- solve(Sigma_MM) - 2 / (pi + (pi - 2) * alpha_ast_st) * tcrossprod(psi_st)
        Delta_ini <- solve(Delta_inv)
        omega_st <- sqrt(diag(Delta_ini))
        alpha_ini <- psi_st * omega_st
        h_ini <- as.vector(Delta_ini %*% psi_st) / sqrt(1 + alpha_ast_st)
        D_ini <- Delta_ini - tcrossprod(h_ini)
        xi_ini <- mu0 + mu_rsd - sqrt(2 / pi) * h_ini
    } else { 
        stop("Inappropriate choice of method.\n", call. = FALSE)
    }
    if (rank == p) { 
        B_ini <- Beta_st_P1
        C_ini <- Beta_st_P2 %*% solve(Delta_ini) %*% B_ini
        U_ini <- solve(D_ini) %*% Delta_ini %*% B_ini
    } else if (rank < p) { 
        C_ini <- Beta_st_P1
        B_ini <- solve(Delta_ini) %*% Beta_st_P2
        U_ini <- solve(D_ini) %*% Beta_st_P2
    }
    ini_param_list <- list(xi = xi_ini, D = D_ini, B = B_ini, U = U_ini, C = C_ini, 
                           h = h_ini, Delta = Delta_ini, alpha = alpha_ini, 
                           Full_Rank = (rank == p), method = method)
    return(ini_param_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
############################## Initialisation for the Selection Equation ##############################
#######################################################################################################
init_slceq <- function(index, Fz) { 
                  if (length(unique(index)) != 2) 
                      stop("Incorrect structure of the selection/missing index.\n", 
                           call. = FALSE)
                  if (!is.null(dim(index))) 
                      stop("Incorrect dimension of the selection/missing index.\n", 
                           call. = FALSE)
                  if (!is.factor(index)) index <- as.factor(index)
                  cs_ini <- rep(0, ncol(Fz))
                  ini_param_list <- list(cs = cs_ini)
                  return(ini_param_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
############################ Integrated Skewness Initialisation Procedure #############################
#######################################################################################################
init_est <- function(X_obs, Foy, Fsy, index, rank, method) { 
    Fy <- Foy; Fz <- Fsy
    ### The modell is full-rank if rank == p, otherwise it's reduced-rank ###
    ### Usually it's set p == d for reduced-rank case if d is already known ###
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
    if (rank <= 0 | rank > p)
        stop("Inappropriate choice of rank for the inverse regression.\n", call. = FALSE)
    ### method = c("MM", "ML", "FOBI", "JADE", "MODE", "MaxSkew", "Heckman-2S") ###
    methods_list <- c("MM", "ML", "FOBI", "JADE", "MODE", "MaxSkew", "Heckman-2S")
    if (!(method %in% methods_list)) 
        stop("Incorrect choice of the initialisation method\n", call. = FALSE)
    if (method == "Heckman-2S") { 
        ### 1st Step - Standard Probit Regression ###
        ini01_est <- init_slceq(index, Fz)
        cs_st <- ini01_est$cs
        ### Obtain the observed inverse Mills ratio ###
        iMr <- zeta1(Fz_o %*% cs_st)
        ### 2nd Step - MLR or RRR ###
        if (rank == p) { 
            fit02 <- lm.fit(x = cbind(1, Fy_o, iMr), y = X_obs, method = "qr")
            Coef02 <- as.matrix(coef(fit02))
            Beta_st2 <- Coef02[2:(1 + r1), ]
            xi_st <- as.vector(Coef02[1, ])
            h_st <- as.vector(Coef02[2 + r1, ])
            res02 <- resid(fit02)
            Sigma_res02 <- cov(res02)
            tmp00 <- mean(zeta2(Fz_o %*% cs_st))
            Delta_st <- Sigma_res02 - tmp00 * tcrossprod(h_st)
            tmp01 <- as.numeric(t(h_st) %*% solve(Sigma_res02) %*% h_st) * (1 + tmp00)
            if (tmp01 >= 1) { 
                h_st <- h_st * (0.995 / sqrt(tmp01))
                Delta_st <- Sigma_res02 - tmp00 * tcrossprod(h_st)
            }
            omega_st <- sqrt(diag(Delta_st))
            D_st <- Delta_st - tcrossprod(h_st)
            alpha_st <- (as.vector(solve(D_st) %*% h_st) * omega_st) / 
                            sqrt(1 + as.numeric(t(h_st) %*% solve(D_st) %*% h_st))
            xi_ini <- xi_st
            Delta_ini <- Delta_st
            D_ini <- D_st
            h_ini <- h_st
            alpha_ini <- alpha_st
            B_ini <- diag(p)
            Gamma_ini <- Delta_ini %*% B_ini
            C_ini <- Beta_st2 %*% solve(Delta_ini)
            U_ini <- solve(D_ini) %*% Gamma_ini
            cs_ini <- cs_st
        } else { 
            cvr <- cbind(Fy_o, iMr)
            cov_cvr <- cov(cvr)
            cov_rsp <- cov(X_obs)
            cov_rc <- cov(X_obs, cvr)
            cov_cr <- t(cov_rc)
            H_matrix <- solve(cov_rsp)
            sqrtm <- pdm_power(H_matrix, 0.5)
            weighted_matrix <- sqrtm %*% cov_rc %*% solve(cov_cvr) %*% cov_cr %*% sqrtm
            eigens <- eigen(weighted_matrix)
            eigenvalues <- eigens[["values"]]
            V_t0 <- eigens[["vectors"]][, 1:(rank + 1)]
            V_t <- as.matrix(V_t0, ncol = rank + 1)
            A1_t <- solve(sqrtm) %*% V_t
            A2_t <- t(V_t) %*% sqrtm %*% cov_rc %*% solve(cov_cvr)
            Beta_st3 <- A1_t %*% A2_t
            mu_rsp <- colMeans(X_obs)
            mu_cvr <- colMeans(cvr)
            mu_t <- as.vector(mu_rsp - Beta_st3 %*% mu_cvr)
            res03 <- t(t(X_obs) - mu_t - Beta_st3 %*% t(cvr))
            Sigma_res03 <- cov(res03)
            xi_st <- mu_t
            h_st <- Beta_st3[, 1 + r1]
            C_st <- t(matrix(A2_t[1:rank, 1:r1], nrow = rank))
            Gamma_st <- as.matrix(A1_t[, 1:rank], nrow = p)
            tmp00 <- mean(zeta2(Fz_o %*% cs_st))
            Delta_st <- Sigma_res03 - tmp00 * tcrossprod(h_st)
            tmp01 <- as.numeric(t(h_st) %*% solve(Sigma_res03) %*% h_st) * (1 + tmp00)
            if (tmp01 >= 1) { 
                h_st <- h_st * (0.995 / sqrt(tmp01))
                Delta_st <- Sigma_res03 - tmp00 * tcrossprod(h_st)
            }
            omega_st <- sqrt(diag(Delta_st))
            D_st <- Delta_st - tcrossprod(h_st)
            alpha_st <- (as.vector(solve(D_st) %*% h_st) * omega_st) / 
                            sqrt(1 + as.numeric(t(h_st) %*% solve(D_st) %*% h_st))
            xi_ini <- xi_st
            Delta_ini <- Delta_st
            D_ini <- D_st
            h_ini <- h_st
            alpha_ini <- alpha_st
            psi_ini <- as.vector(alpha_ini) / sqrt(diag(Delta_ini))
            Gamma_ini <- Gamma_st
            B_ini <- solve(Delta_ini) %*% Gamma_st
            C_ini <- C_st
            U_ini <- solve(D_ini) %*% Gamma_ini
            cs_ini <- cs_st
        } 
    } else { 
        ### 1st Step - Selection Equation ###
        ini01_est <- init_slceq(index, Fz)
        cs_st <- ini01_est$cs
        ### 2nd Step - Outcome Equation ###
        ini02_est <- init_outeq(X_obs, Fy_o, rank = rank, method = method)
        xi_ini <- ini02_est$xi
        D_ini <- ini02_est$D
        B_ini <- ini02_est$B
        U_ini <- ini02_est$U
        C_ini <- ini02_est$C
        h_ini <- ini02_est$h
        Delta_ini <- ini02_est$Delta
        alpha_ini <- ini02_est$alpha
        psi_ini <- as.vector(alpha_ini) / sqrt(diag(Delta_ini))
        Gamma_ini <- Delta_ini %*% B_ini
        cs_ini <- cs_st
    }
    ini_param_list <- list(Gamma = Gamma_ini, xi = xi_ini, Delta = Delta_ini,
                           psi = psi_ini, Co = C_ini, cs = cs_ini, D = D_ini, 
                           h = h_ini, alpha = alpha_ini, B = B_ini, U = U_ini, 
                           Full_Rank = (rank == p), method = method)
    return(ini_param_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
######################## E-Step of ECM Algorithm for Initialisation Procedure #########################
#######################################################################################################
LSMIR_init_estep <- function(pre_param_list, X_obs, Fy, Fz, index) { 
    ### Extraction ###
    xi <- pre_param_list$xi
    D <- pre_param_list$D
    Gamma <- pre_param_list$Gamma
    C <- pre_param_list$C
    h <- pre_param_list$h
    cs <- pre_param_list$cs
    ### Preparation ###
    n <- nrow(Fy); n1 <- nrow(X_obs); n2 <- n - n1
    p <- ncol(X_obs); r1 <- ncol(Fy); r2 <- ncol(Fz)
    if (is.factor(index)) { labels <- levels(index) } else { labels <- unique(index) }
    Fy_m <- as.matrix(Fy[which(index == labels[1]), ])
    Fy_o <- as.matrix(Fy[which(index == labels[2]), ])
    Fz_m <- as.matrix(Fz[which(index == labels[1]), ])
    Fz_o <- as.matrix(Fz[which(index == labels[2]), ])
    ### Compute E(x_i | S_i <= 0), E(x_i %*% t(x_i) | S_i <= 0) ###
    muX_m <- xi + Gamma %*% t(C) %*% t(Fy_m) - h %*% t(zeta1(-Fz_m %*% cs))
    Ff_m <- cbind(Fy_m, Fz_m)
    VX_m <- apply(Ff_m, 1, function(f) { 
                      fi <- f[1:r1]; zi <- f[-(1:r1)]
                      T1 <- (1 + sum(zi * cs) * zeta1(-sum(zi * cs))) * (h %*% t(h))
                      T2 <- tcrossprod(xi + Gamma %*% t(C) %*% fi)
                      T31 <- (xi + Gamma %*% t(C) %*% fi) %*% t(h)
                      T3 <- -zeta1(-sum(zi * cs)) * (T31 + t(T31))
                      D + T1 + T2 + T3
                  })
    ### Compute E(S_i | S_i <= 0), E(S_i^2 | S_i <= 0) ###
    muS_m <- as.vector(Fz_m %*% cs - zeta1(-Fz_m %*% cs))
    vS_m <- as.vector(1 + (Fz_m %*% cs)^2 - (Fz_m %*% cs) * zeta1(-Fz_m %*% cs))
    ### Compute E(S_i | S_i > 0, x_i), E(S_i^2 | S_i > 0, xi) ###
    te01 <- sqrt(1 + as.numeric(t(h) %*% solve(D) %*% h))
    Te02 <- t(t(X_obs) - xi - Gamma %*% t(C) %*% t(Fy_o))
    a_o <- te01 * as.vector(Fz_o %*% cs + Te02 %*% solve(D) %*% h / te01^2)
    muS_o <- (a_o + zeta1(a_o)) / te01
    vS_o <- (1 + a_o^2 + a_o * zeta1(a_o)) / te01^2
    ### Compute E(x_i * S_i | S_i <= 0) ###
    kappaXS_m <- t(as.vector(Fz_m %*% cs - zeta1(-Fz_m %*% cs)) * 
                       t(xi + Gamma %*% t(C) %*% t(Fy_m))) + h
    ### Return Results ###
    Es_list <- list(muX_m = muX_m, VX_m = VX_m, muS_m = muS_m, vS_m = vS_m, 
                    kappaXS_m = kappaXS_m, muS_o = muS_o, vS_o = vS_o)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
######################## CM-Step of ECM Algorithm for Initialisation Procedure ########################
#######################################################################################################
LSMIR_init_cmstep <- function(Es_list, pre_param_list, X_obs, Fy, Fz, index, d) { 
    ### Extraction ###
    xi0 <- pre_param_list$xi
    D0 <- pre_param_list$D
    Gamma0 <- pre_param_list$Gamma
    C0 <- pre_param_list$C
    h0 <- pre_param_list$h
    cs0 <- pre_param_list$cs
    muX_m <- Es_list$muX_m
    VX_m <- Es_list$VX_m
    muS_m <- Es_list$muS_m
    vS_m <- Es_list$vS_m
    kappaXS_m <- Es_list$kappaXS_m
    muS_o <- Es_list$muS_o
    vS_o <- Es_list$vS_o
    ### Preparation ###
    n <- nrow(Fy); n1 <- nrow(X_obs); n2 <- n - n1
    p <- ncol(X_obs); r1 <- ncol(Fy); r2 <- ncol(Fz)
    if (is.factor(index)) { labels <- levels(index) } else { labels <- unique(index) }
    Fy_m <- as.matrix(Fy[which(index == labels[1]), ])
    Fy_o <- as.matrix(Fy[which(index == labels[2]), ])
    Fz_m <- as.matrix(Fz[which(index == labels[1]), ])
    Fz_o <- as.matrix(Fz[which(index == labels[2]), ])
    Ff_m <- cbind(Fy_m, Fz_m); Ff_o <- cbind(Fy_o, Fz_o)
    fy_bar <- colMeans(Fy); fz_bar <- colMeans(Fz)
    ### Update xi ###
    MU01 <- t(X_obs) - Gamma0 %*% t(C0) %*% t(Fy_o) - h0 %*% t(muS_o - Fz_o %*% cs0)
    MU02 <- muX_m - Gamma0 %*% t(C0) %*% t(Fy_m) - h0 %*% t(muS_m - Fz_m %*% cs0)
    xi <- rowMeans(cbind(MU01, MU02))
    ### Update C and Gamma ###
    MG01 <- (t(X_obs) - h0 %*% t(muS_o - Fz_o %*% cs0)) %*% Fy_o
    MG02 <- (muX_m - h0 %*% t(muS_m - Fz_m %*% cs0)) %*% Fy_m
    MG03 <- -n * (xi %*% t(fy_bar))
    C <- solve(t(Fy) %*% Fy) %*% t(MG01 + MG02 + MG03) %*% 
             solve(D0) %*% Gamma0 %*% solve(t(Gamma0) %*% solve(D0) %*% Gamma0)
    Gamma <- (MG01 + MG02 + MG03) %*% C %*% solve(t(C) %*% t(Fy) %*% Fy %*% C)
    ### Update cs ###
    Mcs01 <- (t(X_obs) - Gamma %*% t(C) %*% t(Fy_o)) %*% Fz_o
    Mcs02 <- (muX_m - Gamma %*% t(C) %*% t(Fy_m)) %*% Fz_m
    Mcs03 <- -n * (xi %*% t(fz_bar))
    Mcs0 <- (Mcs01 + Mcs02 + Mcs03) / (1 + as.numeric(t(h0) %*% solve(D0) %*% h0))
    cs <- solve(t(Fz) %*% Fz) %*% 
              (t(Fz_o) %*% muS_o + t(Fz_m) %*% muS_m - t(Mcs0) %*% solve(D0) %*% h0)
    cs <- as.vector(cs)
    ### Update h and D ###
    MD11 <- tcrossprod(t(X_obs) - xi - Gamma %*% t(C) %*% t(Fy_o))
    MD12 <- matrix(rowSums(VX_m), nrow = p, ncol = p)
    MD131 <- -muX_m %*% t(xi + Gamma %*% t(C) %*% t(Fy_m))
    MD13 <- MD131 + t(MD131)
    MD14 <- tcrossprod(xi + Gamma %*% t(C) %*% t(Fy_m))
    MD01 <- MD11 + MD12 + MD13 + MD14
    mD21 <- sum(vS_o) + sum(vS_m)
    mD22 <- -2 * as.numeric(t(cs) %*% (t(Fz_o) %*% muS_o + t(Fz_m) %*% muS_m))
    mD23 <- as.numeric(t(cs) %*% t(Fz) %*% Fz %*% cs)
    mD02 <- mD21 + mD22 + mD23
    vh31 <- colSums(as.vector(muS_o - Fz_o %*% cs) * t(t(X_obs) - xi - Gamma %*% t(Fy_o %*% C)))
    vh32 <- rowSums(kappaXS_m)
    vh33 <- -colSums(as.vector(muS_m) * t(xi + Gamma %*% t(Fy_m %*% C)))
    vh34 <- -colSums(as.vector(Fz_m %*% cs) * t(muX_m - xi - Gamma %*% t(Fy_m %*% C)))
    vh03 <- vh31 + vh32 + vh33 + vh34
    h <- vh03 / mD02
    MD02 <- mD02 * (h %*% t(h))
    MD03 <- -vh03 %*% t(h) - h %*% t(vh03)
    D <- (MD01 + MD02 + MD03) / n
    ### Update Delta ###
    Delta <- D + h %*% t(h)
    ### Update the actual constrained log-Likelihood values ###
    qtl_out <- X_obs - t(xi + Gamma %*% t(Fy_o %*% C))
    qtl_slc_obs <- as.vector(Fz_o %*% cs + qtl_out %*% solve(Delta) %*% h) / 
                       as.numeric(sqrt(1 - t(h) %*% solve(Delta) %*% h))
    qtl_slc_mis <- as.vector(-Fz_m %*% cs)
    lL_c <- sum(dmvnorm(qtl_out, mean = rep(0, p), sigma = Delta, log = TRUE)) + 
                sum(pnorm(qtl_slc_obs, mean = 0, sd = 1, log.p = TRUE)) + 
                sum(pnorm(qtl_slc_mis, mean = 0, sd = 1, log.p = TRUE))
    ### Return Results ###
    param_list <- list(xi = xi, D = D, Gamma = Gamma, C = C, h = h, cs = cs, 
                       log_Lik = lL_c)
    return(param_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

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

#######################################################################################################
########################### E-Step of ECM-PX Algorithm for the LSM-IR Model ###########################
#######################################################################################################
LSMIR_main_estep <- function(pre_param_list, X_obs, Foy, Fsy, index) { 
    ### Extraction ###
    xi <- pre_param_list$xi
    Delta <- pre_param_list$Delta
    Gamma <- pre_param_list$Gamma
    Co <- pre_param_list$Co
    psi <- pre_param_list$psi
    cs <- pre_param_list$cs
    sigma_s <- pre_param_list$sigma_s
    ### Preparation ###
    Fy <- Foy; Fz <- Fsy
    n <- nrow(Fy); n1 <- nrow(X_obs); n2 <- n - n1
    p <- ncol(X_obs); r1 <- ncol(Fy); r2 <- ncol(Fz)
    if (is.factor(index)) { labels <- levels(index) } else { labels <- unique(index) }
    Fy_m <- as.matrix(Fy[which(index == labels[1]), ])
    Fy_o <- as.matrix(Fy[which(index == labels[2]), ])
    Fz_m <- as.matrix(Fz[which(index == labels[1]), ])
    Fz_o <- as.matrix(Fz[which(index == labels[2]), ])
    ### Intermediate parameters ###
    te01 <- sqrt(sigma_s^2 + as.numeric(t(psi) %*% Delta %*% psi))
    ###### (p \times n2) ######
    Mu_y_m <- xi + Gamma %*% t(Co) %*% t(Fy_m)
    Mu_y_o <- xi + Gamma %*% t(Co) %*% t(Fy_o)
    ###########################
    te02_m <- as.numeric(Fz_m %*% cs + t(Mu_y_m) %*% psi)
    te02_o <- as.numeric(Fz_o %*% cs + t(Mu_y_o) %*% psi)
    te03_o <- as.numeric(Fz_o %*% cs + X_obs %*% psi) / sigma_s
    tau_z_m <- te02_m / te01
    tau_z_o <- te02_o / te01
    h <- as.numeric(Delta %*% psi) / te01
    ### Compute E(x_i | S_i <= 0), E(x_i %*% t(x_i) | S_i <= 0) ###
    ### (p \times n2) ###
    muX_m <- Mu_y_m - h %*% t(zeta1(-tau_z_m))
    Ff_m <- rbind(Mu_y_m, tau_z_m)
    ### p^2 \times n2 ###
    VX_m <- apply(Ff_m, 2, function(f, Delta, h) { 
                      mui <- f[1:p]; taui <- f[1 + p]
                      T1 <- taui * zeta1(-taui) * (h %*% t(h))
                      T2 <- mui %*% t(mui)
                      T3 <- -zeta1(-taui) * (mui %*% t(h) + h %*% t(mui))
                      Delta + T1 + T2 + T3
                  }, Delta, h)
    ### Compute E(S_i | S_i <= 0), E(S_i^2 | S_i <= 0) ###
    muS_m <- te01 * (tau_z_m - zeta1(-tau_z_m))
    vS_m <- te01^2 * (1 + tau_z_m^2 - tau_z_m * zeta1(-tau_z_m))
    ### Compute E(S_i | S_i > 0, x_i), E(S_i^2 | S_i > 0, xi) ###
    muS_o <- (te03_o + zeta1(te03_o)) * sigma_s
    vS_o <- (1 + te03_o^2 + te03_o * zeta1(te03_o)) * sigma_s^2
    ### Compute E(x_i * S_i | S_i <= 0) ###
    ### p \times n2 ###
    kappaXS_m <- (t((tau_z_m - zeta1(-tau_z_m)) * t(Mu_y_m)) + h) * te01
    ### Return Results ###
    Es_list <- list(muX_m = muX_m, VX_m = VX_m, muS_m = muS_m, vS_m = vS_m, 
                    kappaXS_m = kappaXS_m, muS_o = muS_o, vS_o = vS_o)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
######################## CM(PX)-Step of ECM-PX Algorithm for the LSM-IR Model #########################
#######################################################################################################
LSMIR_main_cmstep <- function(Es_list, pre_param_list, X_obs, Foy, Fsy, index, d, psiA = NULL, 
                              csA = NULL, sigmasA = NULL) { 
 ### psiA, csA, sigmasA represent user-defined values for the parameters in missingness mechanism ###
    ### Extraction ###
    xi0 <- pre_param_list$xi
    Delta0 <- pre_param_list$Delta
    Gamma0 <- pre_param_list$Gamma
    Co0 <- pre_param_list$Co
    psi0 <- pre_param_list$psi
    cs0 <- pre_param_list$cs
    muX_m <- Es_list$muX_m
    VX_m <- Es_list$VX_m
    muS_m <- Es_list$muS_m
    vS_m <- Es_list$vS_m
    kappaXS_m <- Es_list$kappaXS_m
    muS_o <- Es_list$muS_o
    vS_o <- Es_list$vS_o
    ### Preparation ###
    Fy <- Foy; Fz <- Fsy
    n <- nrow(Fy); n1 <- nrow(X_obs); n2 <- n - n1
    p <- ncol(X_obs); r1 <- ncol(Fy); r2 <- ncol(Fz)
    if (is.factor(index)) { labels <- levels(index) } else { labels <- unique(index) }
    pos_m <- which(index == labels[1])
    pos_o <- which(index == labels[2])
    Fy_m <- as.matrix(Fy[pos_m, ]); Fy_o <- as.matrix(Fy[pos_o, ])
    Fz_m <- as.matrix(Fz[pos_m, ]); Fz_o <- as.matrix(Fz[pos_o, ])
    Ff_m <- cbind(Fy_m, Fz_m); Ff_o <- cbind(Fy_o, Fz_o)
    fy_bar <- colMeans(Fy); fz_bar <- colMeans(Fz)
    xo_bar <- colMeans(X_obs); muXm_bar <- rowMeans(muX_m)
    ### Update xi ###
    xi <- (n1 * xo_bar + n2 * muXm_bar) / n - as.vector(Gamma0 %*% t(Co0) %*% fy_bar)
    ### Update Co and Gamma ###
    MG01 <- t(X_obs) %*% Fy_o + muX_m %*% Fy_m - n * (xi %*% t(fy_bar))
    Co <- solve(t(Fy) %*% Fy) %*% t(MG01) %*% solve(Delta0) %*% Gamma0 %*% 
              solve(t(Gamma0) %*% solve(Delta0) %*% Gamma0)
    Gamma <- MG01 %*% Co %*% solve(t(Co) %*% t(Fy) %*% Fy %*% Co)
    Gamma0 <- matrix(Gamma, nrow = p, ncol = d)
    Co0 <- matrix(Co, nrow = r1, ncol = d)
    sgn_Gamma <- apply(Gamma0, 2, function(g0) { sign(g0[which.max(abs(g0))]) })
    Gamma <- t(sgn_Gamma * t(Gamma0))
    Co <- t(sgn_Gamma * t(Co0))
    ### Update Delta ###
    MD01 <- tcrossprod(t(X_obs) - xi - Gamma %*% t(Co) %*% t(Fy_o))
    MD02 <- matrix(rowSums(VX_m), nrow = p, ncol = p)
    MD031 <- -muX_m %*% t(xi + Gamma %*% t(Co) %*% t(Fy_m))
    MD03 <- MD031 + t(MD031)
    MD04 <- tcrossprod(xi + Gamma %*% t(Co) %*% t(Fy_m))
    Delta <- (MD01 + MD02 + MD03 + MD04) / n
    ### Update psi ###
    n_psi_U <- length(which(is.na(psiA)))
    if (is.null(psiA) | n_psi_U == p) { 
        vps01 <- t(X_obs) %*% muS_o + rowSums(kappaXS_m)
        vps02 <- -(t(X_obs) %*% Fz_o + muX_m %*% Fz_m) %*% cs0
        psi <- as.vector(solve(crossprod(X_obs) + MD02 + 1e-08*diag(p)) %*% (vps01 + vps02))
    } else if (!is.null(psiA) & n_psi_U > 0 & n_psi_U < p) { 
        pos_psi_K <- which(!is.na(psiA))
        pos_psi_U <- which(is.na(psiA))
        psi_K <- psiA[pos_psi_K]
        X_obs_K <- matrix(X_obs[, pos_psi_K], nrow = n1)
        X_obs_U <- matrix(X_obs[, pos_psi_U], nrow = n1)
        MD02_UU <- matrix(MD02[pos_psi_U, pos_psi_U], nrow = n_psi_U)
        MD02_UK <- matrix(MD02[pos_psi_U, pos_psi_K], nrow = n_psi_U)
        kappaXS_m_U <- matrix(kappaXS_m[pos_psi_U, ], nrow = n_psi_U)
        muX_m_U <- matrix(muX_m[pos_psi_U, ], nrow = n_psi_U)
        vpsiA1_U <- t(X_obs_U) %*% muS_o + rowSums(kappaXS_m_U)
        vpsiA2_U <- -(t(X_obs_U) %*% Fz_o + muX_m_U %*% Fz_m) %*% cs0
        vpsiA3_U <- -(t(X_obs_U) %*% X_obs_K + MD02_UK) %*% psi_K
        MpsiA4_U <- crossprod(X_obs_U) + MD02_UU
        psi_U <- as.vector(solve(MpsiA4_U) %*% (vpsiA1_U + vpsiA2_U + vpsiA3_U))
        psi <- rep(psi_K, len = p)
        psi[pos_psi_K] <- psi_K
        psi[pos_psi_U] <- psi_U
    } else if (!is.null(psiA) & n_psi_U == 0) { 
        psi <- psiA
    }
    ### Update cs ###
    if (!is.null(csA)) { 
        if (length(csA) != ncol(Fz)) 
           stop("Incorrect input of the argument of csA.\n", call. = FALSE)
        cs <- as.vector(csA)
        edf_cs <- length(cs)
    } else { 
        vcs01 <- t(Fz_o) %*% (muS_o - X_obs %*% psi)
        vcs02 <- t(Fz_m) %*% (muS_m - t(muX_m) %*% psi)
        cs <- as.vector(solve(t(Fz) %*% Fz) %*% (vcs01 + vcs02))
    }
    ### PX-Step ###
    ### Update sigma_s ###
    term_s01 <- sum(vS_o) + sum((X_obs %*% psi)^2) + sum((Fz_o %*% cs)^2) - 
                    2 * sum((X_obs %*% psi) * muS_o) - 2 * sum((Fz_o %*% cs) * muS_o) + 
                    2 * sum((X_obs %*% psi) * (Fz_o %*% cs))
    term_s02 <- sum(vS_m) + as.numeric(t(psi) %*% MD02 %*% psi) + sum((Fz_m %*% cs)^2) - 
                    2 * sum(psi * kappaXS_m) - 2 * sum((Fz_m %*% cs) * muS_m) + 
                    2 * sum((t(muX_m) %*% psi) * (Fz_m %*% cs))
    sigma_0 <- sqrt((term_s01 + term_s02) / n)
    if (!is.null(sigmasA)) { sigma_s <- sigmasA } else { 
        sigma_s <- max(sigma_0, 1e-06)
    }
    ### Scaling after updating the PX-parameter ###
    psi <- psi / sigma_s; cs <- cs / sigma_s
    sigma_px <- sigma_s; sigma_s <- 1
    ### Update the actual-constrained log-likelihood values ###
    ### Intermediate parameters ###
    te01 <- sqrt(sigma_s^2 + as.numeric(t(psi) %*% Delta %*% psi))
    ###### (p \times n2) ######
    Mu_y_m <- xi + Gamma %*% t(Co) %*% t(Fy_m)
    Mu_y_o <- xi + Gamma %*% t(Co) %*% t(Fy_o)
    ###########################
    te02_m <- as.numeric(Fz_m %*% cs + t(Mu_y_m) %*% psi)
    te02_o <- as.numeric(Fz_o %*% cs + t(Mu_y_o) %*% psi)
    te030_o <- as.numeric(Fz_o %*% cs + X_obs %*% psi)
    te03_o <- te030_o / sigma_s
    tau_z_m <- te02_m / te01
    tau_z_o <- te02_o / te01
    qtl_out <- X_obs - t(Mu_y_o)
    qtl_slc_obs <- te03_o
    qtl_slc_mis <- -tau_z_m
    lL_ac <- sum(dmvnorm(qtl_out, mean = rep(0, p), sigma = Delta, log = TRUE)) + 
                 sum(pnorm(qtl_slc_obs, mean = 0, sd = 1, log.p = TRUE)) + 
                 sum(pnorm(qtl_slc_mis, mean = 0, sd = 1, log.p = TRUE))
    ### Update the complete-data log-likelihood values ###
    Q11 <- -n*(p+1)*log(2*pi)/2 - n / 2 * determinant(Delta, logarithm = TRUE)[[1]][[1]] - 
               n / 2 * log(sigma_s^2) - sum(diag(solve(Delta) %*% crossprod(X_obs - t(Mu_y_o)))) / 2
    MQ12 <- MD02 - muX_m %*% t(Mu_y_m) - Mu_y_m %*% t(muX_m) + tcrossprod(Mu_y_m)
    Q12 <- -sum(diag(solve(Delta) %*% MQ12)) / 2
    Q13 <- -(sum(vS_o) + sum(te030_o^2) - 2 * sum(te030_o * muS_o)) / (2 * sigma_s^2)
    Q14 <- -(sum(vS_m) + as.numeric(t(psi) %*% MD02 %*% psi) + sum((Fz_m %*% cs)^2) - 
                 2 * sum(psi * kappaXS_m) - 2 * sum((Fz_m %*% cs) * muS_m) + 
                 2 * sum((t(muX_m) %*% psi) * (Fz_m %*% cs))) / (2 * sigma_s^2)
    lL_cd <- Q11 + Q12 + Q13 + Q14
    ### Return Results ###
    param_list <- list(xi = xi, Delta = Delta, Gamma = Gamma, Co = Co, psi = psi, cs = cs, 
                       sigma_s = sigma_s, log_Lik_ac = lL_ac, log_Lik_cd = lL_cd, 
                       comp_aux = te03_o, tau_o = tau_z_o, tau_m = tau_z_m, sigma0 = sigma_px)
    return(param_list)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

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
###########
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

#######################################################################################################
###################################### Wraped-up LSM-IR Function ######################################
#######################################################################################################
LSM_IR <- function(X_obs, X_mis = NULL, y, Foy, Fsy, index, d, psiA = NULL, csA = NULL, 
                   sigmasA = NULL, method = method, maxit = 1e04, eps = 1e-04, graph = FALSE) { 
    methods_list <- c("MM", "ML", "FOBI", "JADE", "MODE", "MaxSkew", "Heckman-2S")
    if (!(method %in% methods_list)) 
        stop("Incorrect choice of the initialisation method\n", call. = FALSE)
    n1 <- nrow(X_obs)
    n <- length(y)
    n2 <- n - n1
    p <- ncol(X_obs)
    Fy <- Foy; Fz <- Fsy
    ini_param_list <- init_est(X_obs, Fy, Fz, index, rank = d, method = method)
    tmpA1 <- LSMIR_init_ecm(X_obs, X_mis, Fy, Fz, index, d, ini_param_list, 
                        maxit = maxit, eps = 1e-03, graph = FALSE)
    init00 <- list()
    init00$Gamma <- tmpA1$Gamma_hat
    init00$xi <- tmpA1$xi_hat
    init00$Delta <- tmpA1$Delta_hat
    init00$Co <- tmpA1$C_hat
    psi00 <- (abs(diag(tmpA1$Delta_hat)) + 1e-16)^-0.5 * tmpA1$alpha_hat
    init00$psi <- psi00
    a00 <- sqrt(1 + abs(as.numeric(t(psi00) %*% tmpA1$Delta_hat %*% psi00)))
    muX00 <- t(tmpA1$xi_hat + tmpA1$Gamma_hat %*% t(tmpA1$C_hat) %*% t(Fy))
    cs00 <- solve(crossprod(Fz)) %*% t(Fz) %*% 
                (Fz %*% tmpA1$cs_hat * a00 - muX00 %*% psi00)
    init00$cs <- as.numeric(cs00)
    init00$Full_Rank <- ini_param_list$Full_Rank
    init00$method <- ini_param_list$method
    Xi00_y <- init00$xi + init00$Gamma %*% t(init00$Co) %*% t(Fy)
    S00_y_hat <- Fz %*% cs00 + t(Xi00_y) %*% psi00
 ### Avoid that the probability of complete-case units being technically zero. ###
 ### Scale the initial guesses of parameters psi and cs when ... ###
 ### ... the ''initial-guessed'' probabilities are too close to zero. ###
    if (min(S00_y_hat) <= -36) { 
        init00$psi <- -36 / min(S00_y_hat) * psi00
        init00$cs <- -36 / min(S00_y_hat) * cs00
    }
    n_psi_U <- length(which(is.na(psiA)))
    if (is.null(psiA) | n_psi_U == p) { 
        init00$psi <- init00$psi
    } else if (!is.null(psiA) & n_psi_U > 0 & n_psi_U < p) { 
        pos_psi_K <- which(!is.na(psiA))
        psi_K <- psiA[pos_psi_K]
        init00$psi[pos_psi_K] <- psi_K
    } else if (!is.null(psiA) & n_psi_U == 0) { 
        init00$psi <- psiA
    }
    if (is.null(csA)) { 
        init00$cs <- init00$cs
    } else if (length(csA) == ncol(Fz) & is.numeric(csA) & all(!is.na(csA))) { 
        init00$cs <- csA
    } else { 
        stop("Incorrect setup of the vector 'csA'.\n", call. = FALSE)
    }
    init00$sigma_s <- 1
    tmpA2 <- LSMIR_main_ecm(X_obs, X_mis, y, Fy, Fz, index, d, init00, psiA = psiA, csA = csA, 
                          sigmasA = sigmasA, maxit = maxit, eps = eps, graph = graph)
    return(tmpA2)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
