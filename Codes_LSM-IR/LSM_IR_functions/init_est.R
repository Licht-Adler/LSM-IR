library(Matrix)
library(sn)
library(mgcv)

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
