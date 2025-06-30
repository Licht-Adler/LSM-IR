library(Matrix)
library(sn)
library(mgcv)

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
