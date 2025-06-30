library(Matrix)
library(sn)
library(mgcv)

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
