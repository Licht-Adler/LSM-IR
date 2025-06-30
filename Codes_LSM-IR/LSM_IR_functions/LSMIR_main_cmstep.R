library(Matrix)
library(sn)
library(mgcv)

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
