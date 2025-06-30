library(Matrix)
library(sn)
library(mgcv)

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
