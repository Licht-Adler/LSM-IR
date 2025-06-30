library(Matrix)
library(sn)
library(mgcv)

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
