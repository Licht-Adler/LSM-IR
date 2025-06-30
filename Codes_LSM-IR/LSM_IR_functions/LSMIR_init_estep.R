library(Matrix)
library(sn)
library(mgcv)

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
