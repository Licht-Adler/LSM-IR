library(Matrix)
library(sn)
library(mgcv)

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
