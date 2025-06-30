library(Matrix)
library(sn)
library(mgcv)

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
