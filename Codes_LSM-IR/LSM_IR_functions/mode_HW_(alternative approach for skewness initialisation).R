library(Matrix)
library(sn)
library(mgcv)

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
