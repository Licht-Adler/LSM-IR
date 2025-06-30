library(Matrix)
library(sn)
library(mgcv)

#######################################################################################################
################ Generate the basis functions $F_o(y)$ in the inverse outcome equation ################
#######################################################################################################
bF_o <- function(y, case = c("poly", "categ", "fourier", "pdisc", "pcont", "nps", "BS"), 
                 degree = 1, nslices = 3, center = FALSE, scale = FALSE) { 
            require(splines)
            require(mgcv)
            Indicator <- function(H, y) { as.numeric(y %in% H) }
            rightends <- function(n, nslices) { 
                             if ((n - round(n) != 0) | (nslices - round(nslices) != 0))
                                 stop("'n' and 'nslices' must be both integers.\n", call. = FALSE)
                             quotient <- floor(n / nslices)
                             multiples <- (1:nslices) * quotient
                             remainder <- n %% nslices
                             ends <- rep(0, nslices)
                             if (remainder > 0) { 
                                 rest <- c(1:remainder, rep(remainder, nslices - remainder)) 
                                 ends <- multiples + rest
                             } else { 
                                 ends <- multiples
                             }
                             if (prod(ends) <= 0) 
                                 stop("Error happens in choosing endpoints.\n", call. = FALSE)
                             return(ends)
            }
            slicing <- function(y, nslices = 3) { 
                           n <- length(y)
                           y_asce <- sort(y)
                           ends_r <- rightends(n, nslices)
                           ends_l <- c(0, ends_r[-nslices]) + 1
                           ends <- cbind(ends_l, ends_r)
                           if (nslices == 1) { 
                               bins_y <- list(y_asce)
                           } else if (nslices == n) { 
                               bins_y <- as.list(y_asce)
                           } else {
                               bins_y <- apply(ends, 1, function(x, y) { list(y[x[1]:x[2]]) }, y_asce)
                               bins_y <- lapply(bins_y, unlist)
                           }
                           return(bins_y)
            }
            case <- match.arg(case)
            basis_method <- c("poly", "categ", "fourier", "pdisc", "pcont", "nps", "BS")
            if (!(case %in% basis_method)) stop("Incorrect choice of basis functions.\n", call. = FALSE)
            n <- length(y)
            if (nslices > n) 
                stop("'nslices' must not be larger than 'n'.\n", call. = FALSE)
            if (nslices <= 0) 
                stop("'nslices' must be positive.\n", call. = FALSE)
            if (case == "categ") { 
                if (is.factor(y)) y <- as.character.factor(y)
                bins_y <- unique(sort(y))
                r <- length(bins_y) - 1
                fy0 <- sapply(bins_y, Indicator, y)
                fy <- fy0[, 1:r]
            }
            else if (case == "fourier") { 
                fy <- array(0, c(n, 2 * degree))
                fy[, 2 * (1:degree) - 1] <- sapply(1:degree, function(dg, y) { cos(2 * pi * y * dg) }, y)
                fy[, 2 * (1:degree)] <- sapply(1:degree, function(dg, y) { sin(2 * pi * y * dg) }, y)
            }
            else if (case == "poly") { 
                if (degree == 0) stop("This case is not defined.\n", call. = FALSE)
                fy <- sapply(1:degree, function(dg, y) { y^dg }, y)
            }
            else if (case == "pdisc") { 
                if ((nslices == 0) | (nslices == 1)) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                r <- (degree + 1) * nslices - 1
                bins_y <- slicing(y, nslices)
                F0 <- sapply(bins_y, Indicator, y)
                Y0 <- matrix(y, nrow = n, ncol = nslices) - 
                          t(matrix(sapply(bins_y, function(x) { x[1] }), nrow = nslices, ncol = n))
                fy0 <- array(0, c(n, r + 1))
                if (degree == 0) { 
                    fy0 <- F0
                    fy <- fy0[, 1:r]
                } else if (degree == 1) { 
                    fy0[, 2 * (1:nslices) - 1] <- F0
                    fy0[, 2 * (1:nslices)] <- F0 * Y0
                    fy <- fy0[, -(r + 1 - degree)]
                } else if (degree == 2) { 
                    fy0[, 3 * (1:nslices) - 2] <- F0
                    fy0[, 3 * (1:nslices) - 1] <- F0 * Y0
                    fy0[, 3 * (1:nslices)] <- F0 * Y0^2
                    fy <- fy0[, -(r + 1 - degree)]
                } else if (degree == 3) { 
                    fy0[, 4 * (1:nslices) - 3] <- F0
                    fy0[, 4 * (1:nslices) - 2] <- F0 * Y0
                    fy0[, 4 * (1:nslices) - 1] <- F0 * Y0^2
                    fy0[, 4 * (1:nslices)] <- F0 * Y0^3
                    fy <- fy0[, -(r + 1 - degree)]
                }
            }
            else if (case == "pcont") { 
                if ((nslices == 0) | (nslices == 1)) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                if (degree == 0) 
                    stop("Piecewise constant continuous is not defined.\n", call. = FALSE)
                r <- nslices * degree
                fy <- array(0, c(n, r))
                bins_y <- slicing(y, nslices)
                F0 <- sapply(1:nslices, function(s, L) { rowSums(cbind(L[, s:nslices], 0)) }, 
                                 sapply(bins_y, Indicator, y, simplify = "matrix"), 
                                 simplify = "matrix")
                Y0 <- matrix(y, nrow = n, ncol = nslices) - 
                          t(matrix(sapply(bins_y, function(x) { x[1] }), nrow = nslices, ncol = n))
                Y0[, 1] <- y
                if (degree == 1) { 
                    fy <- F0 * Y0
                } else if (degree == 2) { 
                    fy[, 2 * (1:nslices) - 1] <- F0 * Y0
                    fy[, 2 * (1:nslices)] <- F0 * Y0^2
                } else if (degree == 3) { 
                    fy[, 3 * (1:nslices) - 2] <- F0 * Y0
                    fy[, 3 * (1:nslices) - 1] <- F0 * Y0^2
                    fy[, 3 * (1:nslices)] <- F0 * Y0^3
                }
            }
            if (case == "nps") { 
                if ((nslices == 0) | (nslices == 1)) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                if (degree < 3) 
                    stop("The minimum degree of power is 3.\n", call. = FALSE)
                r <- nslices + degree - 1
                fy <- array(0, c(n, r))
                bins_y <- slicing(y, nslices)
                F0 <- sapply(1:nslices, function(s, L) { rowSums(cbind(L[, s:nslices], 0)) }, 
                                 sapply(bins_y, Indicator, y, simplify = "matrix"), 
                                 simplify = "matrix")
                Y0 <- matrix(y, nrow = n, ncol = nslices) - 
                          t(matrix(sapply(bins_y, function(x) { x[1] }), nrow = nslices, ncol = n))
                fy[, 1:degree] <- sapply(1:degree, function(dg, y) { y^dg }, y)
                fy[, (1 + degree):(nslices + degree - 1)] <- F0[, -1] * Y0[, -1]^degree
            }
            else if (case == "BS") { 
                if (nslices <= 1) 
                    stop("The minimum number of slices is 2.\n", call. = FALSE)
                if (degree <= 0) 
                    stop("The polynomials must have positive powers.\n", call. = FALSE)
                df <- nslices + degree
                Knots <- smoothCon(s(y, bs = "bs", k = df, m = c(degree, 2)), data = data.frame(y), 
                                   knots = NULL, absorb.cons = TRUE)[[1]]$knots
                Knots[(degree + 1):(df + 1)] <- quantile(y, seq(from = 0, to = 1, 
                                                                length.out = nslices + 1))
                Knots <- data.frame(y = Knots)
                fy <- smoothCon(s(y, bs = "bs", k = df, m = c(degree, 2)), data = data.frame(y), 
                                knots = Knots, absorb.cons = TRUE)[[1]]$X
            }

            Fy0 <- Re(fy)
            Foy <- scale(Fy0, center = center, scale = scale)
            return(Foy)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
################### Generate the basis functions $F_s(y)$ in the selection equation ###################
#######################################################################################################
bF_s <- function(z, case = c("const", "poly", "categ", "fourier", "pcont", "pdisc","nps", "BS"), 
                 degree = 1, nslices = 1, center = FALSE, scale = FALSE) { 
            require(splines)
            require(mgcv)
            basis_method <- c("const", "poly", "categ", "fourier", "pdisc", "pcont", "nps", "BS")
            if (!(case %in% basis_method)) 
                stop("Incorrect choice of basis functions.\n", call. = FALSE)
            if (case != "const") {
                Fz0 <- bF_o(z, case = case, degree = degree, nslices = nslices, 
                            center = center, scale = scale)
                Fz1 <- scale(Fz0, center = center, scale = scale)
                Fz <- cbind(1, Fz1)
                return(Fz)
            } else {
                Fz <- matrix(1, nrow = length(z), ncol = 1)
                return(Fz)
            }
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
