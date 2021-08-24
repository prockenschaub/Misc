library(mice)

riwish <- function(v, S) {
  #' Draw a random number from an inverse Wishart distribution
  #' 
  #' @param v Degrees of freedom 
  #' @param S Scale matrix (p \times p)
  #' 
  #' @return A single draw (scalar) from the specified Wishart distribution
  #' 
  #' Adapted from `MCMCpack::riwish`.
  X_inv <- rWishart(1, v, solve(S * (v - ncol(S) - 1)))[, , 1]
  solve(X_inv)
}

rmatnorm <- function (M, U, V) {
  #' Draw a random number from a matrix normal distribution
  #' 
  #' @param M Mean matrix (n \times p)
  #' @param U Individual scale (n \times n)
  #' @param V Parameter scale (p \times p)
  #' 
  #' @return A single draw (n \times p-matrix) from the specified matrix normal
  #'   distribution
  #' 
  #' Adapted from `matrixNormal::rmatnorm`
  n <- nrow(M)
  p <- ncol(M)
  
  if (ncol(U) != nrow(U) || ncol(U) != n)
    stop("U must be a symmetric positive-definite matrix of size n")
  if (ncol(V) != nrow(V) || ncol(V) != p)
    stop("V must be a symmetric positive-definite matrix of size p")
  
  Sigma <- kronecker(V, U)
  vec.X <- MASS::mvrnorm(1, as.vector(M), Sigma)
  matrix(vec.X, nrow = n, ncol = p)
}

.mv.norm.draw <- function (y, ry, x, ...) {
  #' This function draws random values of beta and sigma under the multivariate 
  #' Bayesian linear regression model. The univariate case was described in 
  #' Rubin (1987, p. 167) and implemented in \code{mice::.norm.draw}. This 
  #' implementation is a generalisation of this approach to the multivariate
  #' case following Box and Tiao (1973, p. 460).
  #' 
  #' @param y Incomplete data matrix (n \times d) where d is the number of 
  #' outcomes imputed together in the block
  #' @param ry Vector of missing data pattern (FALSE=missing, TRUE=observed)
  #' @param x Matrix (n \times p) of complete covariates.
  #' @param ... Other named arguments
  #' @return A \code{list} containing components \code{coef} (least squares 
  #' estimate), \code{var} (estimated outcome variance matrix),  \code{beta} 
  #' (drawn regression weights) and \code{sigma} (drawn value of the
  #' residual standard deviation).
  #' 
  #' @note The algorithm currently assumes that a row in \code{y} is either 
  #' missing completely or observed completely.
  #' 
  #' @seealso \code{mice::.norm.draw}
  #' @references
  #' Rubin, D.B. (1987). \emph{Multiple imputation for nonresponse in surveys}. 
  #' New York: Wiley.
  #' Box, G. E. P. and Tiao, G. C. (1973). Bayesian Inference in Statistical 
  #' Analysis. Wiley.
  x_ry = x[ry, , drop = FALSE]
  y_ry = y[ry, ]
  
  # Fit the multivariate model
  fit <- lm.fit(x = x_ry, y = y_ry, ...)
  coef <- fit$coefficients
  var <- crossprod(fit$residuals) / fit$df.residual
  
  # Draw from the posterior of the model parameters
  sigma.star <- riwish(fit$df.residual, var)
  beta.star <- rmatnorm(coef, solve(t(x_ry) %*% x_ry), sigma.star)
  
  parm <- list(coef, var, beta.star, sigma.star, "mv")
  names(parm) <- c("coef", "var", "beta", "sigma", "estimation")
  parm
}

matchrann <- function(D, T, k = 5L) {
  #' Find the index of matched donor units.
  #' 
  #' k-nearest matches using the Approximate Near Neighbor (ANN) C++ 
  #' library.
  #' 
  #' @param D Numeric matrix with values from donor cases.
  #' @param T Numeric matrix with values from target cases
  #' @param k Integer, number of unique donors from which a random draw is made. 
  #' For k = 1 the function returns the index in d corresponding to the closest 
  #' unit. For multiple imputation, the advice is to set values in the range of 
  #' k = 5 to k = 10.
  #' @return An integer vector with length(t) elements. Each element is an index 
  #' in the matrix D.
  #' 
  #' @note Does not yet explicitly deal with ties in the predictions.
  mice:::install.on.demand("RANN")
  n <- nrow(T)
  idx <- RANN::nn2(D, T, k)$nn.idx
  choice <- sample(seq_len(k), size = n, replace = TRUE)
  as.vector(t(idx))[k * (seq_len(n) - 1) + choice]
}

mice.impute.mvpmm <- function(data, formula, type, format = "list", donors = 5L, 
                              matchtype = 1L, scale = TRUE, ...) {
  #' Multivariate joint imputation by predictive mean matching, as suggested in 
  #' van Buuren, S. (2018, chapter 4.7).
  #' 
  #' @param data A data frame containing incomplete and auxiliary variables.
  #' @param formula A formula specifying the role of each variable
  #' in the imputation model. 
  #' @param type Not currently used.
  #' @param format A character vector specifying the type of object that should
  #' be returned. The default is \code{format = "list"}. No other formats are
  #' currently supported.
  #' @param donors The size of the donor pool among which a draw is made.
  #' The default is \code{donors = 5L}. Setting \code{donors = 1L} always 
  #' selects the closest match, but is not recommended. Values between 3L and 
  #' 10L provide the best results in most cases (Morris et al, 2015).
  #' @param matchtype Type of matching distance. The default choice
  #' (\code{matchtype = 1L}) calculates the distance between
  #' the \emph{predicted} value of \code{yobs} and
  #' the \emph{drawn} values of \code{ymis} (called type-1 matching).
  #' Other choices are \code{matchtype = 0L}
  #' (distance between predicted values) and \code{matchtype = 2L}
  #' (distance between drawn values).
  #' @param scale Flag indicating if the predicted values for each outcome 
  #' should be scaled by the outcome variance prior to matching.
  #' @return A list of imputations for all incomplete variables in the model, 
  #' that can be stored in the the imp component of the mids object.
  #' 
  #' @note The algorithm currently assumes that a row in \code{y} is either 
  #' missing completely or observed completely.
  #' 
  #' @references 
  #' Van Buuren, S. (2018).
  #' \href{https://stefvanbuuren.name/fimd/sec-FCS.html#sec:MICE}
  #' {\emph{Flexible Imputation of Missing Data. Second Edition.}}
  #' Chapman & Hall/CRC. Boca Raton, FL.
  
  if (is.null(formula)) {
    stop(
      "mice.impute.mvpmm currently requires the use of formulas to define ",
      "the imputation model"
    )
  }
  
  mf <- model.frame(formula, data)  # TODO: find a more elegant solution to get y
  y <- as.matrix(data[, colnames(model.response(mf, "numeric"))])
  x <- mice:::obtain.design(data, formula)
  
  ry <- !is.na(y[, 1]) # TODO: assumes that all elements are jointly present or absent
  
  parm <- .mv.norm.draw(y, ry, x, ...)
  
  if (scale) {
    var_inv <- t(chol(solve(parm$var)))
    sigma_inv <- t(chol(solve(parm$sigma)))
  } else {
    var_inv <- sigma_inv <- diag(ncol(parm$var))
  }
  
  if (matchtype == 0L) {
    yhatobs <- x[ry, , drop = FALSE] %*% parm$coef %*% var_inv
    yhatmis <- x[!ry, , drop = FALSE] %*% parm$coef %*% var_inv
  }
  if (matchtype == 1L) {
    yhatobs <- x[ry, , drop = FALSE] %*% parm$coef %*% var_inv
    yhatmis <- x[!ry, , drop = FALSE] %*% parm$beta %*% sigma_inv
  }
  if (matchtype == 2L) {
    yhatobs <- x[ry, , drop = FALSE] %*% parm$beta %*% sigma_inv
    yhatmis <- x[!ry, , drop = FALSE] %*% parm$beta %*% sigma_inv
  }
  idx <- matchrann(yhatobs, yhatmis, donors)
  
  return(as.list(as.data.frame(y[ry, ][idx, ])))
}






set.seed(42)
n <- 1000
d <- 6
x <- cbind(1, MASS::mvrnorm(n, rep(0, d), 0.2 + 0.8 * diag(d)))
y <- as.vector(x %*% c(1, -0.5, 2, -1, 0.25, -0.75, 1.5) + rnorm(n))

data <- as.data.frame(x[, -1])
names(data) <- paste0("x", seq_len(d))
data$y <- y

data_mis <- data
data_mis[sample(c(FALSE, TRUE), n, TRUE, c(0.8, 0.2)), c("x1", "x2")] <- NA
data_mis[sample(c(FALSE, TRUE), n, TRUE, c(0.9, 0.1)), c("x4")] <- NA
 


blocks <- list(
  x1_x2 = c("x1", "x2"),
  x4 = "x4"
)
attr(blocks, "calltype") <- c(x1_x2 = "formula", x4 = "formula")

method <- c(x1_x2 = c("mvpmm"), x4 = "pmm")
formulas <- list(
  x1_x2 = cbind(x1, x2) ~ x3 + x4 + x5 + x6 + y,
  x4 = x4 ~ x1 + x2 + x3 + x5 + x6 + y
)

imp <- mice(data_mis, method = method, blocks = blocks, formulas = formulas)
densityplot(imp)

