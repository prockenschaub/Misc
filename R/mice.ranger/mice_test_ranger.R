# Reproduces the simulation presented in Tables 1&1 of 
#
# Doove, L. L., S. Van Buuren, and E. Dusseldorp. 2014. “Recursive 
# Partitioning for Missing Data Imputation in the Presence of Interaction 
# Effects.” Computational Statistics & Data Analysis 72: 92–104.
# 
# for the `ranger` package. All code below is based on earlier code written
# by L.L. Doove.

library(MASS)
library(mice)
library(microbenchmark)


# Simulate the data -------------------------------------------------------

createdata <- function(n=1000, 
                       B=rep(1., 7),
                       interact = c(3, 3)){
  # Simulate the dataset with outcome y and covariates x1-x10 (only x1, x2, x3, x8, and x9 are used)
  #
  # Parameters
  # ----------
  # n : integer
  #   number of observations
  # B : numeric vector
  #   true model coefficients used to simulate y 
  # interact : integer vector
  #   vector of length 2 indicating the covariates (1-10) used to create 
  #   the interactions of different strength. 
  #     - High correlation: c(3, 3)
  #     - Medium correlation: c(1, 2)
  #     - Low correlation: c(8, 9)
  #
  # Returns
  # -------
  # data.frame
  Param <- as.vector(B)
  x <- round(mvrnorm(n=n, mu=c(rep(0,10)), Sigma=matrix(c(1,.5,.5,.5,.0,.0,.0,.0,.0,.0,
                                                          .5,1,.5,.5,.0,.0,.0,.0,.0,.0,
                                                          .5,.5,1,.5,.0,.0,.0,.0,.0,.0,
                                                          .5,.5,.5,1,.0,.0,.0,.0,.0,.0,
                                                          .0,.0,.0,.0,1,.3,.3,.3,.3,.3,
                                                          .0,.0,.0,.0,.3,1,.3,.3,.3,.3,
                                                          .0,.0,.0,.0,.3,.3,1,.3,.3,.3,
                                                          .0,.0,.0,.0,.3,.3,.3,1,.3,.3,
                                                          .0,.0,.0,.0,.3,.3,.3,.3,1,.3,
                                                          .0,.0,.0,.0,.3,.3,.3,.3,.3,1), nrow=10, byrow=T)),2)
  colnames(x) <- paste0("x", 1:10)
  
  Intercept <- as.vector(rep(0,n))
  X <- matrix(c(Intercept, x[,1], x[,2] , x[,3], x[,8], x[,9], (x[,interact[1]] * x[,interact[2]])),nrow=n)
  eps <- rnorm(n, mean=0, sd=1)
  y <- round((X%*%Param)+eps,2)
  return(cbind(as.data.frame(x), data.frame(xint=x[,interact[1]] * x[,interact[2]], y=y)))
}

makeMARmid <- function(data) {
  # Introduce missingness into the outcome y
  #
  # Parameters
  # ----------
  # data : data.frame
  #   data.frame returned by createdata()
  #
  # Returns
  # -------
  # data.frame
  logistic <- function(x) exp(x)/(1+exp(x))
  x9 <- data[,"x9"]
  x10 <- data[,"x10"]
  for (i in c("y"))
  {
    p.marmid <- 1-logistic(-.7+abs((x9*x10)-mean(x9*x10)))                      # -.7 gives 50% missing
    r.marmid <- rbinom(nrow(data),1,p.marmid)                                   # -3 gives 10% missing
    data[r.marmid==0,i] <- NA                                                   # -1.6 gives 30% missing
  }
  return(data)
}



# Custom imputation based on ranger package -------------------------------

mice.impute.rf.custom <- function (y, ry, x, wy = NULL, ntree = 10, ...) {
  mice:::install.on.demand("randomForest", ...)
  if (is.null(wy)) 
    wy <- !ry
  onetree <- function(xobs, xmis, yobs, ...) {
    fit <- randomForest::randomForest(x = xobs, y = yobs, 
                                      ntree = 1, ...)
    leafnr <- predict(object = fit, newdata = xobs, nodes = TRUE)
    leafnr <- as.vector(attr(leafnr, "nodes"))
    
    nodes <- predict(object = fit, newdata = xmis, nodes = TRUE)
    nodes <- as.vector(attr(nodes, "nodes"))
    
    donor <- lapply(nodes, function(s) yobs[leafnr == s])
    return(donor)
  }
  ntree <- max(1, ntree)
  nmis <- sum(wy)
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- y[ry]
  forest <- sapply(seq_len(ntree), FUN = function(s) onetree(xobs, 
                                                             xmis, yobs, ...))
  if (nmis == 1) 
    forest <- array(forest, dim = c(1, ntree))
  impute <- apply(forest, MARGIN = 1, FUN = function(s) sample(unlist(s), 
                                                               1))
  return(impute)
}


mice.impute.ranger.onetree <- function (y, ry, x, wy = NULL, ntree = 10, ...) {
  mice:::install.on.demand("ranger", ...)
  if (is.null(wy)) 
    wy <- !ry
  onetree <- function(xobs, xmis, yobs, ...) {
    fit <- ranger::ranger(x = xobs, y = yobs, num.trees = 1)
    
    leafnr <- predict(object = fit, data = xobs, type = "terminalNodes")
    leafnr <- as.vector(ranger::predictions(leafnr))  
    
    nodes <- predict(object = fit, data = xmis, type = "terminalNodes")
    nodes <- as.vector(ranger::predictions(nodes))  
    
    donor <- lapply(nodes, function(s) yobs[leafnr == s])
    return(donor)
  }
  ntree <- max(1, ntree)
  nmis <- sum(wy)
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- y[ry]
  forest <- sapply(seq_len(ntree), FUN = function(s) onetree(xobs, 
                                                             xmis, yobs, ...))
  if (nmis == 1) 
    forest <- array(forest, dim = c(1, ntree))
  impute <- apply(forest, MARGIN = 1, FUN = function(s) sample(unlist(s), 
                                                               1))
  return(impute)
}


mice.impute.ranger <- function (y, ry, x, wy = NULL, ntree = 10, ...) {
  mice:::install.on.demand("ranger", ...)
  if (is.null(wy)) 
    wy <- !ry
  
  ntree <- max(1, ntree)
  nmis <- sum(wy)
  xobs <- x[ry, , drop = FALSE]
  xmis <- x[wy, , drop = FALSE]
  yobs <- y[ry]
  
  # Fit forest and determine donors for imputation
  fit <- ranger::ranger(x = xobs, y = yobs, num.trees = ntree)
  
  nodes <- predict(object = fit, data = rbind(xobs, xmis), 
                   type = "terminalNodes", predict.all = TRUE)
  nodes <- ranger::predictions(nodes)
  nodes_obs <- nodes[1:nrow(xobs), ]
  nodes_mis <- nodes[(nrow(xobs)+1):nrow(nodes), ]
  
  select_donors <- function(i){
    donors <- split(yobs, nodes_obs[, i])
    donors[as.character(nodes_mis[, i])]
  }
  
  forest <- sapply(seq_len(ntree), FUN = select_donors)
  
  if (nmis == 1) 
    forest <- array(forest, dim = c(1, ntree))
  impute <- apply(forest, MARGIN = 1, FUN = function(s) sample(unlist(s), 
                                                               1))
  return(impute)
}


# Define experiment -------------------------------------------------------

run_and_eval_imputation <- function(data, method = "pmm", m = 10, maxit = 1, ntree = 100) {
  # Impute the data with the chosen parameters and extract the model estimates
  #
  # Parameters
  # ----------
  # data : data.frame
  #   the simulated data
  # method : character
  #   the method to use for imputation (pmm, rf.custom, ranger)
  # m : integer
  #   Number of imputated datasets
  # maxit : integer
  #   Number of iterations run by mice
  # ntree : integer
  #   Number of trees used for fitting the random forests
  #
  # Returns
  # -------
  # matrix (#coefs x #measurements)
  #   Results of imputation, where rows are the coefficients and columns 
  #   are measurements like mean, std. error, quantiles, fmi, and lambda
  init <- mice(data, maxit = 0, print = FALSE)
  predMat <- init$predictorMatrix
  predMat[, "xint"] <- 0
  
  imp <- mice(data, method = method, predictorMatrix = predMat,
              m = m, maxit = maxit, ntree = ntree, print=FALSE)
  fit <- with(imp, lm(y~x1+x2+x3+x8+x9+xint))
  est <- pool(fit)
  tab <- summary(est, conf.int = TRUE)
  return(
    as.matrix(
      cbind(
        tab[,c("estimate","std.error","2.5 %", "97.5 %")], 
        est$pooled[, c("fmi", "lambda")]
  )))
}


create_tables <- function(s, B = rep(1., 7)) {
  # Take the coefficient summaries returned by run_and_eval_imputation()
  # and convert them into the tables presented in the paper.
  #
  # Parameters
  # ----------
  # s : 3D array (# coefs x # measurements x # of simulations)
  #   Results from a call of run_and_eval_imputation()
  # B : numeric vector
  #   True coefficients used to simulate the data
  #
  # Returns
  # -------
  # matrix
  true <- as.matrix(B,,1)
  bias <- round(rowMeans(apply(s[,1,],2,function(s) s-true)),3)
  isinint <- apply(s[,3,],2,function(s) s < true) & apply(s[,4,],2,function(s) true < s)
  cov <- round(rowMeans(isinint),3)
  intwidth <- s[,4,] - s[,3,]
  aiw <- round(rowMeans(intwidth),3)
  fmi <- round(rowMeans(s[,5,]),3)
  lambda <- round(rowMeans(s[,6,]),3)
  prop_table <- cbind(bias,cov,aiw,fmi,lambda)
  rownames(prop_table) <- names(coef(lm(y~x1+x2+x3+x8+x9+xint, data=createdata(n=100))))
  colnames(prop_table) <- c("Bias", "Coverage", "CI width", "fmi", "Lambda")
  return(prop_table)
}


run_simulation <- function(nsim = 10, nobs = 100, 
                           B = rep(1., 7), interact = c(3, 3),
                           m = 10, maxit = 1, ntree = 100) {
  # Define and run the full experiment.
  #
  # Paramters
  # ---------
  # nsim : integer
  #   Number of simulations to run
  # nobs : integer
  #   Number of observations per simulation
  # B : numeric vector
  #   True coefficients
  # interact : numeric vector
  #   Index of the two variables that define the interaction
  # m : integer
  #   Number of imputated datasets
  # maxit : integer
  #   Number of iterations run by mice
  # ntree : integer
  #   Number of trees used for fitting the random forests
  #
  # Returns
  # -------
  # list of matrices
  set.seed(41872)
  
  create_and_ampute <- function() {
    # Simulate study data and ampute the outcome
    makeMARmid(createdata(n = nobs, B = B, interact = interact))
  }
  
  data <- replicate(nsim, create_and_ampute(), simplify = FALSE)
  
  lapply(
    list(
      pmm = sapply(data, simplify = "array",
                   run_and_eval_imputation, 
                   m = m, method="pmm", maxit = maxit), 
      rfor = sapply(data, simplify = "array",
                    run_and_eval_imputation, 
                    m = m, method="rf.custom", maxit = maxit, ntree = ntree),
      rang = sapply(data, simplify = "array",
                    run_and_eval_imputation, 
                    m = m, method="ranger", maxit = maxit, ntree = ntree)
    ),
    create_tables, B = B
  )
}



# Compare speed -----------------------------------------------------------

set.seed(1234)
B <- c(b0=0, b1=.31, b2=.31, b3=.31, b4=.31, b5=.31, b6=.28)
data <- replicate(
  5,
  makeMARmid(createdata(n = 1000, B = B, interact = c(3, 3))),
  simplify = FALSE
)

# PMM
microbenchmark(
  sapply(data, simplify = "array",
         run_and_eval_imputation,
         m = 2, method="pmm", maxit = 1),
  times = 10L
)

# randomForest
microbenchmark(
  sapply(data, simplify = "array",
         run_and_eval_imputation,
         m = 2, method="rf.custom", maxit = 1, ntree = 100),
  times = 10L
)

# ranger (implementation using the same onetree function as randomForest)
# NOTE: this implementation is much slower than randomForest, perhaps due
#       to a bottleneck in switching to C
microbenchmark(
  sapply(data, simplify = "array",
         run_and_eval_imputation,
         m = 2, method="ranger.onetree", maxit = 1, ntree = 100),
  times = 10L
)

# ranger (implementation fitting entire forest at once)
microbenchmark(
  sapply(data, simplify = "array",
         run_and_eval_imputation,
         m = 2, method="ranger", maxit = 1, ntree = 100),
  times = 10L
)



# High correlation (r = 1.0) --------------------------------------------

high <- list(
  # Small effect size (f^2 = 0.02, Table B.1)
  small = run_simulation(
    nsim = 200, nobs = 1000, interact = c(3, 3),
    B = c(b0=0, b1=.34, b2=.34, b3=.34, b4=.34, b5=.34, b6=.11), 
    m = 20, maxit = 1, ntree = 100
  ),

  # Medium effect size (f^2 = 0.15, Table B.1) == results in Table 1
  medium = run_simulation(
    nsim = 200, nobs = 1000, interact = c(3, 3),
    B = c(b0=0, b1=.31, b2=.31, b3=.31, b4=.31, b5=.31, b6=.28), 
    m = 20, maxit = 1, ntree = 100
  ),
  
  # Large effect size (f^2 = 0.35, Table B.1)
  large = run_simulation(
    nsim = 200, nobs = 1000, interact = c(3, 3),
    B = c(b0=0, b1=.28, b2=.28, b3=.28, b4=.28, b5=.28, b6=.42), 
    m = 20, maxit = 1, ntree = 100
  )
)


# Medium correlation (r = 0.5) --------------------------------------------

medium <- list(
  # Small effect size (f^2 = 0.02, Table B.1)
  small = run_simulation(
    nsim = 200, nobs = 1000, interact = c(1, 2),
    B = c(b0=0, b1=.34, b2=.34, b3=.34, b4=.34, b5=.34, b6=.13), 
    m = 20, maxit = 1, ntree = 100
  ),
  
  # Medium effect size (f^2 = 0.15, Table B.1)
  medium = run_simulation(
    nsim = 200, nobs = 1000, interact = c(1, 2),
    B = c(b0=0, b1=.31, b2=.31, b3=.31, b4=.31, b5=.31, b6=.35), 
    m = 20, maxit = 1, ntree = 100
  ),
  
  # Large effect size (f^2 = 0.35, Table B.1)
  large = run_simulation(
    nsim = 200, nobs = 1000, interact = c(1, 2),
    B = c(b0=0, b1=.28, b2=.28, b3=.28, b4=.28, b5=.28, b6=.53), 
    m = 20, maxit = 1, ntree = 100
  )
)



# Low correlation (r = 0.3) --------------------------------------------

medium <- list(
  # Small effect size (f^2 = 0.02, Table B.1)
  small = run_simulation(
    nsim = 200, nobs = 1000, interact = c(8, 9),
    B = c(b0=0, b1=.34, b2=.34, b3=.34, b4=.34, b5=.34, b6=.13), 
    m = 20, maxit = 1, ntree = 100
  ),
  
  # Medium effect size (f^2 = 0.15, Table B.1)
  medium = run_simulation(
    nsim = 200, nobs = 1000, interact = c(8, 9),
    B = c(b0=0, b1=.31, b2=.31, b3=.31, b4=.31, b5=.31, b6=.37), 
    m = 20, maxit = 1, ntree = 100
  ),
  
  # Large effect size (f^2 = 0.35, Table B.1)
  large = run_simulation(
    nsim = 200, nobs = 1000, interact = c(8, 9),
    B = c(b0=0, b1=.28, b2=.28, b3=.28, b4=.28, b5=.28, b6=.57), 
    m = 20, maxit = 1, ntree = 100
  )
)