library(mice)


icet <- function(df, impvars, compvars, powers = NULL, m = 5, maxit = 10, knn = 10L, verbose = FALSE){

  imputations <- rep(list(initialise(df)), m)
  i_miss <- lapply(subset(df, select = impvars), is.na)
  
  for(i in seq_len(m)){
    for(j in seq_len(maxit)){
      for(col in impvars){
        remain <- names(impvars)[names(impvars) != col]
        
        imputations[[i]][i_miss[[col]], col] <- tuni(imputations[[i]], col, c(remain, compvars), i_miss[[col]], powers, knn, verbose)
      }
    }
  }
  
  imputations
}


initialise <- function(df){
  
  miss <- colSums(is.na(df)) > 0
  impvars <- names(df)[miss]
  
  for(var in impvars){
    i_miss <- is.na(df[[var]])
    n_miss <- sum(i_miss)
    df[i_miss, var] <- sample(df[!i_miss, var], n_miss, replace = TRUE)
  }
  
  df
}


icet.mids <- function(df, impvars, compvars, powers = NULL, m = 5, maxit = 10, knn = 10L, verbose = FALSE){
  
  n_cols <- ncol(df)
  ignore <- names(df)[!names(df) %in% c(impvars, compvars)]
  r <- !is.na(df)
  
  predictorMatrix <- 1 - diag(1, nrow = n_cols, ncol = n_cols)
  row.names(predictorMatrix) <- colnames(predictorMatrix) <- names(df)
  predictorMatrix[ignore, ] <- predictorMatrix[, ignore] <- 0
  
  method <- setNames(rep("", ncol(df)), names(df))
  method[impvars] <- "pmm"
  
  visitSequence <- c(impvars, compvars, ignore)
  visitSequence <- setNames(nm = visitSequence)
  
  where <- mice::make.where(df, keyword = "missing")
  blocks <- setNames(as.list(names(df)), names(df))
  attr(blocks, "calltype") <- setNames(rep("type", ncol(df)), names(df))
  
  nmis <- colSums(is.na(df))
  
  df_tmp <- df
  imp <- mice:::initialize.imp(df, m, where, blocks, visitSequence, 
                        method, nmis, NULL)
  chainMean <- chainVar <- mice:::initialize.chain(blocks, maxit, m)
  
  for(k in seq_len(maxit)){
    for(i in seq_len(m)){
      
      for (j in visitSequence) {
    
          ry <- r[, j]
          wy <- where[, j]
          df_tmp[(!ry) & wy, j] <- imp[[j]][(!ry)[wy], i]
      }
      
      for(j in impvars){
        remain <- names(impvars)[names(impvars) != j]
        imp[[j]][, i] <- tuni(df_tmp, j, c(remain, compvars), where[, j], powers, knn, verbose)
      }
    }
    
    for (j in impvars) {
        chainVar[j, k, ] <- apply(imp[[j]], 2L,  var, na.rm = TRUE)
        chainMean[j, k, ] <- colMeans(as.matrix(imp[[j]]), na.rm = TRUE)
    }
  }
  
  formulas <- lapply(names(df), function(x) as.formula(paste(x, "~ 0")))
  names(formulas) <- names(df)
  
  for(j in c(impvars, compvars)){
    predictors <- c(impvars, compvars)[c(impvars, compvars) != j]
    formulas[[j]] <- mice:::extend.formula(formulas[[j]], predictors)
  }
  
  post <- setNames(vector("character", n_cols), names(df))
  blots <- setNames(vector("list", n_cols), names(df))
  
  midsobj <- list(data = df, imp = imp, m = m, where = where, 
                  blocks = blocks, call = NULL, nmis = nmis, method = method, 
                  predictorMatrix = predictorMatrix, visitSequence = visitSequence, 
                  formulas = formulas, post = post, blots = blots, seed = NULL, 
                  iteration = maxit, lastSeedValue = .Random.seed, 
                  chainMean = chainMean, chainVar = chainVar, loggedEvents = NULL, 
                  version = packageVersion("mice"), date = Sys.Date())
  oldClass(midsobj) <- "mids"
  
  return(midsobj)
}

abb <- function(df, var){
  
  cc <- which(complete.cases(df[[var]]))
  df[sample(cc, length(cc), replace = TRUE), ]
}


tuni <- function(df, impvar, compvar, i_miss, powers = NULL, knn = 10L, verbose = FALSE){
  
  if(is.null(powers)){
    powers <- seq(-2, 3, 0.2)
  }
  
  trans <- function(x, p) if(abs(p) < 10e-9) log(x) else x ^ p
  inv_trans <- function(x, p) if(abs(p) < 10e-9) exp(x) else x ^ (1 / p)
  
  lljmax <- -Inf
  
  bootstrp <- abb(df, impvar)
  
  for(p in powers){
    xp <- trans(bootstrp[[impvar]], p)
    reg <- lm(xp ~ ., data = subset(bootstrp, select = compvar))
    ll <- logLik(reg)
    
    if(abs(p) < 10e-9){
      jacobian <- sum(-log(bootstrp[[impvar]]))
    } else {
      jacobian <- sum(log(abs(p)) + (p - 1) * log(bootstrp[[impvar]]))
    }
    
    llj <- ll + jacobian
    
    if(llj > lljmax) {
      lljmax = llj
      pstar = p
		}
  }
  
  if(verbose)
    cat(sprintf("In ABB loop for %s, p* = %f", impvar, pstar), "\n")
  
  xp <- trans(df[[impvar]], pstar)
  
  xp_imp <- mice.impute.pmm(xp, !i_miss, model.matrix(~ 0 + ., data = subset(df, select = compvar)), donors = knn)
  x_imp <- inv_trans(xp_imp, pstar)
  
  x_imp
}
