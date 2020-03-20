mfpmi <- function (formula, mids, family = gaussian, alpha = 0.05, select = 1, maxits = 20, verbose = FALSE) {
  
  Terms <- terms(formula, specials = "fp")
  m <- model.frame(formula, complete(imp2)[0, ])
  
  nx <- ncol(m) - 1
  nobs <- nrow(mids$data)
  df.list <- rep(1, nx)
  alpha.list <- rep(alpha, nx)
  select.list <- rep(select, nx)
  fp.mpos <- attr(Terms, "specials")$fp
  fp.pos <- fp.mpos -1
  fp.def <- m[, fp.mpos, drop = FALSE]
  df.list[fp.pos] <- unlist(lapply(fp.def, attr, "df"))
  
  alpha.list[fp.pos] <- unlist(lapply(fp.def, attr, "alpha"))
  alpha.list[sapply(alpha.list, is.na)] <- alpha
  select.list[fp.pos] <- unlist(lapply(fp.def, attr, "select"))
  select.list[sapply(select.list, is.na)] <- select
  
  fit <- mfpmi.fit(formula, mids, df.list, alpha.list, select.list, 
                   verbose = verbose, family = family, maxits = maxits)
  attr(fit, "class") <- c("mfpmi", "mira", "matrix")
}


mfpmi.fit <- function(formula, mids, family, alpha, select, verbose = TRUE, maxits = 20, ...){
  int <- as.numeric(!cox) # intercept
  nx <- ncol(x) - int
  nobs <- nrow(x)
  x.names <- dimnames(x)[[2]][int + seq(nx)]
  #    if(sum((df == 1 | df == 2 | df == 4), na.rm=TRUE) != nx)
  if(sum(df %in% c(1,2,4), na.rm=TRUE) != nx)
    stop("df is invalid")
  if(sum((alpha > 0 & alpha <= 1), na.rm=TRUE) != nx)
    stop("alpha is invalid")
  if(sum((select > 0 & select <= 1), na.rm=TRUE) != nx) stop(
    "select is invalid")    
  if(is.null(xnames)) xnames <- x.names
  #
  # Step 1: Order variables by LR test
  #  via one-step backward selection
  #
  x.order <- fp.order(formula, mids, family = family, ...)
  x.names <- x.names[x.order]
  df <- df[x.order]
  alpha <- alpha[x.order]
  select <- select[x.order]
  scaling <- scaling[x.order]
}


formula <- y ~ fp(x, df = 4, alpha = 0.1) + fp(z, df = 2, select = 0.5)
mids <- imp2
alpha <- 0.05
select <- 1

formula <- y ~ I(x^2) + I(log(z))

add <- as.formula(paste0("~ . + ", fpmi.gen("x", c(1, 1))))
update.formula(formula, add)

fpmi.fit <- function(formula, mids, variable, df, ...){
  
  add <- function(formula, term){
    term_as_form <- as.formula(paste0("~ . + ", term))
    update.formula(formula, term_as_form)
  }
  
  pwrs <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
  npwrs <- length(pwrs)
  
  # fit the linear model
  
  linear <- add(formula, variable)
  fit_lin <- with(mids, glm(formula(format(linear)), ...))

  fit_fp1 <- NULL
  p_fp1 <- 1
  fit_fp2 <- NULL
  p_fp2 <- 1
  
  if(df > 1){
    # Find best single power transformation 
    for(i in 1:npwrs) {
      
      if(i != 1){ # (but ignore 1, which is the same as linear)
        fp1 <- add(formula, fpmi.gen(variable, i))
        fit1 <- with(mids, glm(formula(format(fp1)), ...))
        p <- D1(fit1)$result[, "P(>F)"]
        
        if(p < p_fp1) {
          fit_fp1 <- fit1
          p_fp1 <- p
        }
      }
      
      if(df == 4) {
        # Find best two power transformation
        j <- i
        while(j <= npwrs) {
          
          fp2 <- add(formula, fpmi.gen(variable, c(i, j)))
          fit2 <- with(mids, glm(formula(format(fp2)), ...))
          p <- D1(fit2)$result[, "P(>F)"]
          
          if(p < p_fp2) {
            fit_fp2 <- fit2
            p_fp2 <- p
          }
          
          j <- j + 1
        }
      }
    }
    
  }
  
  res <- list(lin = fit_lin, fp1 = fit_fp1, fp2 = fit_fp2)
  res
}

fpmi.gen <- function(variable, pwrs){
  
  if(length(pwrs) > 2){
    stop("Only powers of degree 2 or less are supported. Please provide a maximum of 2 powers.")
  }
  
  gen <- function(pwr){
    if(pwr == 0){
      paste0("log(", variable, ")")
    } else {
      paste0(variable, "^", pwr)
    }
  }
  
  wrap <- function(x)  paste0("I(", x, ")")
  
  pwr1 <- gen(pwrs[1])
  
  if(length(pwrs) == 2){
    pwr2 <- gen(pwrs[2])
    
    if(pwrs[1] == pwrs[2]){
      pwr2 <- paste0("log(", variable, ") * ", pwr2)
    }
    
    return(paste0(wrap(pwr1), " + ", wrap(pwr2)))
  }
  
  return(wrap(pwr1))
}


fpmi.sel <- function(){
  
}


fpmi.order <- function(formula, mids, ...){

  fit <- with(mids, glm(formula = formula(format(formula)), ...))
  
  xnames <- labels(terms(formula))
  ns <- length(xnames)
  p.value <- numeric(ns)
  
  for(i in xnames){
    nested_formula <- update(formula, as.formula(paste(". ~ . -", i)))
    nested_fit <- with(mids, glm(formula = formula(format(nested_formula)), family = gaussian))
    
    p.value[xnames == i] <- D1(fit, nested_fit)$result[, "P(>F)"]
  }
  
  x.order <- order(p.value)
  
  return(x.order)
}
  

pool(with(mids, glm(formula = y ~ fpmi.tr(x, 2), family = gaussian)))


fpmi.tr <- function(x, p, log_x = FALSE){

  mult <- log(x) ^ log_x
  
  if(p == 0){
    x <- log(x)
  } else {
    x <- x ^ p
  }
  
  return(x * mult)
}






