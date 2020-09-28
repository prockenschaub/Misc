#### Functions and Libraries ####
library(mice)
expit = function(x){
  return(exp(x)/(1+exp(x)))
}
#### Initial Setup for Runs ####
xbeta <- 0.45
zbeta <- 0.55
xzdelta <- 0.6
prx <- 0.75
prz_x0 <- 0.6
prz_x1 <- 0.533
R <- 1000
numsub <- 2000
trials <- 1
coef.fl <- matrix(NA, nrow=R, ncol=3)
SE.fl <- matrix(NA, nrow=R, ncol=3)
lower.fl <- matrix(NA, nrow=R, ncol=3)
upper.fl <- matrix(NA, nrow=R, ncol=3)
coef.cra <- matrix(NA, nrow=R, ncol=3)
SE.cra <- matrix(NA, nrow=R, ncol=3)
lower.cra <- matrix(NA, nrow=R, ncol=3)
upper.cra <- matrix(NA, nrow=R, ncol=3)
coef.st <- matrix(NA, nrow=R, ncol=3)
SE.st <- matrix(NA, nrow=R, ncol=3)
lower.st <- matrix(NA, nrow=R, ncol=3)
upper.st <- matrix(NA, nrow=R, ncol=3)
coef.y <- matrix(NA, nrow=R, ncol=3)
SE.y <- matrix(NA, nrow=R, ncol=3)
lower.y <- matrix(NA, nrow=R, ncol=3)
upper.y <- matrix(NA, nrow=R, ncol=3)
coef.al <- matrix(NA, nrow=R, ncol=3)
SE.al <- matrix(NA, nrow=R, ncol=3)
lower.al <- matrix(NA, nrow=R, ncol=3)
upper.al <- matrix(NA, nrow=R, ncol=3)
#### ####
#### Imputation For Loop Mechanism B ####
for(i in 1:R){
  # PR 2020-09-28 --------------------------------------------------- start
  #    Add seeds for reproducibility
  set.seed(i)
  seeds <- runif(6, min = 1, max = 2^20)
  set.seed(seeds[1])  
  # PR 2020-09-28 ----------------------------------------------------- end
  
  e <- rnorm(numsub, mean=0, sd=1)
  ## Generate x and z
  x <-rbinom(n=numsub, size=trials, prob=prx)
  z <- ifelse(x==1, rbinom(n=numsub, size=trials, prob=prz_x1),
              rbinom(n=numsub, size=trials, prob= prz_x0))
  ## Generate y and full dataframe
  y = xbeta*x + zbeta*z + xzdelta*x*z + e
  data.full <-data.frame(x,y,z)
  
  ## Generate missingness
  z.misstop <- z[1:(numsub/2)]
  z.missbot <- z[(numsub/2+1):numsub]
  z.misstop[rbinom(n=numsub/2, size=trials, prob=expit(-3+1.5*x[1:(numsub/2)])) ==0] <- NA
  z.miss <- c(z.misstop, z.missbot)
  y.misstop <- y[1:(numsub/2)]
  y.missbot <- y[(numsub/2+1):numsub]
  y.missbot[rbinom(n=numsub/2, size=trials, prob=expit(-2+2.5*z[((numsub/2)+1):numsub])) ==0] <- NA
  y.miss <- c(y.misstop, y.missbot)
  
  ## Generate missingness dataset
  data.missing <- data.frame(x,y.miss,z.miss)
  
  ## Analysis 1 Full data analysis ".fl"
  fit.fl = lm(y~x+z+x*z, data=data.full)
  
  est.fl <- summary(fit.fl)
  
  coef.fl[i,1] <- est.fl$coeff[2,1]
  coef.fl[i,2] <- est.fl$coeff[3,1]
  coef.fl[i,3] <- est.fl$coeff[4,1]
  
  SE.fl[i,1] <- est.fl$coeff[2,2]
  SE.fl[i,2] <- est.fl$coeff[3,2]
  SE.fl[i,3] <- est.fl$coeff[4,2]
  
  confint.fl = confint(fit.fl)
  lower.fl[i,1] <- confint.fl[2,1]
  lower.fl[i,2] <- confint.fl[3,1]
  lower.fl[i,3] <- confint.fl[4,1]
  
  upper.fl[i,1] <- confint.fl[2,2]
  upper.fl[i,2] <- confint.fl[3,2]
  upper.fl[i,3] <- confint.fl[4,2]
  
  ## Analysis 2 CRA ".cra"
  fit.cra <- with(data.missing, lm(y.miss~x+z.miss+x*z.miss))
  est.cra <- summary(fit.cra)
  
  coef.cra[i,1] <- est.cra$coeff[2,1]
  coef.cra[i,2] <- est.cra$coeff[3,1]
  coef.cra[i,3] <- est.cra$coeff[4,1]
  
  SE.cra[i,1] <- est.cra$coeff[2,2]
  SE.cra[i,2] <- est.cra$coeff[3,2]
  SE.cra[i,3] <- est.cra$coeff[4,2]
  
  confint.cra = confint(fit.cra)
  lower.cra[i,1] <- confint.cra[2,1]
  lower.cra[i,2] <- confint.cra[3,1]
  lower.cra[i,3] <- confint.cra[4,1]
  
  upper.cra[i,1] <- confint.cra[2,2]
  upper.cra[i,2] <- confint.cra[3,2]
  upper.cra[i,3] <- confint.cra[4,2]
  
  ## Analysis 3 Standard No interaction ".st"
  imp.st <- mice(data.missing, print=FALSE, m=5, maxit=10, seed = seeds[2])
  fit.st <- with(imp.st, lm(y.miss~x*z.miss))
  est.st <- pool(fit.st)
  
  coef.st[i,1] <- est.st$pooled$estimate[2]
  coef.st[i,2] <- est.st$pooled$estimate[3]
  coef.st[i,3] <- est.st$pooled$estimate[4]
  
  SE.st[i,1] <- summary(est.st)[2,3]
  SE.st[i,2] <- summary(est.st)[3,3]
  SE.st[i,3] <- summary(est.st)[4,3]
  
  lower.st[i,1] <- summary(est.st, conf.int = T)[2,7]
  lower.st[i,2] <- summary(est.st, conf.int = T)[3,7]
  lower.st[i,3] <- summary(est.st, conf.int = T)[4,7]
  
  upper.st[i,1] <- summary(est.st, conf.int = T)[2,8]
  upper.st[i,2] <- summary(est.st, conf.int = T)[3,8]
  upper.st[i,3] <- summary(est.st, conf.int = T)[4,8]
  
  ## Analysis 4 y interactions ".y"
  
  ini <- mice(data.missing, maxit=0, print=FALSE)
  
  myformulas <-ini$formulas
  myformulas$x <- x ~ y.miss + z.miss
  myformulas$y.miss <- y.miss ~ x + z.miss + x*z.miss
  myformulas$z.miss <- z.miss ~ x + y.miss
  
  # PR 2020-09-28 --------------------------------------------------- start
  #    Add visit sequence and hange into 2 steps to make identical 
  #    to the predmat approach
  myVisSeq = c("x", "z.miss", "y.miss")
  
  imp.y <- mice(data.missing, print = FALSE, m =5, maxit = 0,
                 formulas = myformulas, visitSequence = myVisSeq,
                 seed = seeds[3])
  imp.y <- mice.mids(imp.y, print = FALSE, maxit = 10, seed = seeds[4])
  # PR 2020-09-28 ----------------------------------------------------- end
  
  fit.y <- with(imp.y, lm(y.miss~x*z.miss))
  est.y <- pool(fit.y)
  
  coef.y[i,1] <- est.y$pooled$estimate[2]
  coef.y[i,2] <- est.y$pooled$estimate[3]
  coef.y[i,3] <- est.y$pooled$estimate[4]
  
  SE.y[i,1] <- summary(est.y)[2,3]
  SE.y[i,2] <- summary(est.y)[3,3]
  SE.y[i,3] <- summary(est.y)[4,3]
  
  lower.y[i,1] <- summary(est.y, conf.int = T)[2,7]
  lower.y[i,2] <- summary(est.y, conf.int = T)[3,7]
  lower.y[i,3] <- summary(est.y, conf.int = T)[4,7]
  
  upper.y[i,1] <- summary(est.y, conf.int = T)[2,8]
  upper.y[i,2] <- summary(est.y, conf.int = T)[3,8]
  upper.y[i,3] <- summary(est.y, conf.int = T)[4,8]
  
  ## Analysis 5 All Interations ".al"
  ini <- mice(data.missing, maxit=0, print=FALSE)
  
  myformulas <-ini$formulas
  
  myformulas$x <- x ~ y.miss + z.miss + y.miss*z.miss
  myformulas$y.miss <- y.miss ~ x + z.miss + x*z.miss
  myformulas$z.miss <- z.miss ~ x + y.miss + x*y.miss
  
  # PR 2020-09-28 --------------------------------------------------- start
  #    Add visit sequence and hange into 2 steps to make identical 
  #    to the predmat approach
  myVisSeq = c("x", "z.miss", "y.miss")
  imp.al <- mice(data.missing, print = FALSE, m = 5, maxit = 0,
                 formulas = myformulas, visitSequence = myVisSeq, 
                 seed = seeds[5])
  imp.al <- mice.mids(imp.al, print = FALSE, maxit = 10, seed = seeds[6])
  # PR 2020-09-28 ----------------------------------------------------- end
  
  fit.al <- with(imp.al, lm(y.miss~x*z.miss))
  est.al <- pool(fit.al)
  
  coef.al[i,1] <- est.al$pooled$estimate[2]
  coef.al[i,2] <- est.al$pooled$estimate[3]
  coef.al[i,3] <- est.al$pooled$estimate[4]
  
  SE.al[i,1] <- summary(est.al)[2,3]
  SE.al[i,2] <- summary(est.al)[3,3]
  SE.al[i,3] <- summary(est.al)[4,3]
  
  lower.al[i,1] <- summary(est.al, conf.int = T)[2,7]
  lower.al[i,2] <- summary(est.al, conf.int = T)[3,7]
  lower.al[i,3] <- summary(est.al, conf.int = T)[4,7]
  
  upper.al[i,1] <- summary(est.al, conf.int = T)[2,8]
  upper.al[i,2] <- summary(est.al, conf.int = T)[3,8]
  upper.al[i,3] <- summary(est.al, conf.int = T)[4,8]
  
  if(i%%(R/20)==0) print(paste0((i/R)*100,"%"))
}

#### Saving imputed Mice objects ####
save(coef.fl, file = "form/coef.fl.RData")
save(SE.fl, file = "form/SE.fl.RData")
save(lower.fl, file = "form/lower.fl.RData")
save(upper.fl, file = "form/upper.fl.RData")
save(coef.cra, file = "form/coef.cra.RData")
save(SE.cra, file = "form/SE.cra.RData")
save(lower.cra, file = "form/lower.cra.RData")
save(upper.cra, file = "form/upper.cra.RData")
save(coef.st, file = "form/coef.st.RData")
save(SE.st, file = "form/SE.st.RData")
save(lower.st, file = "form/lower.st.RData")
save(upper.st, file = "form/upper.st.RData")
save(coef.y, file = "form/coef.y.RData")
save(SE.y, file = "form/SE.y.RData")
save(lower.y, file = "form/lower.y.RData")
save(upper.y, file = "form/upper.y.RData")
save(coef.al, file = "form/coef.al.RData")
save(SE.al, file = "form/SE.al.RData")
save(lower.al, file = "form/lower.al.RData")
save(upper.al, file = "form/upper.al.RData")

