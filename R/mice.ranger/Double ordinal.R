library(mice)

### Interactie: exp(b12)=1.5, b12=0.4
or_or_small <- properties(nsim=200, n=1000, m=20, maxit=1, ntrees=100, randomsel=3,
  B=c(Intercept=0, x1B=.5, x2B=1, x3B=-1.1, x4B=-.5, x8B=.4, x9B=-.5, x1Bx2B=.4))
or_or_small

### Interactie: exp(b12)=3, b12=1.1
or_or_large <- properties(nsim=200, n=1000, m=20, maxit=1, ntrees=100, randomsel=3,
  B=c(Intercept=0, x1B=.5, x2B=1, x3B=-1, x4B=-.5, x8B=.5, x9B=-1, x1Bx2B=1.1))
or_or_large

### Create binomial data
createdata <- function(n=1000, B=c(Intercept=0, x1B=-2.5, x2B=3, x3B=2, x4B=1.5, x8B=-1, x9B=-1.5, x1Bx2B=-2.5))
{
  Param <- as.vector(B)                                                         # Parameters logistic regression
  x1 <- factor(sample(LETTERS[1:2], n, replace=T))                              # Sample data for factor with 3 levels
  x2 <- factor(sample(LETTERS[1:2], n, replace=T))
  x3 <- factor(sample(LETTERS[1:2], n, replace=T))
  x4 <- factor(sample(LETTERS[1:2], n, replace=T))
  x5 <- factor(sample(LETTERS[1:3], n, replace=T))
  x6 <- factor(sample(LETTERS[1:3], n, replace=T))
  x7 <- factor(sample(LETTERS[1:3], n, replace=T))
  x8 <- factor(sample(LETTERS[1:2], n, replace=T))
  x9 <- factor(sample(LETTERS[1:2], n, replace=T))
  x10<- factor(sample(LETTERS[1:3], n, replace=T))
  X  <- model.matrix(~x1+x2+x3+x4+x8+x9+x1*x2)                               # Create data frame with dummys in order to perform matrix multiplication with the parameters (B)
  binomial_probabilities <- plogis(X%*%Param) ## plogis(x) = exp(x)/(1+exp(x))  # Results in probability of 'success'
  y <- factor(rbinom(n,1,binomial_probabilities))                               # Random generation for the binomial distribution
  return(data.frame(y=y,x1=x1,x2=x2,x3=x3,x4=x4,x5=x5,x6=x6,x7=x7,x8=x8,x9=x9,x10=x10))# Returns a data frame with factors
}

### Create missings
makeMARmid <- function(data)
{
  level9  <- levels(data[,'x9'])                                                    # Create missingness according to variable 9
  level10 <- levels(data[,'x10'])
  r.mcar <- rbinom(sum(data[,'x9']==level9[1]&data[,'x10']==level10[1]), 1, prob=.4)
  data[data[,'x9']==level9[1]&data[,'x10']==level10[1],][r.mcar==1,1] <- NA
  r.mcar <- rbinom(sum(data[,'x9']==level9[1]&data[,'x10']==level10[2]), 1, prob=.7)
  data[data[,'x9']==level9[1]&data[,'x10']==level10[2],][r.mcar==1,1] <- NA
  r.mcar <- rbinom(sum(data[,'x9']==level9[1]&data[,'x10']==level10[3]), 1, prob=.4)
  data[data[,'x9']==level9[1]&data[,'x10']==level10[3],][r.mcar==1,1] <- NA
  r.mcar <- rbinom(sum(data[,'x9']==level9[2]&data[,'x10']==level10[1]), 1, prob=.4)
  data[data[,'x9']==level9[2]&data[,'x10']==level10[1],][r.mcar==1,1] <- NA
  r.mcar <- rbinom(sum(data[,'x9']==level9[2]&data[,'x10']==level10[2]), 1, prob=.7)
  data[data[,'x9']==level9[2]&data[,'x10']==level10[2],][r.mcar==1,1] <- NA
  r.mcar <- rbinom(sum(data[,'x9']==level9[2]&data[,'x10']==level10[3]), 1, prob=.4)
  data[data[,'x9']==level9[2]&data[,'x10']==level10[3],][r.mcar==1,1] <- NA
  return(data)
}

###------------------------MICE.IMPUTE.CTREE--------------------
mice.impute.ctree <- function(y, ry, x, ...)
{
  library(rpart)
  ## Prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
  ## Fit the tree
  fit <- rpart(yobs~., data=cbind(yobs,xobs),method="class", control=rpart.control(minbucket=5, cp=.0001))
  ## Predict the level-probabilities according to the fitted tree for 'xmis'
  ymismat <- predict(object=fit, newdata=xmis)                                  # Matrix with probabilities per level per n
  ## Random generation of levels to impute, accroding to obtained probabilities
  impute <- apply(ymismat, MARGIN=1, FUN=function(s) sample(colnames(ymismat),size=1, prob=s[])) # Vector with levels to impute per person in xmis
  return(as.matrix(as.vector(impute)))                                          # Returns matrix with level to impute for every n in xmis
}

###------------------------MICE.IMPUTE.Bagging--------------------
mice.impute.bag <- function(y, ry, x, mtry=dim(x)[2], ntrees=100, ...)
{
  library(randomForest)
  ## Prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
  ## Fit the forest
  fit <- randomForest(yobs~., data=cbind(yobs,xobs), ntree=ntrees, mtry=dim(xobs)[2], replace=TRUE, sampsize=length(yobs), nodesize=5)
  leafnr <- attr(predict(object=fit, newdata=xobs, nodes=T),"nodes")            # Matrix with leafnumbers for 'obs' per tree in the forest
  ## Predict nodes
  pred_nodes <- attr(predict(object=fit, newdata=xmis, nodes=T),"nodes")        # Matrix with predicted leafnumbers for xmis per tree in the forest
  ## Draw donors                                                                # List with donors for every xmis (according to matching leafnumbers) per tree
  donors <- sapply(1:ntrees, FUN=function(t) sapply(pred_nodes[,t], function(s) yobs[leafnr[,t]==s[]]))
  ## Random generation of levels to impute, accroding to obtained donors
  impute  <- apply(donors, MARGIN=1, FUN=function(s) sample(unlist(s),1))       # Vector with levels to impute per person in xmis
  return(as.matrix(as.vector(impute)))                                          # Returns matrix with level to impute for every n in xmis
}

###------------------------MICE.IMPUTE.RF--------------------
mice.impute.rf <- function(y, ry, x, ntrees=100, randomsel=(dim(x)[2])/3, ...)
{
  library(randomForest)
  ## Prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
  ## Fit the forest
  fit <- randomForest(yobs~., data=cbind(yobs,xobs), ntree=ntrees, mtry=randomsel, replace=TRUE, sampsize=length(yobs), nodesize=5)
  leafnr <- attr(predict(object=fit, newdata=xobs, nodes=T),"nodes")            # Matrix with leafnumbers for 'obs' per tree in the forest
  ## Predict nodes
  pred_nodes <- attr(predict(object=fit, newdata=xmis, nodes=T),"nodes")        # Matrix with predicted leafnumbers for xmis per tree in the forest
  ## Draw donors                                                                # List with donors for every xmis (according to matching leafnumbers) per tree
  donors <- sapply(1:ntrees, FUN=function(t) sapply(pred_nodes[,t], function(s) yobs[leafnr[,t]==s[]]))
  ## Random generation of levels to impute, accroding to obtained donors
  impute  <- apply(donors, MARGIN=1, FUN=function(s) sample(unlist(s),1))       # Vector with levels to impute per person in xmis
  return(as.matrix(as.vector(impute)))                                          # Returns matrix with level to impute for every n in xmis
}

###------------------------Simulation study--------------------
test.impute <- function(data, m=20, method=method, maxit=1, ntrees=ntrees, randomsel=randomsel)
{
  imp <- mice(as.data.frame(data),method=method,m=m,print=FALSE,maxit=maxit)
  fit <- with(imp, glm(y~x1+x2+x3+x4+x8+x9+x1*x2, family='binomial'))
  est <- pool(fit)
  tab <- summary(est)
  return(tab[,c("est","se","lo 95", "hi 95", "fmi", "lambda")])
}

simulate <- function(nsim=nsim, n=n, m=m,maxit=maxit, ntrees=ntrees, randomsel=randomsel, B=B)
{
#  set.seed(41872)
  logreg <- array(NA,dim=c(length(B),6,nsim))
  ctree <- array(NA,dim=c(length(B),6,nsim))
  bag <- array(NA,dim=c(length(B),6,nsim))
  rfor <- array(NA,dim=c(length(B),6,nsim))
  data <- apply(array(,dim=c(n=1000,11,nsim)),MARGIN=3, FUN=function(s) s <- makeMARmid(createdata(n=n, B=B)))
  logreg[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="logreg",maxit=maxit))
  ctree[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="ctree",maxit=maxit))
  bag[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="bag",maxit=maxit, ntrees=ntrees))
  rfor[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="rf",maxit=maxit, ntrees=ntrees, randomsel=randomsel))
  list(logreg, ctree, bag, rfor)
}

maketable <- function(s, B=B)
{
  true <- as.matrix(B,,1)
  bias <- round(rowMeans(apply(s[,1,],2,function(s) s-true)),3)
  isinint <- apply(s[,3,],2,function(s) s < true) & apply(s[,4,],2,function(s) true < s)
  cov <- round(rowMeans(isinint),3)
  intwidth <- s[,4,] - s[,3,]
  aiw <- round(rowMeans(intwidth),3)
  fmi <- round(rowMeans(s[,5,]),3)
  lambda <- round(rowMeans(s[,6,]),3)
  prop_table <- cbind(bias,cov,aiw,fmi,lambda)
  rownames(prop_table) <- names(coef(glm(y~x1+x2+x3+x4+x8+x9+x1*x2, data=createdata(n=100), family='binomial')))
  colnames(prop_table) <- c("Bias", "Coverage", "CI width", "fmi", "Lambda")
  return(prop_table)
}

properties <- function(nsim=nsim, n=n, m=m, maxit=maxit, ntrees=ntrees, randomsel=randomsel, B=B)
{
  res <- simulate(nsim,n=n,m=m,maxit=maxit,ntrees=ntrees, randomsel=randomsel, B=B)
  list(maketable(res[[1]], B=B), maketable(res[[2]], B=B), maketable(res[[3]], B=B), maketable(res[[4]], B=B))
}
