library(MASS)

ptm <- proc.time()
x8x9_small <- properties(nsim=200,n=1000,m=20,maxit=1, ntrees=100, randomsel=3,
        B=c(beta0=0, beta1=.34, beta2=.34, beta3=.34, beta4=.34, beta5=.34, beta8=.13))
proc.time() - ptm
x8x9_small

ptm <- proc.time()
x8x9_medium <- properties(nsim=200,n=1000,m=20,maxit=1, ntrees=100, randomsel=3,
        B=c(beta0=0, beta1=.31, beta2=.31, beta3=.31, beta4=.31, beta5=.31, beta8=.37))
proc.time() - ptm
x8x9_medium

ptm <- proc.time()
x8x9_large <- properties(nsim=200,n=1000,m=20,maxit=1, ntrees=100, randomsel=3,
        B=c(beta0=0, beta1=.28, beta2=.28, beta3=.28, beta4=.28, beta5=.28, beta8=.57))
proc.time() - ptm
x8x9_large

###
createdata <- function(n=1000, B=c(beta0=0, beta1=.5, beta2=.5, beta3=.5, beta4=.5, beta5=.5, beta6=.5))
{
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
  Intercept <- as.vector(rep(0,n))
  X <- matrix(c(Intercept, x[,1], x[,2] , x[,3], x[,8], x[,9], (x[,8]*x[,9])),nrow=n)
  eps <- rnorm(n, mean=0, sd=1)
  y <- round((X%*%Param)+eps,2)
  return(data.frame(y=y, x1=x[,1], x2=x[,2], x3=x[,3], x4=x[,4], x5=x[,5], x6=x[,6], x7=x[,7], x8=x[,8], x9=x[,9], x10=x[,10]))
}

### Make 50% univariate missing
makeMARmid <- function(data)
{
  logistic <- function(x) exp(x)/(1+exp(x))
  x9 <- data[,"x9"]
  x10 <- data[,"x10"]
  for (i in c("y"))
  {
    p.marmid <- 1-logistic(-.7+abs((x9*x10)-mean(x9*x10)))                      # -.7 geeft 50%
    r.marmid <- rbinom(nrow(data),1,p.marmid)                                   # -3 geeft 10%
    data[r.marmid==0,i] <- NA                                                   # -1.6 geeft 30%
  }
  return(data)
}

###------------------------MICE.IMPUTE.RTREE--------------------
mice.impute.rtree <- function(y, ry, x, ...)
{
  library(rpart)
  ## prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
  ## fit the tree
  fit <- rpart(yobs~., data=cbind(yobs,xobs),method="anova", control=rpart.control(minbucket=5,cp=.0001))
  leafnr  <-floor(as.numeric(row.names(fit$frame[fit$where,])))                 # vector with leafnumbers for 'obs'
  pred_yobs <- predict(fit)                                                     # vector with predicted values for 'obs'
  ## predict nodes
  ymismat <- predict(object=fit, newdata=xmis)                                  # vector with predicted values for xmis
  pred_nodes <- sapply(ymismat,FUN=function(s) unique(leafnr[s[]==pred_yobs]))  # vector with predicted nodes for xmis
  ## draw donors
  donor <- sapply(pred_nodes,function(s) yobs[leafnr==s[]])                     # list with donors for every xmis
  ## impute
  impute  <- sapply(rownames(xmis),function(s) sample(donor[[s[]]],1))          # vector with imputed values for xmis
  return(as.matrix(as.vector(impute)))
}

###------------------------MICE.IMPUTE.CTREE--------------------
mice.impute.ctree <- function(y, ry, x, ...)
{
  library(rpart)
  ## prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
  ## fit the tree model
  fit <- rpart(yobs~., data=cbind(yobs,xobs),method="class", control=rpart.control(minbucket=5, cp=.0001))
  leafnr  <-floor(as.numeric(row.names(fit$frame[fit$where,])))                 # vector with leafnumbers for 'obs'
  pred_yobs <- predict(fit)                                                     # vector with predicted values for 'obs'
  ## predict the nodes
  ymismat <- predict(object=fit, newdata=xmis)                                  # vector with predicted values for xmis
  pred_nodes <- sapply(ymismat,FUN=function(s) unique(leafnr[s[]==pred_yobs]))  # vector with predicted nodes for xmis
  ## draw donors
  donor <- sapply(pred_nodes,function(s) yobs[leafnr==s[]])                     # list with donors for every xmis
  ## impute
  impute  <- sapply(rownames(xmis),function(s) sample(donor[[s[]]],1))          # vector with imputed values for xmis
  return(as.matrix(as.vector(impute)))
}

###------------------------MICE.IMPUTE.Bagging--------------------
mice.impute.bag <- function(y, ry, x, mtry=dim(x)[2], ntrees=100, nodesize=5, ...)
{
  library(randomForest)
  ## prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
 ## fit one tree from the forest and return donors
  onetree <- function(xobs, xmis, yobs)
  {
    ## fit the tree
    fit <- randomForest(yobs~., data=cbind(yobs,xobs), ntree=1, mtry=dim(xobs)[2], replace=TRUE, type=regression, sampsize=length(yobs), nodesize=5)
    leafnr <- predict(object=fit, newdata=xobs, nodes=T)                        # matrix with leafnumbers for 'obs'
    pred_yobs <- predict(object=fit, newdata=xobs, nodes=F)                     # matrix with predicted values for 'obs'
    ## predict nodes
    pred_nodes <- predict(object=fit, newdata=xmis, nodes=T)                    # vector with predicted nodes for xmis
    ## draw donors
    donor <- sapply(pred_nodes,function(s) yobs[leafnr==s[]])                   # list with donors for every xmis
    return(donor)
  }
  ## fit the forest and combine the donors
  rflist <- sapply(1:(ntrees=ntrees), FUN=function(s) onetree(xobs, xmis, yobs)) # matrix of lists with donors for every tree per person
  ## Random generation of levels to impute, accroding to obtained donors
  impute  <- apply(rflist, MARGIN=1, FUN=function(s) sample(unlist(s),1))
  return(as.matrix(as.vector(impute)))
}

###------------------------MICE.IMPUTE.RF--------------------
mice.impute.rf <- function(y, ry, x, mtry=(dim(x)[2])/3, ntrees=100, randomsel=(dim(x)[2])/3, nodesize=5, ...)
{
  library(randomForest)
  ## prepare the data
  xobs <- x[ry,]
  xmis <- x[!ry,]
  yobs <- y[ry]
 ## fit one tree from the forest and return donors
  onetree <- function(xobs, xmis, yobs)
  {
    ## fit the tree
    fit <- randomForest(yobs~., data=cbind(yobs,xobs), ntree=1, mtry=randomsel, replace=TRUE, type=regression, sampsize=length(yobs), nodesize=5)
    leafnr <- predict(object=fit, newdata=xobs, nodes=T)                        # matrix with leafnumbers for 'obs'
    pred_yobs <- predict(object=fit, newdata=xobs, nodes=F)                     # matrix with predicted values for 'obs'
    ## predict nodes
    pred_nodes <- predict(object=fit, newdata=xmis, nodes=T)                    # vector with predicted nodes for xmis
    ## draw donors
    donor <- sapply(pred_nodes,function(s) yobs[leafnr==s[]])                   # list with donors for every xmis
    return(donor)
  }
  ## fit the forest and combine the donors
  rflist <- sapply(1:(ntrees=ntrees), FUN=function(s) onetree(xobs, xmis, yobs)) # matrix of lists with donors for every tree per person
  ## Random generation of levels to impute, accroding to obtained donors
  impute  <- apply(rflist, MARGIN=1, FUN=function(s) sample(unlist(s),1))
  return(as.matrix(as.vector(impute)))
}

###------------------------pilot--------------------
test.impute <- function(data, m=20, method="rtree", maxit=1, ntrees=ntrees, randomsel=randomsel)
{
  library(mice)
  imp <- mice(data,method=method,m=m,print=FALSE,maxit=maxit)
  fit <- with(imp, lm(y~x1+x2+x3+x8+x9+(x8*x9)))
  est <- pool(fit)
  tab <- summary(est)
  return(tab[,c("est","se","lo 95", "hi 95", "fmi", "lambda")])
}

simulate <- function(nsim=nsim, n=n, m=m,maxit=maxit, ntrees=ntrees, randomsel=randomsel, B=B)
{
#  set.seed(41872)
  pmm <- array(NA,dim=c(length(B),6,nsim))
  rtree <- array(NA,dim=c(length(B),6,nsim))
  bag <- array(NA,dim=c(length(B),6,nsim))
  rfor <- array(NA,dim=c(length(B),6,nsim))
  data <- apply(array(,dim=c(n=1000,11,nsim)),MARGIN=3, FUN=function(s) s <- makeMARmid(createdata(n=n, B=B)))
  pmm[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="pmm",maxit=maxit))
  rtree[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="rtree",maxit=maxit))
  bag[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="bag",maxit=maxit, ntrees=ntrees))
  rfor[,,] <- sapply(data, FUN=function(s) test.impute(s, m, method="rf",maxit=maxit, ntrees=ntrees, randomsel=randomsel))
  list(pmm, rtree, bag, rfor)
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
  rownames(prop_table) <- names(coef(lm(y~x1+x2+x3+x8+x9+(x8*x9), data=createdata(n=100))))
  colnames(prop_table) <- c("Bias", "Coverage", "CI width", "fmi", "Lambda")
  return(prop_table)
}

properties <- function(nsim=nsim, n=n, m=m, maxit=maxit, ntrees=ntrees, randomsel=randomsel, B=B)
{
  res <- simulate(nsim,n=n,m=m,maxit=maxit,ntrees=ntrees, randomsel=randomsel, B=B)
  list(maketable(res[[1]], B=B), maketable(res[[2]], B=B), maketable(res[[3]], B=B), maketable(res[[4]], B=B))
}