
dir_root <- "mfpmi" 
dir_test <- "testdata"

source(file.path(dir_root, "icet.R"))

# Simulated example -------------------------------------------------------


set.seed(3472)
n <- 10000
data <- mvtnorm::rmvnorm(n = n, 
                         mean = c(0, 3), 
                         sigma = matrix(c(1, 0.7, 0.7, 1), ncol = 2))
data <- as.data.frame(data)
colnames(data) <- c("y", "xp")


p <- 2
data$x <- data$xp ^ (1 / p)
data$z <- runif(n = n)

probs <- ifelse(data$y <= 0, 0.2, 0.6)
data$x <- ifelse(runif(nrow(data)) <= probs, NA, data$x)

imp <- icet(data, "x", c("y", "z"), m = 5, maxit = 10, verbose = TRUE)
imp2 <- icet.mids(data, "x", c("y", "z"), m = 5, maxit = 10, verbose = TRUE)


std_imp <- mice(data, maxit = 0, print = FALSE)
pred <- std_imp$predictorMatrix
pred[, "xp"] <- 0
pred["xp", ] <- 0

meth <- std_imp$method

std_imp <- mice(data, m = 5, maxit = 10, pred = pred, meth = meth, print = FALSE)


res <- with(imp2, lm(y ~ x + z))
std_res <- with(std_imp, lm(y ~ x))


# Real-world example ------------------------------------------------------


df <- read.csv(file.path(dir_root, dir_test, "gbsg_br_ca.csv"))
df[runif(nrow(df)) < 0.2, "pgr"] <- NA
df$pgr <- df$pgr + 0.1
df[runif(nrow(df)) < 0.2, "nodes"] <- NA
