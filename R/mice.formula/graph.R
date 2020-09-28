require(ggplot2)
## Load all data objects
library(glue)
path <- "form"

load(glue("{path}/coef.al.RData"))
load(glue("{path}/coef.y.RData"))
load(glue("{path}/coef.st.RData"))
load(glue("{path}/coef.cra.RData"))
load(glue("{path}/coef.fl.RData"))
load(glue("{path}/SE.al.RData"))
load(glue("{path}/SE.y.RData"))
load(glue("{path}/SE.st.RData"))
load(glue("{path}/SE.cra.RData"))
load(glue("{path}/SE.fl.RData"))
load(glue("{path}/lower.al.RData"))
load(glue("{path}/lower.y.RData"))
load(glue("{path}/lower.st.RData"))
load(glue("{path}/lower.cra.RData"))
load(glue("{path}/lower.fl.RData"))
load(glue("{path}/upper.al.RData"))
load(glue("{path}/upper.y.RData"))
load(glue("{path}/upper.st.RData"))
load(glue("{path}/upper.cra.RData"))
load(glue("{path}/upper.fl.RData"))
# Confidence Interval Width X
mean(upper.fl[,1] - lower.fl[,1])
mean(upper.cra[,1] - lower.cra[,1])
mean(upper.st[,1] - lower.st[,1])
mean(upper.y[,1] - lower.y[,1])
mean(upper.al[,1] - lower.al[,1])
# Confidence Interval Width Z
mean(upper.fl[,2] - lower.fl[,2])
mean(upper.cra[,2] - lower.cra[,2])
mean(upper.st[,2] - lower.st[,2])
mean(upper.y[,2] - lower.y[,2])
mean(upper.al[,2] - lower.al[,2])
# Confidence Interval Width XZ
mean(upper.fl[,3] - lower.fl[,3])
mean(upper.cra[,3] - lower.cra[,3])
mean(upper.st[,3] - lower.st[,3])
mean(upper.y[,3] - lower.y[,3])
mean(upper.al[,3] - lower.al[,3])
xbeta <- 0.45
zbeta <- 0.55
xzdelta <- 0.6
R = 1000
qmin = function(X){
  return(quantile(X, .05))
}
qmax = function(X){
  return(quantile(X, .95))
}
coverage.fl <- matrix(NA, nrow=R, ncol=3)
coverage.cra <- matrix(NA, nrow=R, ncol=3)
coverage.st <- matrix(NA, nrow=R, ncol=3)
coverage.y <- matrix(NA, nrow=R, ncol=3)
coverage.al <- matrix(NA, nrow=R, ncol=3)
for(i in 1:R){
  coverage.fl[i,1] <- (lower.fl[i,1] < xbeta & upper.fl[i,1] > xbeta)
  coverage.fl[i,2] <- (lower.fl[i,2] < zbeta & upper.fl[i,2] > zbeta)
  coverage.fl[i,3] <- (lower.fl[i,3] < xzdelta & upper.fl[i,3] > xzdelta)
  coverage.cra[i,1] <- (lower.cra[i,1] < xbeta & upper.cra[i,1] > xbeta)
  coverage.cra[i,2] <- (lower.cra[i,2] < zbeta & upper.cra[i,2] > zbeta)
  coverage.cra[i,3] <- (lower.cra[i,3] < xzdelta & upper.cra[i,3]> xzdelta)
  coverage.st[i,1] <- (lower.st[i,1] < xbeta & upper.st[i,1] > xbeta)
  coverage.st[i,2] <- (lower.st[i,2] < zbeta & upper.st[i,2] > zbeta)
  coverage.st[i,3] <- (lower.st[i,3] < xzdelta & upper.st[i,3] > xzdelta)
  coverage.y[i,1] <- (lower.y[i,1] < xbeta & upper.y[i,1] > xbeta)
  coverage.y[i,2] <- (lower.y[i,2] < zbeta & upper.y[i,2] > zbeta)
  coverage.y[i,3] <- (lower.y[i,3] < xzdelta & upper.y[i,3] > xzdelta)
  coverage.al[i,1] <- (lower.al[i,1] < xbeta & upper.al[i,1] > xbeta)
  coverage.al[i,2] <- (lower.al[i,2] < zbeta & upper.al[i,2] > zbeta)
  coverage.al[i,3] <- (lower.al[i,3] < xzdelta & upper.al[i,3] > xzdelta)
  print(paste0((i/R)*100,"%"))
}
# Setup for X data
est.x <- as.data.frame(c(coef.fl[,1], coef.cra[,1], coef.st[,1],
                         coef.y[,1], coef.al[,1]))
names(est.x) = c("est.x")
summary(est.x)
se.x <- as.data.frame(c(SE.fl[,1], SE.cra[,1], SE.st[,1],
                        SE.y[,1], SE.al[,1]))
names(se.x) = c("se.x")
summary(se.x)
cov.x <- as.data.frame(c(coverage.fl[,1], coverage.cra[,1],
                         coverage.st[,1], coverage.y[,1], coverage.al[,1]))
names(cov.x) = c("cov.x")
summary(cov.x)
method.x <- as.data.frame(c(rep("Full Data", R), rep("Complete Records", R), rep("MI (Standard)", R), rep("MI (y interactions)", R), rep("MI (all interactions)", R)))
names(method.x) <- c("method.x")
summary(method.x)
covp.x <- as.data.frame(c(rep(mean(coverage.fl[,1]),R),
                          rep(mean(coverage.cra[,1]),R), 
                          rep(mean(coverage.st[,1]),R),
                          rep(mean(coverage.y[,1]),R), 
                          rep(mean(coverage.al[,1]),R)))
names(covp.x) <- c("covp.x")
summary(covp.x)
datx <- data.frame(est.x, se.x, cov.x, method.x, covp.x)
## X data
fl <- subset(datx$est.x, datx$method.x == "Full Data")
fl1 <- subset(fl, fl >= qmin(fl) & fl <= qmax(fl))
cra <- subset(datx$est.x, datx$method.x == "Complete Records")
cra1 <- subset(cra, cra >= qmin(cra) & cra <= qmax(cra))
st <- subset(datx$est.x, datx$method.x == "MI (Standard)")
st1 <- subset(st, st >= qmin(st) & st <= qmax(st))
yi <- subset(datx$est.x, datx$method.x == "MI (y interactions)")
yi1 <- subset(yi, yi >= qmin(yi) & yi <= qmax(yi))
al <- subset(datx$est.x, datx$method.x == "MI (all interactions)")
al1 <- subset(al, al >= qmin(al) & al <= qmax(al))
meanx <- as.data.frame(c(fl, cra, st, yi, al))
mean(fl)
mean(cra)
mean(st)
mean(yi)
mean(al)
meanx <- c(rep(mean(fl), R*.9), rep(mean(cra), R*.9),
           rep(mean(st), R*.9), rep(mean(yi), R*.9), rep(mean(al), R*.9))
plotx <-as.data.frame(c(fl1, cra1, st1, yi1, al1))
names(plotx) <- c("plotx")
mx <- as.data.frame(c(rep("Full Data", R*.9), rep("Complete Records", R*.9), rep("MI (Standard)", R*.9), rep("MI (y interactions)", R*.9), rep("MI (all interactions)", R*.9)))
names(mx) <- c("mx")
cperx <- as.data.frame(c(rep(mean(coverage.fl[,1]),R*.9),
                         rep(mean(coverage.cra[,1]),R*.9),
                         rep(mean(coverage.st[,1]),R*.9), rep(mean(coverage.y[,1]),R*.9),
                         rep(mean(coverage.al[,1]),R*.9)))
names(cperx) <- c("cperx")
plx <- data.frame(plotx, mx, cperx, meanx)
## X Box Plots
p <- ggplot(plx, aes(x = mx, y = plotx))
p + geom_boxplot(outlier.shape = NA, width=0) +
  geom_hline(yintercept=xbeta) +
  geom_point(aes(x=5, y=meanx[1]), colour="blue") +
  geom_point(aes(x=4, y=meanx[901]), colour="blue") +
  geom_point(aes(x=3, y=meanx[1801]), colour="blue") +
  geom_point(aes(x=2, y=meanx[2701]), colour="blue") +
  geom_point(aes(x=1, y=meanx[3601]), colour="blue") +
  scale_x_discrete(limits = c("MI (all interactions)", "MI (y interactions)", "MI (Standard)", "Complete Records", "Full Data")) +
  scale_y_continuous(breaks = c(-0.25, 0.00, 0.25, 0.50, 0.75, 1.00, 1.25)) +
  coord_flip(ylim = c(-0.5, 1.5)) +
  ggtitle("<Name of Method>") +
  theme_classic() + 
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border = element_rect(colour = "black", fill=NA,
                                    size=2)) +
  ylab("Estimated Values for X") +
  xlab("Missing Data Method") # +
  # annotate(geom = "text", x=5, y=1.5, label=round(mean(coverage.fl[,1])*100), size=4.5) +
  # annotate(geom = "text", x=4, y=1.5, label=round(mean(coverage.cra[,1])*100), size=4.5) +
  # annotate(geom = "text", x=3, y=1.5, label=round(mean(coverage.st[,1])*100), size=4.5) +
  # annotate(geom = "text", x=2, y=1.5, label=round(mean(coverage.y[,1])*100), size=4.5) +
  # annotate(geom = "text", x=1, y=1.5, label=round(mean(coverage.al[,1])*100), size=4.5) +
  # annotate(geom = "text", x=5, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=4, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=3, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=2, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=1, y=1.59, label="%", size=4.5) +
  # geom_vline(xintercept = 4.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 3.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 2.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 1.5, fill = NA, size = .5)

# Setup for z data
est.z <- as.data.frame(c(coef.fl[,2], coef.cra[,2], coef.st[,2],
                         coef.y[,2], coef.al[,2]))
names(est.z) = c("est.z")
summary(est.z)
se.z <- as.data.frame(c(SE.fl[,2], SE.cra[,2], SE.st[,2],
                        SE.y[,2], SE.al[,2]))
names(se.z) = c("se.z")
summary(se.z)
cov.z <- as.data.frame(c(coverage.fl[,2], coverage.cra[,2],
                         coverage.st[,2], coverage.y[,2], coverage.al[,2]))
names(cov.z) = c("cov.z")
summary(cov.z)
method.z <- as.data.frame(c(rep("Full Data", R), rep("Complete Records", R), rep("MI (Standard)", R), rep("MI (y interactions)", R), rep("MI (all interactions)", R)))
names(method.z) <- c("method.z")
summary(method.z)
covp.z <- as.data.frame(c(rep(mean(coverage.fl[,2]),R),
                          rep(mean(coverage.cra[,2]),R), rep(mean(coverage.st[,2]),R),
                          rep(mean(coverage.y[,2]),R), rep(mean(coverage.al[,2]),R)))
names(covp.z) <- c("covp.z")
summary(covp.z)
datz <- data.frame(est.z, se.z, cov.z, method.z, covp.z)
## z data
fl <- subset(datz$est.z, datz$method.z == "Full Data")
fl1 <- subset(fl, fl >= qmin(fl) & fl <= qmax(fl))
cra <- subset(datz$est.z, datz$method.z == "Complete Records")
cra1 <- subset(cra, cra >= qmin(cra) & cra <= qmax(cra))
st <- subset(datz$est.z, datz$method.z == "MI (Standard)")
st1 <- subset(st, st >= qmin(st) & st <= qmax(st))
yi <- subset(datz$est.z, datz$method.z == "MI (y interactions)")
yi1 <- subset(yi, yi >= qmin(yi) & yi <= qmax(yi))
al <- subset(datz$est.z, datz$method.z == "MI (all interactions)")
al1 <- subset(al, al >= qmin(al) & al <= qmax(al))
meanz <- as.data.frame(c(fl, cra, st, yi, al))
mean(fl)
mean(cra)
mean(st)
mean(yi)
mean(al)
meanz <- c(rep(mean(fl), R*.9), rep(mean(cra), R*.9),
           rep(mean(st), R*.9), rep(mean(yi), R*.9), rep(mean(al), R*.9))
plotz <-as.data.frame(c(fl1, cra1, st1, yi1, al1))
names(plotz) <- c("plotz")
mz <- as.data.frame(c(rep("Full Data", R*.9), rep("Complete Records", R*.9), rep("MI (Standard)", R*.9), rep("MI (y interactions)", R*.9), rep("MI (all interactions)", R*.9)))
names(mz) <- c("mz")
cperz <- as.data.frame(c(rep(mean(coverage.fl[,2]),R*.9),
                         rep(mean(coverage.cra[,2]),R*.9),
                         rep(mean(coverage.st[,2]),R*.9), rep(mean(coverage.y[,2]),R*.9),
                         rep(mean(coverage.al[,2]),R*.9)))
names(cperz) <- c("cperz")
plz <- data.frame(plotz, mz, cperz, meanz)
## Z Box Plots
p <- ggplot(plz, aes(x = mz, y = plotz))
p + geom_boxplot(outlier.shape = NA, width=0) +
  geom_hline(yintercept=zbeta) + 
  geom_point(aes(x=5, y=meanz[1]), colour="blue") +
  geom_point(aes(x=4, y=meanz[901]), colour="blue") +
  geom_point(aes(x=3, y=meanz[1801]), colour="blue") +
  geom_point(aes(x=2, y=meanz[2701]), colour="blue") +
  geom_point(aes(x=1, y=meanz[3601]), colour="blue") +
  scale_x_discrete(limits = c("MI (all interactions)", "MI (y interactions)", "MI (Standard)", "Complete Records", "Full Data")) +
  scale_y_continuous(breaks = c(-0.25, 0.00, 0.25, 0.50, 0.75, 1.00, 1.25)) +
  coord_flip(ylim = c(-0.5, 1.5)) +
  ggtitle("<Name of Method>") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab("Estimated Values for Z") +
  xlab("Missing Data Method") # +
  # annotate(geom = "text", x=5, y=1.5,
  #        label=round(mean(coverage.fl[,2])*100), size=4.5) +
  # annotate(geom = "text", x=4, y=1.5,
  #          label=round(mean(coverage.cra[,2])*100), size=4.5) +
  # annotate(geom = "text", x=3, y=1.5,
  #          label=round(mean(coverage.st[,2])*100), size=4.5) +
  # annotate(geom = "text", x=2, y=1.5,
  #          label=round(mean(coverage.y[,2])*100), size=4.5) +
  # annotate(geom = "text", x=1, y=1.5,
  #          label=round(mean(coverage.al[,2])*100), size=4.5) +
  # annotate(geom = "text", x=5, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=4, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=3, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=2, y=1.59, label="%", size=4.5) +
  # annotate(geom = "text", x=1, y=1.59, label="%", size=4.5) +
  # geom_vline(xintercept = 4.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 3.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 2.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 1.5, fill = NA, size = .5)
# Setup for xz data
est.xz <- as.data.frame(c(coef.fl[,3], coef.cra[,3],
                          coef.st[,3], coef.y[,3], coef.al[,3]))
names(est.xz) = c("est.xz")
summary(est.xz)
se.xz <- as.data.frame(c(SE.fl[,3], SE.cra[,3], SE.st[,3],
                         SE.y[,3], SE.al[,3]))
names(se.xz) = c("se.xz")
summary(se.xz)
cov.xz <- as.data.frame(c(coverage.fl[,3], coverage.cra[,3],
                          coverage.st[,3], coverage.y[,3], coverage.al[,3]))
names(cov.xz) = c("cov.xz")
summary(cov.xz)
method.xz <- as.data.frame(c(rep("Full Data", R), rep("Complete Records", R), rep("MI (Standard)", R), rep("MI (y interactions)", R), rep("MI (all interactions)", R)))
names(method.xz) <- c("method.xz")
summary(method.xz)
covp.xz <- as.data.frame(c(rep(mean(coverage.fl[,3]),R),
                           rep(mean(coverage.cra[,3]),R), rep(mean(coverage.st[,3]),R),
                           rep(mean(coverage.y[,3]),R), rep(mean(coverage.al[,3]),R)))
names(covp.xz) <- c("covp.xz")
summary(covp.xz)
datxz <- data.frame(est.xz, se.xz, cov.xz, method.xz, covp.xz)
## xz data
fl <- subset(datxz$est.xz, datxz$method.xz == "Full Data")
fl1 <- subset(fl, fl >= qmin(fl) & fl <= qmax(fl))
cra <- subset(datxz$est.xz, datxz$method.xz == "Complete Records")
cra1 <- subset(cra, cra >= qmin(cra) & cra <= qmax(cra))
st <- subset(datxz$est.xz, datxz$method.xz == "MI (Standard)")
st1 <- subset(st, st >= qmin(st) & st <= qmax(st))
yi <- subset(datxz$est.xz, datxz$method.xz == "MI (y interactions)")
yi1 <- subset(yi, yi >= qmin(yi) & yi <= qmax(yi))
al <- subset(datxz$est.xz, datxz$method.xz == "MI (all interactions)")
al1 <- subset(al, al >= qmin(al) & al <= qmax(al))
meanxz <- as.data.frame(c(fl, cra, st, yi, al))
mean(fl)
mean(cra)
mean(st)
mean(yi)
mean(al)
meanxz <- c(rep(mean(fl), R*.9), rep(mean(cra), R*.9),
            rep(mean(st), R*.9), rep(mean(yi), R*.9), rep(mean(al), R*.9))
plotxz <-as.data.frame(c(fl1, cra1, st1, yi1, al1))
names(plotxz) <- c("plotxz")
mxz <- as.data.frame(c(rep("Full Data", R*.9), rep("Complete Records", R*.9), rep("MI (Standard)", R*.9), rep("MI (y interactions)", R*.9), rep("MI (all interactions)", R*.9)))
names(mxz) <- c("mxz")
cperxz <- as.data.frame(c(rep(mean(coverage.fl[,3]),R*.9),
                          rep(mean(coverage.cra[,3]),R*.9),
                          rep(mean(coverage.st[,3]),R*.9), rep(mean(coverage.y[,3]),R*.9),
                          rep(mean(coverage.al[,3]),R*.9)))
names(cperxz) <- c("cperxz")
plxz <- data.frame(plotxz, mxz, cperxz, meanxz)
## xz Box Plots
p <- ggplot(plxz, aes(x = mxz, y = plotxz))
p + geom_boxplot(outlier.shape = NA, width=0) +
  geom_hline(yintercept=xzdelta) +
  geom_point(aes(x=5, y=meanxz[1]), colour="blue")+
  geom_point(aes(x=4, y=meanxz[901]), colour="blue") +
  geom_point(aes(x=3, y=meanxz[1801]), colour="blue") +
  geom_point(aes(x=2, y=meanxz[2701]), colour="blue") +
  geom_point(aes(x=1, y=meanxz[3601]), colour="blue") +
  scale_x_discrete(limits = c("MI (all interactions)", "MI (y interactions)", "MI (Standard)", "Complete Records", "Full Data")) +
  scale_y_continuous(breaks = c(-0.25, 0.00, 0.25, 0.50, 0.75, 1.00, 1.25)) +
  coord_flip(ylim = c(-0.5, 1.5)) +
  ggtitle("<Name of Method>") +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border = element_rect(colour = "black", fill=NA, size=2)) +
  ylab("Estimated Values for xz") +
  xlab("Missing Data Method") +
  annotate(geom = "text", x=5, y=1.4, label=round(mean(coverage.fl[,3])*100), size=4.5) +
  annotate(geom = "text", x=4, y=1.4, label=round(mean(coverage.cra[,3])*100), size=4.5) +
  annotate(geom = "text", x=3, y=1.4, label=round(mean(coverage.st[,3])*100), size=4.5) +
  annotate(geom = "text", x=2, y=1.4, label=round(mean(coverage.y[,3])*100), size=4.5) +
  annotate(geom = "text", x=1, y=1.4, label=round(mean(coverage.al[,3])*100), size=4.5) +
  annotate(geom = "text", x=5, y=1.5, label="%", size=4.5) +
  annotate(geom = "text", x=4, y=1.5, label="%", size=4.5) +
  annotate(geom = "text", x=3, y=1.5, label="%", size=4.5) +
  annotate(geom = "text", x=2, y=1.5, label="%", size=4.5) +
  annotate(geom = "text", x=1, y=1.5, label="%", size=4.5) #+
  # geom_vline(xintercept = 4.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 4.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 3.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 2.5, fill = NA, size = .5) +
  # geom_vline(xintercept = 1.5, fill = NA, size = .5)

ggsave("xz_by_method.png", dpi = 300, width = 7, height = 4)

