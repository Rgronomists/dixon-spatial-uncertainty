# R code extracted from uncert.Rmd

library(spup)    # key package: SPatial Uncertainty Propagation
library(raster)  # raster graphics
library(purrr)   # library with similar functionality to dplyr

# plot the C and N maps
data(OC, OC_sd, TN, TN_sd)
par(mfrow=c(2,2), mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
plot(OC); legend('top', bty='n', legend='mean C', inset=-0.05)
plot(OC_sd); legend('top', bty='n', legend='sd C', inset=-0.05)
plot(TN); legend('top', bty='n', legend='mean total N', inset=-0.05)
plot(TN_sd); legend('top', bty='n', legend='sd total N', inset=-0.05)

# calculate and plot estimated map of C/N ratio 
CN <- OC/TN
par(mfrow=c(1,1))
plot(CN); legend('top', bty='n', legend='Estimated C/N', inset=-0.05)

# simpler uncertainty propagation, regression of meat pH on log time

meat <- data.frame(
  time = rep(c(1,2,4,6,8), rep(2,5)),
  pH = c(7.02, 6.93, 6.42, 6.51, 6.07, 5.99, 
    5.59, 5.80, 5.51, 5.36)
  )
meat$logtime <- log(meat$time)
with(meat, plot(logtime, pH, pch=19, col=4) )
meat.lm <- lm(pH ~ logtime, data=meat)
coef(meat.lm)
abline(coef(meat.lm))
abline(h=5.7, lty=3)

# question is when does average pH = 5.7?
# parametric bootstratp "by hand"

library(mvtnorm)
mub <- coef(meat.lm)   # estimates => mean parameter values
vcb <- vcov(meat.lm)   # variance-covariance matrix
rb <- rmvnorm(1000, mean=mub, sigma=vcb)
head(rb)
b0 <- rb[,1]
b1 <- rb[,2]
estT <- exp((5.7 - b0)/b1)
hist(estT, main='')
sd(estT)
# se of estimate, log Time scale
exp(quantile(estT, c(0.025, 0.975)))
# percentile bootstrap 95% CI for Time

# using propagate
library(propagate)
dimnames(rb)[[2]] <- c('b0','b1')
# set column names so we propagate() knows what we are talking about
That <- expression(exp((5.7-b0)/b1))
estTp <- propagate(expr=That, data=rb, type='sim',  nsim=100000)
estTp

# Spatial uncertainty

library(spup)
# describe the spatial correlation for each variable
Ccrm <- makecrm(acf0 = 0.6, range = 5000, model = "Sph")
Ncrm <- makecrm(acf0 = 0.4, range = 5000, model = "Sph")
par(mfrow=c(1,2), mar=c(3,3,0,0)+0.2, mgp=c(2,0.8,0))
plot(Ccrm); legend('top', bty='n', legend='C')
plot(Ncrm); legend('top', bty='n', legend='N')

Cum <- defineUM(TRUE, distribution = "norm", distr_param = c(OC, OC_sd), crm = Ccrm, id = "OC")
Num <- defineUM(TRUE, distribution = "norm", distr_param = c(TN, TN_sd), crm = Ncrm, id = "TN")
CNum <- defineMUM(UMlist = list(Cum, Num), 
  cormatrix = matrix(c(1,0.7,0.7,1), nrow=2, ncol=2))

CNrand2 <- genSample(UMobject = CNum, n = 2, samplemethod = "ugs", nmax = 20, asList = FALSE)
plot(CNrand2)

CNratio <- function(OC, TN) {
  OC/TN
}

Nmc <- 100
CNrand100 <- genSample(UMobject = CNum, n = Nmc, samplemethod = "ugs",  
  nmax = 20, asList = T)
ratiorand <- spup::propagate(realizations=CNrand100, model=CNratio, n=Nmc)

ratio <- stack(ratiorand)
names(ratio) <- paste("CNratio", 1:nlayers(ratio), sep='.')
par(mfrow=c(1,1))
plot(ratio[[1:4]])

ratiomean <- mean(ratio)
ratiosd <- calc(ratio, fun=sd)
par(mfrow=c(1,2))
plot(ratiomean, main='Mean C/N ratio')
plot(ratiosd, main='sd of C/N ratio')

bigCN <- function(r) {
  sum(as.matrix(r) > 20, na.rm=T)
}
bigratio <- sapply(ratiorand, bigCN)
hist(bigratio, main='', xlab='Number of cells with c/N > 20')

