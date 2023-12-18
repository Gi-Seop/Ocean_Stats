

## Periodicity test, monthly MSL data - MMSL ======================
library(imputeTS)

idata <- read.csv("../data/Busan_MSL_data_1961_2022.csv")
idata <- idata[-1,]
## 1961년 자료는 분석에서 제외

## AMSL, MMSL - annual and monthly MSL data set,

YY <- idata$YEAR
nYY <- length(YY)
month <- seq(1962, 2022, length.out=61*12)
ndata <- nYY*12

MMSL <- as.vector(t(idata[,2:13]))
MMSL <- na_kalman(MMSL)

str(MMSL)

# linear model fitting
fit1 <- lm(MMSL ~ month)
plot(month, MMSL, type="b", xlab="YEAR",ylab="monthly MSL (cm)", cex=1.2)
abline(fit1, col="blue", lwd=2)
abline(h=mean(MMSL), col="red", lwd=2) 

rr <- fit1$residuals
info1 <- influence.measures(fit1)			## 영향자료 진단
info1_summary <- summary(info1)
ploc <- as.numeric(rownames(info1_summary))


plot(month,rr, type="b",xlab="YEAR",ylab = "De-trended M-MSL (cm)", cex=1.2)
abline(h=0,col="red", lwd=2) 
points(month[ploc],rr[ploc], col="red", pch=16, cex=1.2)
points(month[ploc],rr[ploc], col="red", pch=12, cex=2.0)


## install.packages("TSA")
library(TSA)
prd1 <- periodogram(rr)
str(prd1)

# periodogram: manual calculation
ndata <- length(rr)
azero <- (1/ndata)*sum(rr)
hn <- round(ndata/2)
ck <- matrix(0, nrow=hn, ncol=3) ## (ak, bk, Ik), k=1,2,..., hn(=n/2)

# iteration n/2
for (kk in 1:hn) {
  wk <- 2*pi*kk/ndata
  ak <- 0			## Fourier coefficients for cos function.
  bk <- 0			## Fourier coefficients for sin function.
  for (ii in 1:ndata)
  {
    ak <- ak + (2/ndata)*cos(wk*ii)*rr[ii]
    bk <- bk + (2/ndata)*sin(wk*ii)*rr[ii] 
  }
  ck[kk,1] <- ak
  ck[kk,2] <- bk
}
ck[,3] <- (ndata/2)*(ck[,1]^2 + ck[,2]^2)

# visualization for spectrum
plot(prd1$freq, log(prd1$spec), type = "l",
     cex.lab = 1.2, main = "TSA",
     xlab = "Frequency (1/s)", ylab = "Magnitude of Spectrum (m^2-s)")	
plot(seq(0, 0.5, length.out = hn), log(ck[,3]), type = "l",
     cex.lab = 1.2, main = "Manual",
     xlab = "Frequency (1/s)", ylab = "Magnitude of Spectrum (m^2-s)")


spectrum(rr)


# F-statistics

F_stat <- matrix(0,nrow=hn, ncol=1)

TTT <- sum(ck[,3])
for (kk in 1:hn) {
  super <- (ndata-3)*ck[kk,3]
  sub1 <- 2*(TTT - ck[kk,3])
  F_stat[kk,1] <- super/sub1
}

plot(prd1$freq[1:hn], log(F_stat),
     type = "b", pch = 16, xlab = "Frequency", cex.lab = 1.2)

F_crit <- qf(df1 = 2, df2 = ndata-3, p = 0.95)	## p=1-α(significance level)

abline(h = log(F_crit), col = "red", lwd = 2)
sidx <- which(F_stat > F_crit)
points(prd1$freq[sidx], log(F_stat[sidx]), col = "red", pch = 16, cex = 1.0)
text(0.25, 2.0, "Critical F value", cex=1.2)

alpha <- 0.05
g_fisher <- max(ck[,3])/(sum(ck[,3]))
g_crit <- 1 - (alpha/hn)^(1/(hn-1))
ga_crit <- -2*log(1-(1-alpha)^(1/ndata))
g_fisher
g_crit

## install.packages("ptest")
library(ptest)
ptestg(rr, method = "Fisher")

# convert frequency to period
midx <- which.max(prd1$spec)
Tpeak <- 1/prd1$freq[midx]
Tpeak


