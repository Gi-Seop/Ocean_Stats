
library(imputeTS)
library(trend)
library(randtests)
library(nortest)

idata <- read.csv("../data/Busan_MSL_data_1961_2022.csv")

idata <- idata[-1,]
## 1961년 자료는 분석에서 제외

## AMSL, MMSL - annual and monthly MSL data set,

YY <- idata$YEAR
nYY <- length(YY)
month <- seq(1962, 2022, length.out=61*12)

AMSL <- idata$Mean
## 행렬형식 자료를 벡터로 변환하는 경우, 행 연산이 우선, 열 연산을 위해서는 전치행렬 이용
MMSL <- as.vector(t(idata[,2:13]))

## Missing data filling-in
AMSL <- na_kalman(AMSL)
MMSL <- na_kalman(MMSL)

plot(YY, AMSL, type="o")
plot(month,MMSL, type="l", lwd=1.2)


# fittin model and residual analysis
fit1 <- lm(AMSL ~ YY)
r_AMSL <- fit1$residuals 

plot(YY, r_AMSL, type="o")


mk.test(AMSL)
mk.test(r_AMSL)

runs.test(AMSL)
runs.test(r_AMSL)


bartels.rank.test(AMSL)
bartels.rank.test(r_AMSL)

turning.point.test(AMSL)
turning.point.test(r_AMSL)


### Decomposition of the time series data
tAMSL <- ts(AMSL, frequency = 1, start = 1962)
tMMSL <- ts(MMSL,frequency = 12, start = c(1962,1))

dmp1 <- stl(tMMSL, s.window="periodic")  ## {stats} 기본 패키지

plot(dmp1)  
## order - seasonal, trend, remainder

var(dmp1$time.series[,1])
var(dmp1$time.series[,2])
var(dmp1$time.series[,3])

acf(r_AMSL)

## 성분분리 잔차(reaminder)의 분포함수 - 정규분포 적합
rr <- as.numeric(dmp1$time.series[,3])
hist(rr, prob=T)
mr <- mean(rr); msd <- sd(rr); 
mxx <- seq(min(rr), max(rr), 0.1)
fxx <- dnorm(mxx, mean=mr, sd=msd) 
lines(mxx, fxx, col="blue", lwd=2)

shapiro.test(rr)


 