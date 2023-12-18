
library(forecast)

idata <- read.csv("../data/co2_monthly_data.csv")

## Four station - 안면도, 제주 고산, 울릉도, 독도 
str(idata)
ndata <- nrow(idata)
nyear <- ndata/12

YY <- seq(1999,2022,1) 
month <- seq(1999, 2022, length.out=nyear*12)

ACO2 <- idata[,2] 
plot(month, ACO2, type="l")

tACO2 <- ts(ACO2, frequency=12, start=c(1999,1))

dmp1 <- stl(tACO2, s.window="periodic")  ## {stats} 기본 패키지
plot(dmp1)

var(dmp1$time.series[,1])
var(dmp1$time.series[,2])
var(dmp1$time.series[,3])

aarima <- auto.arima(ACO2)
asarima <- auto.arima(tACO2)

## Performance test using 2021 - 2022 data. 
## Model setup using 1999-2020 data

ndata <- length(ACO2)
test_data <-  ACO2[(ndata-24+1):ndata]
model_data <- ts(ACO2[1:(ndata-24)], frequency=12, start=c(1999,1))
asarima <- auto.arima(model_data)

fct1 <- forecast(asarima, 24)
fest <- fct1$mean
plot(test_data, fest, xlab="Observed CO2", ylab="Estimated CO2", type="p")
abline(0,1, col="red", lwd=2)

rmse <- sqrt(sum((test_data-fest)^2)/length(test_data))
estr <- substr(as.character(rmse),1,5)
legend("topleft", paste("RMSE = ",estr,sep=""), cex=1.3)


library(randtests)
library(nortest)
library(trend)

res_diag <- lm(fest ~ test_data)

turning.point.test(res_diag$residuals)
mk.test(res_diag$residuals)
ad.test(res_diag$residuals)



