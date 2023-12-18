
idata <- read.csv("../data/KOEM_BK1417.csv")

str(idata)

vdata <- as.data.frame(idata[, 8:18])

## Multiple linear regression --- Secchi depth (종속변수) 
df1 <- vdata[,-1]
TR <- vdata[,1]
mfit <- lm(TR ~ ., data = df1) # 시간 열을 제외한 모든 변수를 사용
summary(mfit)


## Another try!!! chlorophyll-a 항목 예측으로 변경

cdata <- as.data.frame(idata[, c(8:18, 38)])
cstr <- c("SD", "WT.S", "WT.B", "S.S", "S.B", "pH.S", "pH.B", 
		"DO.S", "DO.B", "COD.S", "COD.B", "chlorophyll-a")
dim(cdata)
TR1 <- as.numeric(cdata[,12])
df2 <- cdata[,-12]
mfit <- lm(TR1 ~ ., data = df2) # 시간 열을 제외한 모든 변수를 사용
summary(mfit)

plot(mfit)

#install.packages("fmsb")
library(fmsb)
VIF(mfit)

