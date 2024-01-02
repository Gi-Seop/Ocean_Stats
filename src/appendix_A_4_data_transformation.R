

## BK1417 지점 표층 수온자료 할당
SWT <- read.csv("../data/KOEM_BK1417.csv")[,9] 
head(SWT) ## 자료의 기본 구조 (설명변수가 변수이름으로 할당)

plot(SWT, xlab="No. of data index", ylab="WT (C)", type="o")
# 그림 생략

qqnorm(SWT)
qqline(SWT, col="red", lwd=2)


library(nortest)
ad.test(SWT)



# install.packages(EnvStats)		## 한번 설치한 경우, 재설치는 필요 없음.
library(EnvStats) 
trn1 <- boxcox(SWT, optimize=TRUE)
trn1$lambda				## 최적 매개변수(lambda) 추정 결과 = -0.10 
# [1] -0.1034598
TBCxx <- boxcoxTransform(SWT, trn1$lambda)			## Box-cox 변환 자료
qqnorm(TBCxx)
qqline(TBCxx, col="red", lwd=2)


ad.test(TBCxx) ## 또는 shapiro.test()


nn <- length(SWT)
## Cumulative distribution function value (probability) assignment, (m-0.5)/n
cf_lower <- (0:(nn-1))/nn
cf_upper <- (1:nn)/nn
cf_center <- (cf_lower+cf_upper)/2

stxx <- sort(SWT)
TLxx <- qnorm(cf_lower)
TUxx <- qnorm(cf_upper)
TMxx <- qnorm(cf_center)		## 완벽한 변환 자료, TMxx

qqnorm(TMxx)
qqline(TMxx, col="red", lwd=2)

plot(stxx, TMxx, main="Point-by-point Transformation", type="o")
