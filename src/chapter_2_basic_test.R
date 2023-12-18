

# 필요한 라이브러리 호출
# install.packages(c("randtests", "trend")) # 설치가 되어있지 않은 경우 실행
library(randtests)
library(trend)

runs.test(nt_all)
runs.test(nt_kor)
turning.point.test(nt_all)
mk.test(nt_all)


# 필요한 라이브러리 호출
# install.packages("nortest") # 설치가 되어있지 않은 경우 실행
library(nortest)

# 자료 불러오기
secchi_depth <- read.csv("../data/KOEM_BK1417.csv")[,8]

shapiro.test(secchi_depth)

ad.test(secchi_depth)


qqnorm(secchi_depth)
qqline(secchi_depth)


# transform
# 필요한 라이브러리 호출
# install.packages("EnvStats") # 설치가 되어있지 않은 경우 실행
library(EnvStats)

bct <- boxcox(secchi_depth, optimize=TRUE)
bct$lambda


# 최적 추정 매개변수를 이용한 Secchi depth 자료변환
lambda_opt <- bct$lambda
Tdata <- boxcoxTransform(secchi_depth, lambda=lambda_opt)
plot(Tdata)


shapiro.test(Tdata)

ad.test(Tdata)


