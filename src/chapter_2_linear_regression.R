


idata <- read.csv("../data/Busan_MSL_data_1961_2022.csv")
year <- idata$YEAR
msl <- idata$Mean

model <- lm(msl ~ year)

summary(model)

op <- par(no.readonly = TRUE)
par(mar = c(5,5,5,5))
plot(x = year, y = msl,
     main = "Year~MSL Linear Regression",
     xlab = "YEAR", ylab = "Mean Sea Level(cm)",
     pch = 16, las = 1,
     cex.main = 3, cex.lab = 2.5, cex.axis = 1.5)
abline(model, lwd = 2)


# confidence interval
confint(model)

newx <- seq(min(year), max(year), 0.1)

conf <- predict(model,
                newdata = data.frame(year = newx),
                interval = "confidence")

pred <- predict(model,
                newdata = data.frame(year = newx),
                interval = "prediction")

# 계산한 신뢰구간의 상한선과 하한선
lines(newx, conf[, 3], lty = 2, lwd = 2, col = 2)
lines(newx, conf[, 2], lty = 2, lwd = 2, col = 2)

# 계산한 예측구간의 상한선과 하한선
lines(newx, pred[, 3], lty = 6, lwd = 1.5, col = 4)
lines(newx, pred[, 2], lty = 6, lwd = 1.5, col = 4)

# 범례
legend("topleft", lty = c(1, 2, 6), lwd = c(2, 2, 1),
       seg.len = 3, text.font = 2, box.lty = 0,
       legend = c("reg.line", "conf.interval", "predicted"),
       col = c(1, 2, 4), cex = 1.2,
       horiz = F)
box()



# 모델 요약 정보를 이용한  분포에서의 상·하한 임계치 계산
# alpha = 0.05 조건일 때
alpha <- 0.05
n <- length(summary(model)$residuals)

# 자유도 54인  분포의 유의수준 5% 범위
# 여기서는 상한 97.5%값을 계산하며, 하한 2.5%는 부호가 음수로 바뀜
threshold <- qt(p = 1-alpha/2, df = n-2)

# 임계치 확인
threshold

# 식 5, 6의 계산: 추정치와 표준오차 정의
b_est <- summary(model)[[4]][1] #  추정치
a_est <- summary(model)[[4]][2] #  추정치
b_se <- summary(model)[[4]][3] #  표준오차
a_se <- summary(model)[[4]][4] #  표준오차

# 식 5, 6의 계산: 의 상·하한 오차 계산
b_upper <- b_est + threshold * b_se
b_lower <- b_est - threshold * b_se
a_upper <- a_est + threshold * a_se
a_lower <- a_est - threshold * a_se

# 결과를 Table 형태로 만들기
CI_result <- data.frame(c(b_lower, a_lower), c(b_upper, a_upper))
colnames(CI_result) <- c("2.5 %", "97.5%")
rownames(CI_result) <- c("(Intercept)", "year")

# confint(model) 실행 결과와 값 비교
CI_result

#                    2.5 %        97.5%
# (Intercept)  -467.3804387 -317.2049682
# year           0.1933944    0.2688833
