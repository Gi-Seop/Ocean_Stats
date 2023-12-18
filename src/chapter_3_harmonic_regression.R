

library(lubridate)
library(tidyverse)

# 초기 결측 자료 제거 및 시간 변수의 자료형 변경
df <- read.csv("../data/irregular_timeseries.csv")[-c(1:2),]
df$time <- as.POSIXct(df$time) %>% decimal_date

# 모델 적합에 사용할 자료와 검증용 자료의 분리
df_fit <- df[df$time <= decimal_date(ymd(20101231)),]
df_val <- df[df$time > decimal_date(ymd(20101231)),]



# Harmonic Regression 함수 정의
h_reg <- function(x, y, m){
  
  yy <- as.numeric(y)
  Tj <- 1/m
  wj <- 2*pi/Tj
  
  constant_term <- rep(1, length(x))
  cos_term <- t(apply(as.matrix(x), 1, function(tt) cos(wj*tt)))
  sin_term <- t(apply(as.matrix(x), 1, function(tt) sin(wj*tt)))
  
  if(length(m) == 1){
    cos_term <- t(cos_term)
    sin_term <- t(sin_term)
  }
  
  colnames(cos_term) <- paste0("coef_", "cos", m)
  colnames(sin_term) <- paste0("coef_", "sin", m)
  
  PP <- cbind(constant_term, cos_term, sin_term) 
  CC <- solve(t(PP)%*%PP)%*%(t(PP)%*%yy)
  eyy <- PP%*%CC
  res <- list(coefs = CC[,1],
              fitted = as.numeric(eyy),
              residual = y - as.numeric(eyy))
  return(res)
  
}


# Harmonic Regression 모델의 예측함수
predict_hreg <- function(model, new_x){
  
  new_w <- gsub(pattern = "[a-zA-Z]|_",
                replacement = "",
                x = names(model$coefs[-1])) %>%
    as.numeric %>% unique %>% { 2 * pi * .}
  
  constant_term <- rep(1, length(new_x))
  cos_term <- t(apply(as.matrix(new_x), 1, function(x) cos(new_w*x)))
  sin_term <- t(apply(as.matrix(new_x), 1, function(x) sin(new_w*x)))
  
  if(length(new_w) == 1){
    cos_term <- t(cos_term)
    sin_term <- t(sin_term)
  }
  
  colnames(cos_term) <- paste0("coef_", "cos", new_w)
  colnames(sin_term) <- paste0("coef_", "sin", new_w)
  
  new_PP <- cbind(constant_term, cos_term, sin_term) 
  CC <- model$coefs
  new_y <- as.numeric(new_PP%*%CC)
  
  return(new_y)
  
}



# `h_reg()` 함수를 이용한 수온 주기 성분 계산
## 일 시 단위를 decimal year로 변환(주기 계산 편의를 위함)
deci_year <- df_fit$time
wt <- df_fit$water_temperature

# 사전 정의한 `h_reg()` 함수를 이용하여 1년 주기 성분 계산
wt_p1 <- h_reg(x = deci_year, y = wt, m = c(1, 0.5))

# 계수 확인
wt_p1$coefs


# 초기 그래픽 설정 저장(그래픽 매개변수 복원 위함)
op <- par(no.readonly = T)

# 그래픽 매개변수 설정(글꼴, 여백 등)
par(family = "serif", mar = c(5, 7, 5, 5))

# 관측 자료 시각화
plot(deci_year, wt, type = "p",
     ylim = c(0, 30),
     xlab = "Decimal Years",
     ylab = "Water Temperature",
     main = "Harmonic Regression Fitted",
     pch = 16, cex = 0.5, col = "gray70",
     las = 1, cex.main = 2, cex.lab = 2, cex.axis = 1.7)

# 적합 결과와 잔차 시계열
lines(deci_year, wt_p1$fitted, lwd = 2, col = 2, lty = 1)
lines(deci_year, wt_p1$residual + mean(wt), lwd = 2, col = 4, lty = 1)

# 범례
legend("topright", legend = c("Observed", "Fitted", "Residuals"),
       cex = 1.5,
       col = c("gray70",2,4), pch = c(16, NA, NA),
       lty = c(NA, 1, 1), lwd = 2,
       horiz = T, text.font = 2)

# prediction
new_x <- df_val$time
new_y <- predict_hreg(model = wt_p1, new_x = new_x)

deci_year_val <- df_val$time
wt_val <- df_val$water_temperature

# test data visualization
plot(deci_year_val, wt_val, type = "p",
     ylim = c(0, 30),
     xlab = "Decimal Years",
     ylab = "Water Temperature",
     main = "Harmonic Regression Fitted",
     pch = 16, cex = 1, col = "gray70",
     las = 1, cex.main = 2, cex.lab = 2, cex.axis = 1.7)

# predicted data visualization
lines(new_x, new_y, lwd = 2)
legend("topright", legend = c("Observed", "Fitted"),
       cex = 1.5,
       col = c("gray70",1), pch = c(16, NA),
       lty = c(NA, 1), lwd = 2,
       horiz = T, text.font = 2)

# 상관계수 확인
cor(wt_val, new_y)
# [1] 0.986438

# root mean square error 계산
sqrt(mean((wt_val - new_y)^2))
# [1] 1.204498


# 관측 정보와 예측 정보의 비교
plot(wt_val, new_y,
     main = "Observation vs. Prediction [WT]",
     xlab = "Observation", ylab = "Prediction",
     pch = 16, cex = 2,
     las = 1, cex.main = 2, cex.lab = 2, cex.axis = 1.7)
abline(a = 0, b = 1, lwd = 2)
# 범례
legend("bottomright", legend = "1:1 Line",
       lwd = 2, cex = 2, text.font = 2)
# 그래픽 매개변수 복원
par(op)



# `h_reg()` 함수를 이용한 질산염 주기 성분 계산
## 결측 구간 제거 및 
complete_no3_id <- complete.cases(df_fit$nitrate)
deci_year_no3 <- deci_year[complete_no3_id]
no3 <- df_fit$nitrate[complete_no3_id]

# 사전 정의한 `h_reg()` 함수를 이용하여 다양한 주기의 결합 신호 계산
no3_p1 <- h_reg(x = deci_year_no3, y = no3, m = 1)

# p2, p3 모델은 후보 주기를 다양하게 조합하여 적합시킴
no3_p2 <- h_reg(x = deci_year_no3, y = no3, m = c(1/4, 1/3, 1/2, 1, 2))
no3_p3 <- h_reg(x = deci_year_no3, y = no3, m = c(1:15))

# 계수 확인
round(no3_p1$coefs, 2)

# constant_term     coef_cos1     coef_sin1 
#      8.585438      2.083307      3.338605 

# 관측 자료 시각화
op <- par(no.readonly = T)
par(mar = c(5, 5, 3, 1))

plot(deci_year_no3, no3, type = "p",
     ylim = c(0, 25),
     xlab = "Decimal Years", ylab = "Nitrate",
     main = "Harmonic Regression Fitted",
     pch = 16, cex = 0.5, col = "gray70",
     las = 1, cex.main = 2, cex.lab = 2, cex.axis = 1.7)

# 적합 결과 시각화
lines(deci_year_no3, no3_p1$fitted, lwd = 2, col = 2, lty = 2)
lines(deci_year_no3, no3_p2$fitted, lwd = 2, col = 3, lty = 2)
lines(deci_year_no3, no3_p3$fitted, lwd = 2, col = 4, lty = 2)

# 범례
legend("topright",
       legend = c("Observed", "Fitted: p1", "Fitted: p2", "Fitted: p3"),
       cex = 1.2, col = c("gray70",2,3,4), pch = c(16, NA, NA, NA),
       lty = c(NA, 1, 1, 1), lwd = 2, horiz = T, text.font = 2)

par(op)





# 각 모델 조합의 AIC를 이용한 최적 모델 구조 선정
# 후보 주기 정의 및 frequency grid 생성
guess_m <- c(1/4, 1/3, 1/2, 1, 2, 15)
t_comb <- expand.grid(rep(list(c(F, T)), length(guess_m)),
                      KEEP.OUT.ATTRS = F)[-1, ]
colnames(t_comb) <- guess_m %>% round(3) %>% paste0("t_",.)

# 후보주기 조합에 따른 AIC 계산
model_diag <- apply(
  X = t_comb,
  MARGIN =  1,
  FUN = function(x){
    
    temp_m <- as.numeric(guess_m[t(x)])
    fit <- h_reg(x = deci_year_no3, y = no3, m = temp_m)
    mse <- mean((no3 - fit[[2]])^2, na.rm = T)
    aic <- 2*2*length(temp_m) + length(no3)*log(mse)
    res <- data.frame(aic = aic, mse = mse, t(x))
    return(res)
    
  }
  
) %>%
  do.call("rbind", .) %>%
  arrange(aic)

model_diag[1:10,]


# 선정된 모델의 후보주기(0.5, 1년)를 매개변수로 선정하여 과거 자료로 모델 적합
selected_no3_model <- h_reg(x = deci_year_no3, y = no3, m = c(0.5, 1))

## 검증 기간의 시간을 새 독립변수로 지정
new_x_no3 <- df_val$time

## 적합된 모델에 예측 대상 시간을 입력하여 예측값 계산
new_y_no3 <- predict_hreg(model = selected_no3_model, new_x = new_x_no3)

# 검증자료를 새 이름으로 정의
no3_val <- df_val$nitrate

# 예측 기간 자료 시각화
par(family = "serif", mar = c(5, 7, 5, 5))

# 검증 기간 자료 시계열
plot(new_x_no3, no3_val, type = "p",
     ylim = c(0, 25),
     xlab = "Decimal Years", ylab = "Nitrate",
     main = "Harmonic Regression Fitted",
     pch = 16, cex = 0.5, col = "gray70",
     las = 1, cex.main = 2, cex.lab = 2, cex.axis = 1.7)

# 예측 자료 시계열
lines(new_x_no3, new_y_no3, lwd = 2)

# 상관계수 확인
cor(no3_val, new_y_no3)

# RMSE 계산
sqrt(mean((no3_val - new_y_no3)^2))

# 검증자료와 예측 자료의 비교
plot(no3_val, new_y_no3,
     main = "Observation vs. Prediction [NO3]",
     xlab = "Observation", ylab = "Prediction",
     pch = 16, cex = 2,
     las = 1, cex.main = 2, cex.lab = 2, cex.axis = 1.7)
abline(a = 0, b = 1, lwd = 2)

# 범례
legend("bottomright", legend = "1:1 Line",
       lwd = 2, cex = 2, text.font = 2)


