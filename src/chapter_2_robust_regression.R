
library(MASS)
idata <- read.csv("../data/Busan_MSL_data_1961_2022.csv")
year <- idata$YEAR
msl <- idata$Mean


op <- par(no.readonly = TRUE)
par(mar = c(5,5,5,5))
plot(x = year, y = msl,
     main = "Year~MSL Robust Linear Regression",
     xlab = "YEAR", ylab = "Mean Sea Level(cm)",
     pch = 16, las = 1,
     cex.main = 2, cex.lab = 2.5, cex.axis = 1.5)

fit1 <- lm(msl ~ year)
fit2 <- rlm(msl ~ year) # 가중계수(범위 0-1) 정보는 fit2$w 
oa <- fit1$coefficients[1]
ob <- fit1$coefficients[2]
ra <- fit2$coefficients[1]
rb <- fit2$coefficients[2]

OLS <- oa + ob*year
RLS <- ra + rb*year


lines(year, OLS, col="red", lwd=2)
lines(year, RLS, col="blue", lwd=2) 


legend("bottomright", pch = c(16, NA, NA),
       lty = c(0, 1, 1), lwd = c(2, 2, 2),
       seg.len = 3, text.font = 2, box.lty = 0,
       legend = c("Data", "Ordinary Curve", "Robust Curve"),
       col = c(1, 2, 4), cex = 1.2,
       horiz = F)

