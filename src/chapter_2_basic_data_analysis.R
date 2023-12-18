

library(robustbase)

ty_all_data <- read.csv("../data/typhoon_no_data_all.csv")
ty_kor_data <- read.csv("../data/typhoon_no_data.csv")

str(ty_all_data)

year <- ty_all_data$YEAR
nt_all <- ty_all_data$TA
nt_kor <- ty_kor_data$TK
ndata <- length(year)

# basic plot
plot(year, nt_all, type="o", ylim=c(0, 40))
points(year, nt_kor, pch=16, col="red")

# histogram and basic statistics visualization
hist(nt_all, prob=TRUE, breaks="FD", 
	xlab="No. of Typhoons per year",ylab="probability",
	cex.lab=1.3, main="")
box(); grid(lty=3)
mc <- mean(nt_all)
nsd <- sd(nt_all)

mm <- median(nt_all)
dst1 <- density(nt_all)
mp <- dst1$x[which.max(dst1$y)]
lines(dst1$x, dst1$y, col="black", lwd=3)

abline(v=c(mc, mm, mp), col=c("red", "cyan", "magenta"), lwd=3)

nxx <- seq(mc-3*nsd, mc+3*nsd, 0.1)
fxx <- dnorm(nxx, mean=mc, sd=nsd)
lines(nxx, fxx, col="blue", lwd=2)
legend("topright", legend=c("mean", "median", "mode"),
	lty=1, lwd=3, col=c("red", "cyan", "magenta"), cex=1.2)

text(32, 0.08, paste("mean = ", substr(as.character(mc),1,4), sep=""), cex=1.2, adj=0)
text(32, 0.075, paste("S.D. = ", substr(as.character(nsd), 1,3), sep=""), cex=1.2, adj=0)


# boxplot and adjusted boxplot
op <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(3,3,1,1))
boxplot(cbind(nt_all, nt_kor),  horizontal=TRUE, names=c("ALL", "KOREA"))

df1 <- as.data.frame(cbind(nt_all, nt_kor))
adjbox(df1, horizontal=TRUE, names=c("ALL", "KOREA"))
par(op)

# proportion test
prop.test(sum(nt_kor), sum(nt_all))
