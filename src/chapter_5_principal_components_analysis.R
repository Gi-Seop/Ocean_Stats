

idata <- read.csv("../data/KOEM_BK1417.csv")

str(idata)

vdata <- as.data.frame(idata[, 8:18])
cstr <- c("SD", "WT.S", "WT.B", "S.S", "S.B", "pH.S", "pH.B", 
		"DO.S", "DO.B", "COD.S", "COD.B")

colnames(vdata) <- cstr
nvar <- ncol(vdata)

svdata <- scale(vdata)
dc1 <- svd(svdata)
sqrt(dc1$d)
cvm <- cov(svdata)

ec1 <- eigen(cvm)
ec1

## Orthogonality check!!! 
p1 <- 2; p2 <- 3
sum(ec1$vectors[,p1] * ec1$vectors[,p2])

round(cumsum(ec1$values)/sum(ec1$values), 2)
plot(1:nvar, ec1$values)


pca1 <- princomp(svdata)

cpp <- cumsum(prop.table(pca1$sdev^2))
xscr <- 1:nvar

plot(xscr, cpp, type = "o")
abline(h=c(0.7, 1.0), col="blue", lwd=2)

summary(pca1)


# pc1, pc2
biplot(pca1)
grid(lty = 3, col = "red")

# pc1, pc3
biplot(pca1, choices = c(1,3))  ## PC 축 번호 선택
grid(lty = 3, col = "red")





