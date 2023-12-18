
idata <- read.csv("../data/KOEM_BK1417.csv")

str(idata)

vdata <- as.data.frame(idata[, 8:18])
cstr <- c("SD", "WT.S", "WT.B", "S.S", "S.B", "pH.S", "pH.B", 
		"DO.S", "DO.B", "COD.S", "COD.B")

colnames(vdata) <- cstr


##-------------------------------------------------
## H-clustering 
dd <- dist(t(vdata))
hc <- hclust(dd, method = "complete", members = NULL)
plot(hc)

## install.packages("factoextra")
library(factoextra)

kmc <- kmeans(vdata, 4)
fviz_cluster(kmc, data = vdata, palette = "jco")


## Dsicriminant analysis... 
library(MASS)
cdata <- as.data.frame(idata[, c(8:18, 38, 39)])
cstr <- c("SD", "WT.S", "WT.B", "S.S", "S.B", "pH.S", "pH.B", 
		"DO.S", "DO.B", "COD.S", "COD.B", "chlorophyll-a", "WQI")

colnames(cdata) <- cstr

model <- lda(as.factor(WQI) ~ ., cdata)
plot(model, cex.lab=1.3, cex=1.3, col="blue")
