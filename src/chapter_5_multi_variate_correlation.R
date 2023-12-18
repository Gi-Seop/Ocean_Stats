
library(corrplot)

idata <- read.csv("../data/KOEM_BK1417.csv")

str(idata)

vdata <- as.data.frame(idata[, 8:18])
cstr <- c("SD", "WT.S", "WT.B", "S.S", "S.B", "pH.S", "pH.B", 
		"DO.S", "DO.B", "COD.S", "COD.B")

colnames(vdata) <- cstr

cmat <- cor(vdata, method="spearman")
corrplot(cmat)
corrplot.mixed(cmat, lower="circle", upper="number")

