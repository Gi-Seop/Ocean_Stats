

# source("chapter_5_principal_components_analysis.R")

fa1 <- factanal(covmat=cov(svdata), factors=6, n.obs=92)
str(fa1)
fa1$PVAL
fa1$loadings
fa1$uniqueness
fa1$correlation

