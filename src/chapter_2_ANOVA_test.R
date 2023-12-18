
## ANOVA CODE
idata <- read.csv("../data/KOEM_BK1417.csv")
aov1 <- aov(idata[,39] ~ as.factor(idata[,4]), data=idata)
## 25  - surface Total inorganic N, 
summary(aov1)
phoc <- TukeyHSD(aov1, conf.level=0.95)
phoc
plot(phoc)

