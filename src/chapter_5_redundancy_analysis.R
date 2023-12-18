
# libraries
library(tidyverse)
library(vegan)


# data path
env_path <- "../data/env_data.csv"
sp_path <- "../data/spe_data.csv"


# data load
xx <- read.csv(env_path, row.names = 1)[1:11,] %>% as.matrix
yy <- read.csv(sp_path, row.names = 1)[1:11,] %>% as.matrix

# number of observations and centered-data transformation 
# in this note, using hellinger transformation(Legendre & Gallagher, 2001)
nn <- nrow(xx)
xs <- scale(xx)
ys <- decostand(yy, 'hellinger')


# using library 'vegan'
op <- par(no.readonly = T)
par(mfrow=c(2,1), mar=c(3,3,1,1))

rda_sp <- vegan::rda(ys ~ ., as.data.frame(xs), scale = T)
summary (rda_sp)

plot(rda_sp, cex.lab=1.2, cex.axis=1.8)

cca_sp <- vegan::cca(ys ~ ., as.data.frame(xs), scale = T)
summary(cca_sp)
plot(cca_sp,  cex.lab=1.2, cex.axis=1.8)

par(op)

# check linear assumption using DCA if Axis length is > < ?
decorana(ys)




