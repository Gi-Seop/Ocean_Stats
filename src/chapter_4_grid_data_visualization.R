

library(raster)

etopo_raster <- raster("../data/ETOPO_kr.tiff")
etopo_mat <- as.matrix(etopo_raster)

class(etopo_raster)
class(etopo_mat)


plot(etopo_raster)


image(etopo_mat)

# image rotation 90 degrees to clockwise
etopo_mat_t <- t(apply(etopo_mat, 2, rev))
image(etopo_mat_t)



library(reshape2)
etopo_melt <- melt(etopo_mat_t)
names(etopo_melt) <- c("x", "y", "z")

library(ggplot2)
# visualize raster object with ggplot2
ggplot(etopo_melt, aes(x, y, fill = z)) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  theme(panel.grid = element_blank())



# 3D data visualization using persp
persp(etopo_mat_t, ticktype = "detailed",
      phi = 40, theta = 40,
      xlab = "x", ylab = "y", zlab = "z")

# 3D data vidualization using plot3D
library(plot3D)
persp3D(z = etopo_mat_t, ticktype = "detailed",
        xlab = "x", ylab = "y", zlab = "z")


# 3D data visualization using rayshader

library(rayshader)
library(ggplot2)
library(viridis)
library(tidyverse)

#Data from Social Security administration
etopo_3d <- ggplot(etopo_melt, aes(x, y, fill = z)) +
  geom_raster() +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  theme(panel.grid = element_blank())

plot_gg(etopo_3d, multicore = TRUE, height = 5, width = 6, scale = 300)


# raster and contour lines using base plot
plot(etopo_raster)
contour(etopo_raster, add = TRUE, nlevels = 20, col = "gray50")
contour(etopo_raster, add = TRUE, nlevels = 10, col = "white")

# add contour lines
ggplot(etopo_melt, aes(x, y, fill = z)) +
  geom_raster() +
  geom_contour(aes(z = z), colour = "white") +
  coord_equal() +
  scale_fill_gradientn(colours = terrain.colors(10)) +
  theme_bw() +
  theme(panel.grid = element_blank())


