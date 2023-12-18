

# step 1: loading data and subroutines
source("chapter_4_coordinate_transform.R")
source("chapter_4_spatial_estimation_sub.R")

sample_2d <- obs_surface[, c("x", "y", "water_temperature")]

# pre defined object
## number of samples
coord_var_names <- c("x", "y")

# scaling observation data.frame coordinates
## for revert scaling values
min_x <- min(sample_2d$x)
max_x <- max(sample_2d$x)
min_y <- min(sample_2d$y)
max_y <- max(sample_2d$y)
range_x <- max(sample_2d$x) - min(sample_2d$x)
range_y <- max(sample_2d$y) - min(sample_2d$y)
sample_2d[, coord_var_names] <- apply(
  sample_2d[, coord_var_names],
  2,
  scale_minmax
)


# scaling target grid coordinates
target_grid <- target_grid_utm[!intersect_points,]
colnames(target_grid) <- c("x", "y")
target_coord <- data.frame(x = (target_grid$x - min_x)/(range_x),
                           y = (target_grid$y - min_y)/(range_y))


idw_res <- idw(xi = sample_2d, x0 = target_coord,
               coord = 1:2, val = 3, order = 2,
               step_size = 10000, n_cores = 8,
               progress = T, silence = F)


# target_mean_fields
res_df_idw <- data.frame(target_grid, idw_res)



# plotting chunk
colr <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlBu')))


ext_res <- extent(res_df_idw[,1:2])
bnd_raster <- raster(ext_res, ncol = 18, nrow = 14) # set resolution as grid size
res_raster <- rasterize(res_df_idw[, 1:2],
                        bnd_raster,
                        res_df_idw[, 3],
                        fun = mean)

plot(res_raster)
plot(st_geometry(coastline_utm), col = "gray80", add = T)
points(obs_surface[, c("x", "y")], pch = 16)



# step 4: 10-fold cross validation using function `idw_cv()`

idw_cv_res <- idw_cv(xi = sample_2d,
                     method = "nfold", n = 10,
                     coord = 1:2, val = 3, order = 2,
                     step_size = 10000, n_cores = 8,
                     progress = F, silence = T)

op  <- par(no.readonly = T)
par(mar = c(5,5,4,2))
plot(idw_cv_res$org_val, idw_cv_res$est_val,
     xlim = c(22, 25), ylim = c(22, 25),
     pch = 16,
     main = "10-fold Cross Validated",
     xlab = "Observed", ylab = "Estimated",
     cex.main = 2, cex.lab = 2, cex.axis = 1.5)
abline(b = 1, a = 0, lwd = 2)
par(op)


sqrt(mean((idw_cv_res$org_val - idw_cv_res$est_val)^2))



