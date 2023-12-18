

# step 1: loading data and subroutines
source("chapter_4_coordinate_transform.R")
source("chapter_4_spatial_estimation_sub.R")
source("chapter_4_variogram_functions.R")

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


# calculating experimental variogram using function `vge()`
ex_vgm <- vge(xi = sample_2d, s_var = 1:2, val = 3, valid_n = 3,
              max_range = 1, lag_buffer = 0.1, overlap_factor = 2,
              plot = TRUE)


# fitting theoretical variogram using function `vgt()`
vg_fit <- vgt(ex_vgm = ex_vgm,
              vgm_list = list(vgm_s = "Sph"),
              par0 = c(0.5, 0.5, 0.5),
              lower = c(0.1, 0.1, 0.1),
              upper = c(1, 2, 1.5),
              plot = TRUE)

krige_2d_res <- krige(xi = sample_2d, x0 = target_coord,
                      coord = 1:2, val = 3,
                      vgm_model = "Sph",
                      vgm_options = list(vgm_s = "Sph",
                                         pars = vg_fit$par),
                      step_size = 10000,
                      n_cores = 8,
                      progress = T, silence = F)


# target_mean_fields
res_df_kriging <- data.frame(target_grid, krige_2d_res)

str(res_df_kriging)


# plotting chunk
colr <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, 'RdYlBu')))


ext_res <- extent(res_df_kriging[,1:2])
bnd_raster <- raster(ext_res, ncol = 18, nrow = 14) # set resolution as grid size
res_raster <- rasterize(res_df_kriging[, 1:2],
                        bnd_raster,
                        res_df_kriging[, 3],
                        fun = mean)

# estimated value
plot(res_raster)
plot(st_geometry(coastline_utm), col = "gray80", add = T)
points(obs_surface[, c("x", "y")], pch = 16)

# estimated variance
res_raster_var <- rasterize(res_df_kriging[, 1:2],
                            bnd_raster,
                            res_df_kriging[, 4],
                            fun = mean)
plot(res_raster_var)
plot(st_geometry(coastline_utm), col = "gray80", add = T)
points(obs_surface[, c("x", "y")], pch = 16)



# step 4: 10-fold cross validation using function `krige_cv()`

kriging_cv_res <- krige_cv(xi = sample_2d,
                           method = "nfold", n = 10,
                           coord = 1:2, val = 3, vgm_model = "Sph",
                           vgm_options = list(vgm_s = "Sph",
                                              pars = vg_fit$par),
                           step_size = 10000, n_cores = 8,
                           progress = F, silence = T)


str(kriging_cv_res)


op  <- par(no.readonly = T)
par(mar = c(5,5,4,2))
plot(kriging_cv_res$org_val, kriging_cv_res$est_val,
     xlim = c(22, 25), ylim = c(22, 25),
     pch = 16,
     main = "10-fold Cross Validated",
     xlab = "Observed", ylab = "Estimated",
     cex.main = 2, cex.lab = 2, cex.axis = 1.5)
abline(b = 1, a = 0, lwd = 2)
par(op)


sqrt(mean((kriging_cv_res$org_val - kriging_cv_res$est_val)^2))


