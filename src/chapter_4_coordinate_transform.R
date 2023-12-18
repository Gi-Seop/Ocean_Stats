

library(tidyverse)
library(sf)

obs_surface <- read.csv("../data/KOEM_Busan_202208.csv")[-11,]
# except BK1424 (35.04833 129.9883)

prj_str_longlat <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
prj_str_utm_52 <- "+proj=utm +zone=52 +datum=WGS84 +units=m +no_defs"

sf_longlat <- st_as_sf(x = obs_surface,
                       coords = c("long","lat"),
                       crs = prj_str_longlat)
coord_lonlat <- st_coordinates(sf_longlat)


sf_utm <- st_transform(x = sf_longlat, crs = prj_str_utm_52)
coord_utm <- st_coordinates(sf_utm)


obs_surface$x <- coord_utm[, 1]
obs_surface$y <- coord_utm[, 2]


# load coastlines
coastline_lonlat <- read_sf("../data/coastline/coastline_lonlat.shp")
coastline_utm <- read_sf("../data/coastline/coastline_utm.shp")


# check result
op <- par(no.readonly = T)
par(mfrow = c(1,2))

plot(st_geometry(coastline_lonlat), axes = T)
points(obs_surface[, c("long", "lat")], pch = 16)

plot(st_geometry(coastline_utm), axes = T)
points(obs_surface[, c("x", "y")], pch = 16)

par(op)



# grid generation
source("chapter_4_coordinate_transform_sub.R")
target_grid_lonlat <- expand.grid(long = seq(129.01, 129.18, 0.01),
                                  lat = seq(35.03, 35.16, 0.01))

target_grid_utm <- longlat_to_utm52(target_grid_lonlat)

str(target_grid_utm)


# plot(st_geometry(coastline_utm), axes = T)
# points(target_grid_utm, pch = 16, col = 2)
# points(obs_surface[, c("x", "y")], pch = 16)


# land boundary intersection clipping
target_grid_utm_sf <- st_as_sf(x = target_grid_utm,
                               coords = c("X","Y"),
                               crs = prj_str_utm_52)

intersect_points <- st_intersects(
  target_grid_utm_sf,
  coastline_utm
) %>%
  apply(., 1, any)


op <- par(no.readonly = T)
par(mfrow = c(1,2))

plot(st_geometry(coastline_utm), axes = T)
points(target_grid_utm, pch = 16, col = 2)

plot(st_geometry(coastline_utm), axes = T)
points(target_grid_utm[!intersect_points,], pch = 16, col = 2)

par(op)




