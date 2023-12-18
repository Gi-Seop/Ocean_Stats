
# require(sf)

longlat_to_utm52 <- function(lonlat){
  
  colnames(lonlat) <- c("x", "y")
  
  temp_lonlat <- st_as_sf(
    as.data.frame(lonlat),
    coords = c("x", "y"),
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
  )
  
  temp_utm <- st_transform(
    temp_lonlat,
    crs = "+proj=utm +zone=52 +datum=WGS84 +units=m +no_defs"
  )
  
  return(as.data.frame(st_coordinates(temp_utm)))
  
}


utm52_to_longlat <- function(utm52){
  
  colnames(utm52) <- c("x", "y")
  
  temp_utm <- st_as_sf(
    as.data.frame(utm52),
    coords = c("x", "y"),
    crs = "+proj=utm +zone=52 +datum=WGS84 +units=m +no_defs"
  )
  
  temp_lonlat <- st_transform(
    temp_utm,
    crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
  )
  
  return(as.data.frame(st_coordinates(temp_lonlat)))
  
}
