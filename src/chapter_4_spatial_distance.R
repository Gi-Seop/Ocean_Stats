
library(globe)
globeearth(eye = list(lon = 128, lat = 36))
flatearth("atlas")

# harversine function
haversine_in_km <- function(lon1, lat1, lon2, lat2){
  
  # 지구 반지름 km
  R <- 6371
  
  dLat <- (lat2 - lat1) * pi / 180
  dLon <- (lon2 - lon1) * pi / 180
  
  a <- sin(dLat/2) * sin(dLat/2) +
    cos(lat1 * pi / 180) * cos(lat2 * pi / 180) *
    sin(dLon/2) * sin(dLon/2)
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  
  return(R * c)
  
}


# 서울시청과 부산시청의 거리
haversine_in_km(126.9783, 37.5667, 129.0747, 35.1794) # [1] 325.0855
# 324.86 Google Earth에서의 거리


library(geosphere)

## Distance calculation between two points
p1 <- c(126.9783, 37.5667)
p2 <- c(129.0747, 35.1794)

distVincentyEllipsoid(p1, p2)
distHaversine(p1, p2)
distVincentySphere(p1, p2)



