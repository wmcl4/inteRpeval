library(rhdf5)
library(terra)
library(geodata)

gpm_to_raster <- function(file, hours_in_month, boundary) {
  
  p <- h5read(file, "Grids/G2/precipTotRate/mean")
  p[p < 0] <- NA
  p2d <- apply(p, c(1, 2), mean, na.rm = TRUE)
  p2d <- p2d * hours_in_month
  p2d <- t(p2d)
  
  r <- rast(p2d)
  ext(r) <- c(-180, 180, -67, 67)
  crs(r) <- "EPSG:4326"
  r <- flip(r, direction = "vertical")
  
  boundary <- project(boundary, crs(r))
  r_ca <- crop(r, boundary, mask = TRUE)
  
  return(r_ca)
}
us <- gadm(country = "USA", level = 1, path = tempdir())
ca_boundary <- us[us$NAME_1 == "California", ]
# January (744 hours)
GPM_jan24 <- gpm_to_raster(
  file = "C:/Users/mclou/Downloads/3B-MO.GPM.DPRGMI.CORRAGM.20240101-S000000-E235959.01.V07C.HDF5",
  hours_in_month = 744,
  boundary = ca_boundary
)
plot(GPM_jan24, main = "GPM Total Precip - California Jan 2024 (mm)")

writeRaster(GPM_jan24,  "GPM_jan24.tif")
