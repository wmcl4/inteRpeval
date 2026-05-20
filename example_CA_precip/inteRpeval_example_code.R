##EXAMPLE:precipitation in California
library(geodata)
library(gstat)

precip_path <- "/combined_precip_2024.csv"
precip_stations <- read.csv(precip_path)

GPM_path <- "/GPM_jan24.tif"
GPM_jan24 <- rast(GPM_path)

us <- gadm(country = "USA", level = 1, path = tempdir())
ca_boundary <- us[us$NAME_1 == "California", ]

sp_list <- prepare_spatial_data(
  points_data = precip_stations,
  lon_col     = "Longitude",        
  lat_col     = "Latitude",        
  boundary    = ca_boundary,
  dep_var     = "JAN2024",     
  raster_res  = 10000
)
nn_jan24 <- NearNeigh(sp_list, plot = TRUE)
idw_jan24 <- InvDistWt(sp_list, plot = TRUE)
ukrig_jan24 <- UniversalKriging(sp_list, plot = TRUE)

compare_rasters(nn_jan24,GPM_jan24,"NN","GPM",normalize = T)
compare_rasters(idw_jan24,GPM_jan24,"IDW","GPM",normalize = T)
compare_rasters(ukrig_jan24,GPM_jan24,"Krig","GPM",normalize = T)
