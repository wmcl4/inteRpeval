# inteRpeval

An R package for basic spatial interpolation and raster comparison, built around the `terra` and `gstat` libraries.

---

## Functions

### `prepare_spatial_data()`
Transforms data frames to a readable format which is used in all other functions. Builds a Voronoi polygon raster of the response variable and returns all spatial objects needed by the interpolation functions as a single named list.

### `RMSE()`
Computes Root Mean Square Error between observed and predicted values. Used internally.

### `NearNeigh()`
Nearest Neighbour interpolation via `gstat`. Assigns each prediction location the value of its closest observed point. Evaluates performance using 5-fold cross-validation and reports relative performance against a null (mean) model.

### `InvDistWt()`
Inverse Distance Weighting interpolation. Predictions are a distance-weighted average of the nearest 12 observations, with an inverse distance power of 2. Also evaluated with 5-fold cross-validation.

### `UniversalKriging()`
Universal Kriging interpolation that accounts for large-scale spatial trends in x and y. Automatically fits a spherical variogram model to the trend residuals.

### `compare_rasters()`
A  tool for comparing any two `SpatRaster`s. Reprojects and resamples both rasters to a common grid, then computes and plots their pixel-wise difference on a Red-Blue color scale. Prints mean difference, MAE, RMSE, and min/max statistics to the console. Optional min-max normalisation for comparing methods with different value ranges.

---

## Basic Workflow

```r
library(inteRpeval)
library(terra)

# 1. Prepare spatial data
sp_list <- prepare_spatial_data(
  points_data = my_data,
  lon_col     = "longitude",
  lat_col     = "latitude",
  boundary    = ca_boundary,   # SpatVector polygon
  dep_var     = "precipitation"
)

# 2. Run interpolations
nn_result  <- NearNeigh(sp_list, count = 5, plot = TRUE)
idw_result <- InvDistWt(sp_list, plot = TRUE)
uk_result  <- UniversalKriging(sp_list, count = 12, plot = TRUE)

# 3. Compare two methods
compare_rasters(
  r1      = idw_result,
  r2      = uk_result,
  r1_name = "IDW",
  r2_name = "Universal Kriging",
  normalize = TRUE
)
```

---

## Dependencies

- [`terra`](https://cran.r-project.org/package=terra) — spatial data handling, rasterization, and Voronoi construction
- [`gstat`](https://cran.r-project.org/package=gstat) — geostatistical modelling and interpolation

---

## Acknowledgements

Special thanks to **Robert Hijmans** and the [**rspat**](https://rspatial.org) project for their spatial interpolation tutorials and the code foundations for this package.
---

## Author

William McLoughlin
