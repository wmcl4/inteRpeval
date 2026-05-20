#' Prepare Spatial Data and Voronoi Raster
#'
#' Cleans and reprojects point data, builds a Voronoi polygon raster, and
#' returns all spatial objects needed for interpolation functions.
#'
#' @param points_data A \code{data.frame} containing coordinates and attributes.
#' @param lon_col Character. Name of the longitude column in \code{points_data}.
#' @param lat_col Character. Name of the latitude column in \code{points_data}.
#' @param boundary A \code{SpatVector} polygon defining the study area boundary.
#' @param dep_var Character. Name of the target/response variable column.
#' @param input_crs Character. CRS of the input points. Defaults to NAD83
#'   geographic (\code{"+proj=longlat +datum=NAD83"}).
#' @param target_crs Character. CRS to project data into. Defaults to California
#'   Albers Equal Area.
#' @param raster_res Numeric. Raster resolution in map units (metres). Default
#'   is \code{10000}.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{dta}{Projected \code{SpatVector} of input points.}
#'   \item{cata}{Projected \code{SpatVector} of the boundary polygon.}
#'   \item{r}{Empty template \code{SpatRaster} at the specified resolution.}
#'   \item{vr}{Voronoi-rasterized \code{SpatRaster} of the response variable.}
#'   \item{df}{A \code{data.frame} of projected coordinates and attributes for use with \code{gstat}.}
#'   \item{dep_var}{Sanitised name of the response variable (via \code{make.names}).}
#' }
#'
#' @importFrom terra vect project voronoi rast rasterize values mask crop
#' @export
prepare_spatial_data <- function(
    points_data,          
    lon_col,              
    lat_col,              
    boundary,             
    dep_var,              
    input_crs  = "+proj=longlat +datum=NAD83", #default North America Datum 
    target_crs = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120
                  +x_0=0 +y_0=-4000000 +datum=WGS84 +units=m", #default California
    raster_res = 10000
) {
  #Preserve column names 
  safe_var <- make.names(dep_var)
  names(points_data)[names(points_data) == dep_var] <- safe_var
  dep_var <- safe_var
  
  #(re)project points 
  dsp <- vect(points_data, c(lon_col, lat_col), crs = input_crs)
  dta <- project(dsp, target_crs)
  cata <- project(boundary, target_crs)
  
  
  #Voronoi raster mask
  v <- voronoi(dta)
  values(v) <- values(dta)
  vca <- crop(v, cata)
  r <- rast(cata, res = raster_res)
  vr <- rasterize(vca, r, field = dep_var) 
  
  #Data frame with projected coords (instead of SpatVector; for gstat)
  df <- data.frame(
    geom(dta)[, c("x", "y")],
    as.data.frame(dta)
  )
  #Return as a list
  list(
    dta      = dta,
    cata     = cata,
    r        = r,
    vr       = vr,
    df       = df,
    dep_var  = dep_var
  )
}

#' Root Mean Square Error
#'
#' Calculates the Root Mean Square Error (RMSE) between observed and predicted
#' values, ignoring \code{NA}s.
#'
#' @param observed Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values (same length as
#'   \code{observed}).
#'
#' @return A single numeric value representing the RMSE.
#'
#' @examples
#' RMSE(c(1, 2, 3), c(1.1, 1.9, 3.2))
#'
#' @export
RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm = TRUE))
}

#' Nearest Neighbor Interpolation
#'
#' Performs Nearest Neighbor spatial interpolation using \code{gstat}, with
#' 5-fold cross-validation to evaluate predictive performance.
#'
#' @param spatial_list A named list produced by \code{\link{prepare_spatial_data}}.
#' @param count Integer. Number of nearest neighbors to use. Default is
#'   \code{5}.
#' @param plot Logical. If \code{TRUE}, plots the interpolated raster. Default
#'   is \code{FALSE}.
#'
#' @return Invisibly returns a masked \code{SpatRaster} of Nearest Neighbor
#'   predictions. Also prints the relative cross-validation performance to the
#'   console (1 - RMSE / null RMSE; higher is better).
#' @importFrom gstat gstat gstat
#' @importFrom terra interpolate mask
#' @export
NearNeigh <- function(spatial_list, count = 5, plot = FALSE) {
  df      <- spatial_list$df
  r       <- spatial_list$r
  vr      <- spatial_list$vr
  dep_var <- spatial_list$dep_var
  
  formula_obj <- as.formula(paste(dep_var, "~ 1"))
  dep_values <- df[[dep_var]]
  
  gs <- gstat(formula = formula_obj, locations = ~x+y,
                 data = df, nmax = count, set = list(idp = 0))
  nn <- interpolate(r, gs, debug.level = 0)[[1]]   # keep var1.pred only
  nmask <- mask(nn, vr)
  
  # Cross-validation
  null <- RMSE(mean(dep_values), dep_values)
  kf <- sample(1:5, nrow(df), replace = TRUE)
  rmsenn <- rep(NA, 5)
  
  for (k in 1:5) {
    test <- df[kf == k, ]
    train <- df[kf != k, ]
    gscv <- gstat(formula = formula_obj, locations = ~x+y,
                   data = train, nmax = count, set = list(idp = 0))
    p <- predict(gscv, test, debug.level = 0)$var1.pred
    rmsenn[k] <- RMSE(test[[dep_var]], p)
  }
  
  cat("Nearest Neighbor - Relative performance:", 1 - mean(rmsenn) / null, "\n")
  if (plot) terra::plot(nmask, main = paste("Nearest Neighbor Prediction:", dep_var))
  
  invisible(nmask)
}

#' Inverse Distance Weighting Interpolation
#'
#' Performs Inverse Distance Weighting (IDW) spatial interpolation using
#' \code{gstat} with a local neighbourhood of 12 stations and an inverse
#' distance power of 2, with 5-fold cross-validation to evaluate performance.
#'
#' @param spatial_list A named list produced by \code{\link{prepare_spatial_data}}.
#' @param plot Logical. If \code{TRUE}, plots the interpolated raster. Default
#'   is \code{FALSE}.
#'
#' @return Invisibly returns a masked \code{SpatRaster} of IDW predictions.
#'   Also prints the relative cross-validation performance to the console
#'   (1 - RMSE / null RMSE; higher is better).
#' @importFrom gstat gstat gstat
#' @importFrom terra interpolate mask
#' @export
InvDistWt <- function(spatial_list, plot = FALSE) {
  df      <- spatial_list$df
  r       <- spatial_list$r
  vr      <- spatial_list$vr
  dep_var <- spatial_list$dep_var
  
  formula_obj <- as.formula(paste(dep_var, "~ 1"))
  dep_values <- df[[dep_var]]
  
  # Restricted to the nearest 12 stations for forced "local"
  gs <- gstat(formula = formula_obj, locations = ~x+y, data = df, 
              nmax = 12, set = list(idp = 2)) #current idp weighting = 2
  idw <- interpolate(r, gs, debug.level = 0)[[1]]
  idwr <- mask(idw, vr)
  
  # Cross-validation
  null <- RMSE(mean(dep_values), dep_values)
  kf <- sample(1:5, nrow(df), replace = TRUE)
  rmsenn <- rep(NA, 5)
  
  for (k in 1:5) {
    test <- df[kf == k, ]
    train <- df[kf != k, ]
    gscv <- gstat(formula = formula_obj, locations = ~x+y, data = train,
                   nmax = 12, set = list(idp = 2))
    p <- predict(gscv, test, debug.level = 0)$var1.pred
    rmsenn[k] <- RMSE(test[[dep_var]], p)
  }
  
  cat("IDW - Relative performance:", 1 - mean(rmsenn) / null, "\n")
  if (plot) terra::plot(idwr, 1, main = paste("IDW Prediction:", dep_var))  
  invisible(idwr)
}

#' Universal Kriging Interpolation
#'
#' Performs Universal Kriging spatial interpolation using \code{gstat}, accounting
#' for large-scale spatial trends (x and y as covariates). Automatically fits a
#' spherical variogram model to the trend residuals and evaluates predictive 
#' performance via 5-fold cross-validation.
#' 
#' @param spatial_list A named list produced by \code{\link{prepare_spatial_data}}.
#' @param count Integer. Maximum number of nearest neighbours to use in the local
#'   Kriging neighborhood. Default is \code{12}.
#' @param plot Logical. If \code{TRUE}, plots the interpolated raster. Default
#'   is \code{FALSE}.
#'
#' @details
#' The trend formula is defined as \code{dep_var ~ x + y}, meaning the mean
#' structure is modeled as a linear function of the projected coordinates.
#' Kriging is then performed on the residuals of this trend.
#'
#' Variogram initial parameters are estimated directly from the empirical
#' variogram:
#' \itemize{
#'   \item \strong{Nugget:} minimum semivariance value.
#'   \item \strong{Sill:} range of semivariance (max minus min).
#'   \item \strong{Range:} one-third of the maximum lag distance, a stable
#'     starting estimate.
#' }
#' These smart initials are recalculated independently for each cross-validation
#' fold to reflect the specific training subset.
#'
#' @return Invisibly returns a masked \code{SpatRaster} of Universal Kriging
#'   predictions. Also prints the relative cross-validation performance to the
#'   console (1 - RMSE / null RMSE; higher is better).
#' @importFrom gstat gstat variogram fit.variogram vgm
#' @importFrom terra interpolate mask
#' @export
UniversalKriging <- function(spatial_list, count = 12, plot = FALSE) {
  df      <- spatial_list$df
  r       <- spatial_list$r
  vr      <- spatial_list$vr
  dep_var <- spatial_list$dep_var
  
  # Define universal formula accounting for spatial trends (x and y)
  formula_obj <- as.formula(paste(dep_var, "~ x + y"))
  dep_values <- df[[dep_var]]
  
  # Fit a sample variogram on the RESIDUALS of the spatial trend
  v_mod <- variogram(formula_obj, locations = ~x+y, data = df)
  
  # SMART INITIALS: Extract directly from the empirical variogram data
  init_nugget <- min(v_mod$gamma)
  init_sill <- max(v_mod$gamma) - init_nugget
  init_range <- max(v_mod$dist) / 3   # 1/3 of the max distance is a stable starting range
  
  f_mod <- fit.variogram(v_mod, vgm(psill = init_sill, model = "Sph", 
                                    range = init_range, nugget = init_nugget))
  
  # Build gstat object using the trend model and local neighborhood restriction
  gs <- gstat(formula = formula_obj, locations = ~x+y, data = df, 
                model = f_mod, nmax = count)
  uk  <- interpolate(r, gs, debug.level = 0)[[1]]
  ukr <- mask(uk, vr)
  
  # Cross-validation
  null <- RMSE(mean(dep_values), dep_values)
  kf <- sample(1:5, nrow(df), replace = TRUE)
  rmseuk <- rep(NA, 5)
  
  for (k in 1:5) {
    test <- df[kf == k, ]
    train <- df[kf != k, ]
    
    # Fit variogram for this specific fold's training residuals
    v_mod_cv <- variogram(formula_obj, locations = ~x+y, data = train)
    
    # SMART INITIALS FOR CV: Dynamically adjust to this loop's specific subset
    train_nugget <- min(v_mod_cv$gamma)
    train_sill <- max(v_mod_cv$gamma) - train_nugget
    train_range <- max(v_mod_cv$dist) / 3
    
    f_mod_cv <- fit.variogram(v_mod_cv, vgm(psill = train_sill, model = "Sph", 
                                            range = train_range, nugget = train_nugget))
    
    gscv <- gstat(formula = formula_obj, locations = ~x+y, data = train, 
                        model = f_mod_cv, nmax = count)
    p <- predict(gscv, test, debug.level = 0)$var1.pred
    rmseuk[k] <- RMSE(test[[dep_var]], p)
  }
  
  cat("Universal Kriging - Relative performance:", 1 - mean(rmseuk) / null, "\n")
  if (plot) terra::plot(ukr, 1, main = paste("Universal Kriging Prediction:", dep_var))  
  invisible(ukr)
}


#' Compare Two Spatial Rasters
#'
#' Reprojects and resamples two \code{SpatRaster} objects onto a common grid,
#' then computes and plots their pixel-wise difference, with summary error
#' statistics printed to the console.
#'
#' @param r1 A \code{SpatRaster}. The reference raster.
#' @param r2 A \code{SpatRaster}. The raster to subtract. Will be
#'   reprojected and resampled to match \code{r1}.
#' @param r1_name Character. Display name for \code{r1} used in plot titles and
#'   console output. Default is \code{"Raster 1"}.
#' @param r2_name Character. Display name for \code{r2} used in plot titles and
#'   console output. Default is \code{"Raster 2"}.
#' @param normalize Logical. If \code{TRUE}, both rasters are min-max normalised
#'   to \code{[0, 1]} before differencing, making methods with different value
#'   ranges directly comparable. Default is \code{FALSE}.
#'
#' @details
#' If either raster has more than one layer (e.g. raw \code{gstat} output), only
#' the first layer is used. \code{r2} is always aligned to \code{r1}'s CRS and
#' resolution via bilinear resampling before any computation.
#'
#' The following statistics are printed for the difference raster
#' (\code{r1 - r2}):
#' \itemize{
#'   \item \strong{Mean difference} — systematic bias between the two surfaces.
#'   \item \strong{Mean Absolute Error (MAE)} — average magnitude of disagreement.
#'   \item \strong{RMSE} — penalizes large local discrepancies more heavily.
#'   \item \strong{Min / Max} — the extremes of the difference surface.
#' }
#'
#' @return Invisibly returns the difference \code{SpatRaster} (\code{r1 - r2}).
#'   A diverging Red-Blue plot of the difference is produced as a side effect.
#' @importFrom terra nlyr project crs resample minmax values
#' @export
compare_rasters <- function(r1, r2, r1_name = "Raster 1", r2_name = "Raster 2", normalize = FALSE) {
  
  # Take prediction layer only if multi-layer (e.g. gstat output)
  if (nlyr(r1) > 1) r1 <- r1[[1]]
  if (nlyr(r2) > 1) r2 <- r2[[1]]
  
  # Align spatial grids
  r2 <- project(r2, crs(r1), method = "bilinear")
  r2 <- resample(r2, r1, method = "bilinear")
  
  # Optional: Min-Max Normalization to [0, 1] range
  if (normalize) {
    # Helper to scale a raster using its min/max values
    min_max_scale <- function(r) {
      rng <- minmax(r) # returns a matrix with min in row 1, max in row 2
      r_min <- rng[1, 1]
      r_max <- rng[2, 1]
      
      if (r_max == r_min) return(r - r_min) # Avoid division by zero if raster is uniform
      return((r - r_min) / (r_max - r_min))
    }
    
    r1 <- min_max_scale(r1)
    r2 <- min_max_scale(r2)
    
    #change title to normalized
    title_suffix <- " (Normalized [0,1])"
  } else {
    title_suffix <- ""
  }
  
  diff_r <- r1 - r2
  diff_vals <- as.vector(values(diff_r, na.rm = TRUE))
  
  terra::plot(diff_r,
              main = paste0("Difference: ", r1_name, " minus ", r2_name, title_suffix),
              col  = hcl.colors(100, "RdBu"))
  
  cat(" Mean diff: ", mean(diff_vals))
  cat(" Mean Absolute error: ", mean(abs(diff_vals)))
  cat(" RMSE: ", sqrt(mean(diff_vals^2)))
  cat(" Min / Max: ", min(diff_vals), max(diff_vals))
  
  
  invisible(diff_r)
}
