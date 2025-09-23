roughness_RMSE <- function(las_object, res = 1, dtm = NULL, filter_height = NULL){
  kaiserschmRn::package.install(c('lidR', 'terra'))
  library(lidR)
  library(terra)

  rmse_plane_function <- function(x, y, z) {
    fit <- lm(z ~ x + y)
    z_fit <- predict(fit)
    residuals <- z - z_fit
    rmse <- sqrt(mean(residuals^2, na.rm = TRUE))
    return(rmse)
    }

  assign("rmse_plane_function", rmse_plane_function, envir = .GlobalEnv)

  las <- las_object
  normalised <- any(head(las$Classification == 2) & (head(las$Z)>2))
  stopifnot("las needs to be a height normalised point cloud. Please provide a path to a dtm over the variable 'dtm'" = normalised && !is.null(dtm))
  if(!is.null(dtm)){
    dtm <- terra::rast(dtm)
    las <- las - dtm
  }
  nlas_filtered <- filter_poi(las, Z < filter_height)
  stopifnot("specify a height above ground with variable 'filter_height' which should be still considered to compute roughness" = !is.null(filter_height))

  rmse_raster <- grid_metrics(nlas_filtered, res = res, rmse_plane_function(X, Y, Z))
  plot(rmse_raster)
  return(rmse_raster)
}


