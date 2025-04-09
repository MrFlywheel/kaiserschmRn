pointCloud_crop <- function(las_path=NULL, shp_path=NULL, out_path=NULL, n_cores=NULL){
  kaiserschmRn::package_install('foreach', 'lidR', 'parallel', 'pizzR')
  library(foreach)

  stopifnot('las_path must be the full path name to a .las file ' = !is.null(las_path) & !is.numeric(las_path))
  stopifnot('shp_path must be the full path name to a .las file ' = !is.null(shp_path) & !is.numeric(las_path))

  las_data <- lidR::readLAScatalog(las_path)
  polygons <- sf::st_read(shp_path)
  if (is.null(out_path)) out_path <- file.path(dirname(las_path), '/export')


  PC_crop  <- function(las_data, polygon) {
    clipped_las <- lidR::clip_roi(las_data, polygon)
    lidR::writeLAS(clipped_las, paste0(polygon$ID_clean, '.las'))
  }

  if (is.null(n_cores)) num_cores <- parallel::detectCores()/2
  cl <- parallel::makeCluster(num_cores)
  registerDoParallel(cl)

  polygon_list <- split(polygons, seq(nrow(polygons)))

  cat(sprintf("%s: Export files to '%s'", pizzR::Systime(), out_path))
  clipped_las_list <- foreach(polygon = polygon_list, .packages = c("lidR", "sf")) %dopar% {
    pizzR::setcreate.wd(out_path)
    PC_crop(las_data, polygon)
    gc(reset = T, full = T)
  }
  stopCluster(cl)
}

