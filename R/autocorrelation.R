autocorrelation <- function(rast, sample_size = 10000, cut_off = 50, width = 1){
  library(gstat)

  rast_len <- 1:ncell(rast)
  pix_nonNA <- rast_len[!is.na(values(rast))]
  varsamps <- sample(pix_nonNA, sample_size)
  sampled_sqr <- data.frame(
    x = xyFromCell(rast, varsamps)[, 1],
    y = xyFromCell(rast, varsamps)[, 2],
    value = terra::extract(rast, varsamps))
  names(sampled_sqr)[3] <- 'rast_value'
  sp::coordinates(sampled_sqr) <- ~ x + y
  varsqr <- variogram(rast_value ~ 1, data=sampled_sqr, cutoff = cut_off, width = width)
  plot(varsqr)
}
