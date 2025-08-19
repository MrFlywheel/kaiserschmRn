rad_Cor_P4M <- function(spath, rpath, spectral_file){
  kaiserschmRn::package.install(c('foreign', 'estimatr', 'piledge/pizzR'))
    library(foreign)
    library(estimatr)

    targets <- terra::vect(spath)
    rst <- terra::rast(rpath)
    rst <- rst[[1:5]]
    cat('Plotting input raster. Is it correct?')
    terra::plot(rst)
    targets[['id']] <- NULL

    extr <- terra::extract(rst, targets, method='simple', fun='mean')
    extr <- cbind(ID=paste0(targets$Farbe, '_', targets$Target_ID), extr[-1])

    spec_res2 <- read.table(spectral_file, dec = ".", sep = ";", header = T)

    coefficients <- data.frame(Band=c('Blue', 'Green', 'Red', 'RedEdge', 'NIR'), Coef. = '')
    #band definition P4M
    band_blue <- seq(450-16,450+16,1)
    band_green <- seq(560-16,560+16,1)
    band_red <- seq(650-16,650+16,1)
    band_re <- seq(730-16,730+16,1)
    band_nir <- seq(840-16,840+16,1)

    #Blue
    refl_blue <- colMeans(spec_res2[spec_res2[,1]%in%band_blue,])
    #plot(extr[match(names(refl_blue[-1]),extr[,1]),2],refl_blue[-1], xlim = c(0, 65535), ylim=c(0, 300))
    lm_blue <- lm(refl_blue[-1]~0+extr[match(names(refl_blue[-1]),extr[,1]),2])
    coefficients$Coef.[1] <- lm_blue$coefficients[1]

    #Green
    refl_green <- colMeans(spec_res2[spec_res2[,1]%in%band_green,])
    lm_green <- lm(refl_green[-1]~0+extr[match(names(refl_green[-1]),extr[,1]),3])
    coefficients$Coef.[2] <- lm_green$coefficients[1]

    #Red
    refl_red <- colMeans(spec_res2[spec_res2[,1]%in%band_red,])
    lm_red <- lm(refl_red[-1]~0+extr[match(names(refl_red[-1]),extr[,1]),4])
    coefficients$Coef.[3] <- lm_red$coefficients[1]

    #Red Edge
    refl_RE <- colMeans(spec_res2[spec_res2[,1]%in%band_re,])
    lm_RE <- lm(refl_RE[-1]~0+extr[match(names(refl_RE[-1]),extr[,1]),5])
    coefficients$Coef.[4] <- lm_blue$coefficients[1]

    #NIR
    refl_NIR <- colMeans(spec_res2[spec_res2[,1]%in%band_nir,])
    lm_NIR <- lm(refl_NIR[-1]~0+extr[match(names(refl_NIR[-1]),extr[,1]),6])
    coefficients$Coef.[5] <- lm_blue$coefficients[1]

    ## Raster correction
    rst_cor <- rst * as.numeric(coefficients$Coef.)
    terra::writeRaster(rst_cor, paste0(pizzR::file_path_sans_ext(rpath), '_rad_cor.tif'))
}
