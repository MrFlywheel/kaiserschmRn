
extract_metrics <- function(data_path, type=NULL, date = NULL, polygons=NULL, stepsize = 0.1, out_path = NULL){
  library(terra)
  library(lidR)

    stopifnot("Type needs to specified. Use 'lidar' for LiDAR input or 'ortho' for orthofoto (5 bands from DJI P4M) or 'texture' for a eight band Haralick texture raster." = !is.null(type))

    metric.name <- function(name, date){
      return(c(paste(date, aufloesung, 'b1', name, sep = "_"), paste(date, aufloesung, 'b2', name, sep = "_"),
               paste(date, aufloesung, 'b3', name, sep = "_"), paste(date, aufloesung, 'b4', name, sep = "_"),
               paste(date, aufloesung, 'b5', name, sep = "_")))
    }

    texture.name <- function(name, date){
      return(c(paste(date, aufloesung, 'energy', name, sep = "_"), paste(date, aufloesung, 'entropy', name, sep = "_"),
               paste(date, aufloesung, 'correlation', name, sep = "_"), paste(date, aufloesung, 'invDifMom', name, sep = "_"),
               paste(date, aufloesung, 'Inertia', name, sep = "_"), paste(date, aufloesung, 'ClusterShade', name, sep = "_"),
               paste(date, aufloesung, 'ClusterProm', name, sep = "_"), paste(date, aufloesung, 'HaralickCorr', name, sep = "_")))
    }

    VI.name <- function(VI, date){
      return(c(paste(date, aufloesung, VI, sep = "_")))
    }

  get.metrics <- function(data, date){




  data$ID <- NULL

  data.anzahl <- matrix(nrow(data))
  colnames(data.anzahl) <- paste(date, aufloesung, 'count', sep = "_")
  data.means <- matrix(unlist(lapply(data, FUN=mean, na.rm = T)), ncol=5, byrow = T)
  colnames(data.means) <- metric.name('mean', date)
  data.sd <- matrix(unlist(lapply(data, FUN=sd, na.rm = T)), ncol=5, byrow = T)
  colnames(data.sd) <- metric.name('sd', date)
  data.var <- matrix(unlist(lapply(data, FUN=var, na.rm = T)), ncol=5, byrow = T)
  colnames(data.var) <- metric.name('var', date)
  data.skew <- matrix(unlist(lapply(data, FUN= skewness, na.rm = T)), ncol=5, byrow = T)
  colnames(data.skew) <- metric.name('skew', date)
  data.kurt <- matrix(unlist(lapply(data, FUN= kurtosis, na.rm = T)), ncol=5, byrow = T)
  colnames(data.kurt) <- metric.name('kurt', date)


  data.q90 <- matrix(unlist(lapply(data, quantile, probs = c(0.90), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q90) <- metric.name('q90', date)
  data.q80 <- matrix(unlist(lapply(data, quantile, probs = c(0.80), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q80) <- metric.name('q80', date)
  data.q70 <- matrix(unlist(lapply(data, quantile, probs = c(0.70), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q70) <- metric.name('q70', date)
  data.q60 <- matrix(unlist(lapply(data, quantile, probs = c(0.60), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q60) <- metric.name('q60', date)
  data.q50 <- matrix(unlist(lapply(data, quantile, probs = c(0.50), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q50) <- metric.name('q50', date)
  data.q40 <- matrix(unlist(lapply(data, quantile, probs = c(0.40), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q40) <- metric.name('q40', date)
  data.q30 <- matrix(unlist(lapply(data, quantile, probs = c(0.30), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q30) <- metric.name('q30', date)
  data.q20 <- matrix(unlist(lapply(data, quantile, probs = c(0.20), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q20) <- metric.name('q20', date)
  data.q10 <- matrix(unlist(lapply(data, quantile, probs = c(0.10), na.rm = T)), ncol=5, byrow = T)
  colnames(data.q10) <- metric.name('q10', date)

  NDVI <- as.data.frame((data.means[, 5] - data.means[, 3]) / (data.means[, 5] + data.means[, 3]))
  colnames(NDVI) <- VI.name('NDVI', date)
  NDVI_NIR_G <- as.data.frame(((data.means[, 5] - data.means[, 2]) / (data.means[, 5] + data.means[, 2])))
  colnames(NDVI_NIR_G) <- VI.name('NDVI_NIR_G', date)
  NDVI_NIR_B <- as.data.frame((data.means[, 5] - data.means[, 1]) / (data.means[, 5] + data.means[, 1]))
  colnames(NDVI_NIR_B) <- VI.name('NDVI_NIR_B', date)
  NDVI_NIR_RE<- as.data.frame((data.means[, 5] - data.means[, 4]) / (data.means[, 5] + data.means[, 4]))
  colnames(NDVI_NIR_RE) <- VI.name('NDVI_NIR_RE', date)
  NDVI_RE_B <- as.data.frame((data.means[, 4] - data.means[, 1]) / (data.means[, 4] + data.means[, 1]))
  colnames(NDVI_RE_B) <- VI.name('NDVI_RE_B', date)
  NDVI_B_R <- as.data.frame((data.means[, 1] - data.means[, 3]) / (data.means[, 1] + data.means[, 3]))
  colnames(NDVI_B_R) <- VI.name('NDVI_B_R', date)
  GNDVI <- as.data.frame((data.means[, 5] - data.means[, 3]) / (data.means[, 5] + data.means[, 3]))
  colnames(GNDVI) <- VI.name('GNDVI', date)
  RENDVI <- as.data.frame((data.means[, 5] - data.means[, 4]) / (data.means[, 5] + data.means[, 4]))
  colnames(RENDVI) <- VI.name('RENDVI', date)
  NDWI <- as.data.frame((data.means[, 2] - data.means[, 5]) / (data.means[, 2] + data.means[, 5]))
  colnames(NDWI) <- VI.name('NDWI', date)
  CVI <- as.data.frame((data.means[, 5] * (data.means[, 3] / (data.means[, 2] ** 2))))
  colnames(CVI) <- VI.name('CVI', date)
  EVI <- as.data.frame((2.5 * ((data.means[, 5] - data.means[, 3]) / ((data.means[, 5] + 6* data.means[, 3] - 7.5 * data.means[, 1]) + 1))))
  colnames(EVI) <- VI.name('EVI', date)
  CCCI <- as.data.frame((((data.means[, 5] - data.means[, 4]) / (data.means[, 5] + data.means[, 4])) / ((data.means[, 5] - data.means[, 3]) / (data.means[, 5] + data.means[, 3]))))
  colnames(CCCI) <- VI.name('CCCI', date)
  CCCI_NIR_B_R <- as.data.frame((((data.means[, 5] - data.means[, 1]) / (data.means[, 5] + data.means[, 1]))/((data.means[, 5] - data.means[, 3])/(data.means[, 5] + data.means[, 3]))))
  colnames(CCCI_NIR_B_R) <- VI.name('CCCI_NIR_B_R', date)


  return(data.frame(cbind(data.anzahl, data.means, data.sd, data.var, data.q90, data.q80,
                          data.q70, data.q60, data.q50, data.q40, data.q30, data.q20, data.q10,
                          NDVI, NDVI_NIR_G, NDVI_NIR_B, NDVI_NIR_RE, NDVI_RE_B, NDVI_B_R, GNDVI,
                          NDWI, RENDVI, CVI, EVI, CCCI, CCCI_NIR_B_R)))
}

  get.texture <- function(data_texture, date){   #die Sache hier ist: Die Textur ist als mehrkanaliges tif zu laden. Dh die funktion metric.name funktioniert hier nicht, da diese f?r 5 b?nder ausgelegt ist.

  data_texture$ID <- NULL

  data_texture.means <- matrix(unlist(lapply(data_texture, FUN=mean, na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.means) <- texture.name('mean', date)
  data_texture.sd <- matrix(unlist(lapply(data_texture, FUN=sd, na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.sd) <- texture.name('sd', date)
  data_texture.var <- matrix(unlist(lapply(data_texture, FUN=var, na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.var) <- texture.name('var', date)

  data_texture.q90 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.90), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q90) <- texture.name('q90', date)
  data_texture.q80 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.80), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q80) <- texture.name('q80', date)
  data_texture.q70 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.70), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q70) <- texture.name('q70', date)
  data_texture.q60 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.60), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q60) <- texture.name('q60', date)
  data_texture.q50 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.50), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q50) <- texture.name('q50', date)
  data_texture.q40 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.40), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q40) <- texture.name('q40', date)
  data_texture.q30 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.30), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q30) <- texture.name('q30', date)
  data_texture.q20 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.20), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q20) <- texture.name('q20', date)
  data_texture.q10 <- matrix(unlist(lapply(data_texture, quantile, probs = c(0.10), na.rm = T)), ncol=8, byrow = T)
  colnames(data_texture.q10) <- texture.name('q10', date)



  return(data.frame(cbind(data_texture.means, data_texture.sd, data_texture.var, data_texture.q90, data_texture.q80,
                          data_texture.q70, data_texture.q60, data_texture.q50, data_texture.q40, data_texture.q30, data_texture.q20, data_texture.q10)))
}

  get.lidar.metrics.RGB <- function(data, stepsize, date, PC_type=2){

  stopifnot("stepsize must be between 0 and 1" = stepsize > 0 & stepsize < 1)
  stopifnot("'PC_type' must be equal to 1 (for .txt file) or 2 (for .las file)" = !is.null(PC_type))
  stopifnot("Input'date' must be a vector of length 1 with the date of image aquisition" = typeof(date) == 'character' & length(date)==1)

  #data <- filter_poi(data, Z > 0)
  quantiles = seq(0 + stepsize, 1 - stepsize, by = stepsize)
  z_q <- matrix(unlist(quantile(data$Z, probs = quantiles)), ncol=length(quantiles), byrow = T)
  colnames(z_q) <- paste0(date, '_z_q',quantiles)

  z_iqr <- as.data.frame(IQR(data$Z))
  colnames(z_iqr) <- paste0(date, '_z_iqr')
  z_mean <- as.data.frame(mean(data$Z))
  colnames(z_mean) <- paste0(date, '_z_mean')
  z_var <- as.data.frame(var(data$Z))
  colnames(z_var) <- paste0(date,'_z_var')
  z_sd <- as.data.frame(sd(data$Z))
  colnames(z_sd) <- paste0(date,'_z_sd')
  z_skew <- as.data.frame(skewness(data$Z))
  colnames(z_skew) <- paste0(date,'_z_skew')
  z_kurt <- as.data.frame(kurtosis(data$Z))
  colnames(z_kurt) <- paste0(date,'_z_kurt')
  z_MAD <- as.data.frame(mad(data$Z))
  colnames(z_MAD) <- paste0(date,'_z_MAD')

  z_VCI <- matrix(lidR::VCI(data$Z + abs(min(data$Z)) , zmax = max(data$Z+ abs(min(data$Z)))), dimnames = list('', paste0('VCI_', date)))
  z_Dens <- matrix(unlist(stdmetrics_z(data$Z + abs(min(data$Z)) , zmin = min(data$Z + abs(min(data$Z)))))[28:36], ncol = 9, byrow = T, dimnames = list('', c(paste0('z_D_', date, 1:9 )))) # zmin wird standardmäßig als 0 angenommen, was nur für normalisierte PW gilt!
  z_LAD <- lidR::LAD(data$Z) #wegen unterschiedlicher Baumhöhen würden features je nach Baum unterschiedlich sein (zb Götterbaum hat kein z25. Relative Baumhöhen zw. 0 und 1 funktieren auch nicht, da keine 1m bins. Daher nur für LAI verwendet)
  z_LAI <- sum(z_LAD$lad, na.rm = T)
  names(z_LAI) <- paste0('z_LAI_', date)

  Ch_to_Cw_X <- max(data$Z) / (max(data$X) - min(data$X))
  Ch_to_Cw_Y <- max(data$Z) / (max(data$Y) - min(data$Y))
  z_RI <- rumple_index(data$X, data$Y, data$Z)

  if(PC_type == 2){

    i_q <- matrix(unlist(quantile(data$Intensity, probs = quantiles)), ncol=length(quantiles), byrow = T)
    colnames(i_q) <- paste0(date, '_i_q',quantiles)
    i_iqr <- as.data.frame(IQR(data$Intensity))
    colnames(i_iqr) <- paste0(date, '_i_iqr')
    i_mean <- as.data.frame(mean(data$Intensity))
    colnames(i_mean) <- paste0(date, '_i_mean')
    i_var <- as.data.frame(var(data$Intensity))
    colnames(i_var) <- paste0(date,'_i_var')
    i_sd <- as.data.frame(sd(data$Intensity))
    colnames(i_sd) <- paste0(date,'_i_sd')
    i_skew <- as.data.frame(skewness(data$Intensity))
    colnames(i_skew) <- paste0(date,'_i_skew')
    i_kurt <- as.data.frame(kurtosis(data$Intensity))
    colnames(i_kurt) <- paste0(date,'_i_kurt')
    i_MAD <- as.data.frame(mad(data$Intensity))
    colnames(i_MAD) <- paste0(date,'_i_MAD')

    R_q <- matrix(unlist(quantile(data$R, probs = quantiles)), ncol=length(quantiles), byrow = T)
    colnames(R_q) <- paste0(date, '_R_q',quantiles)
    R_iqr <- as.data.frame(IQR(data$R))
    colnames(R_iqr) <- paste0(date, '_R_iqr')
    R_mean <- as.data.frame(mean(data$R))
    colnames(R_mean) <- paste0(date, '_R_mean')
    R_var <- as.data.frame(var(data$R))
    colnames(R_var) <- paste0(date,'_R_var')
    R_sd <- as.data.frame(sd(data$R))
    colnames(R_sd) <- paste0(date,'_R_sd')
    R_skew <- as.data.frame(skewness(data$R))
    colnames(R_skew) <- paste0(date,'_R_skew')
    R_kurt <- as.data.frame(kurtosis(data$R))
    colnames(R_kurt) <- paste0(date,'_R_kurt')
    R_MAD <- as.data.frame(mad(data$R))
    colnames(R_MAD) <- paste0(date,'_R_MAD')

    G_q <- matrix(unlist(quantile(data$G, probs = quantiles)), ncol=length(quantiles), byrow = T)
    colnames(G_q) <- paste0(date, '_G_q',quantiles)
    G_iqr <- as.data.frame(IQR(data$G))
    colnames(G_iqr) <- paste0(date, '_G_iqr')
    G_mean <- as.data.frame(mean(data$G))
    colnames(G_mean) <- paste0(date, '_G_mean')
    G_var <- as.data.frame(var(data$G))
    colnames(G_var) <- paste0(date,'_G_var')
    G_sd <- as.data.frame(sd(data$G))
    colnames(G_sd) <- paste0(date,'_G_sd')
    G_skew <- as.data.frame(skewness(data$G))
    colnames(G_skew) <- paste0(date,'_G_skew')
    G_kurt <- as.data.frame(kurtosis(data$G))
    colnames(G_kurt) <- paste0(date,'_G_kurt')
    G_MAD <- as.data.frame(mad(data$G))
    colnames(G_MAD) <- paste0(date,'_G_MAD')

    B_q <- matrix(unlist(quantile(data$B, probs = quantiles)), ncol=length(quantiles), byrow = T)
    colnames(B_q) <- paste0(date, '_B_q',quantiles)
    B_iqr <- as.data.frame(IQR(data$B))
    colnames(B_iqr) <- paste0(date, '_B_iqr')
    B_mean <- as.data.frame(mean(data$B))
    colnames(B_mean) <- paste0(date, '_B_mean')
    B_var <- as.data.frame(var(data$B))
    colnames(B_var) <- paste0(date,'_B_var')
    B_sd <- as.data.frame(sd(data$B))
    colnames(B_sd) <- paste0(date,'_B_sd')
    B_skew <- as.data.frame(skewness(data$B))
    colnames(B_skew) <- paste0(date,'_B_skew')
    B_kurt <- as.data.frame(kurtosis(data$B))
    colnames(B_kurt) <- paste0(date,'_B_kurt')
    B_MAD <- as.data.frame(mad(data$B))
    colnames(B_MAD) <- paste0(date,'_B_MAD')


    VDVI_mean <- (2*data$G - data$R - data$B) / (2*data$G + data$R + data$B)
    VDVI_mean <- kaiserschmRn::remove_Inf_Nan(VDVI_mean)
    VDVI_mean <- as.data.frame(mean(VDVI_mean, na.rm = T))
    colnames(VDVI_mean) <- 'VDVI_mean'
    VDVI_sd <- (2*data$G - data$R - data$B) / (2*data$G + data$R + data$B)
    VDVI_sd <-  kaiserschmRn::remove_Inf_Nan(VDVI_sd)
    VDVI_sd <- as.data.frame(sd(VDVI_sd, na.rm = T))
    colnames(VDVI_sd) <- 'VDVI_sd'
    NGRDI_mean <- (data$G - data$R) / (data$G + data$R)
    NGRDI_mean <- kaiserschmRn::remove_Inf_Nan(NGRDI_mean)
    NGRDI_mean <- as.data.frame(mean(NGRDI_mean, na.rm = T))
    colnames(NGRDI_mean) <- 'NGRDI_mean'
    NGRDI_sd <- (data$G - data$R) / (data$G + data$R)
    NGRDI_sd <- kaiserschmRn::remove_Inf_Nan(NGRDI_sd)
    NGRDI_sd <- as.data.frame(sd(NGRDI_sd, na.rm = T))
    colnames(NGRDI_sd) <- 'NGRDI_sd'
    VARI_mean <- (data$G - data$R) / (data$G + data$R - data$B)
    VARI_mean <- kaiserschmRn::remove_Inf_Nan(VARI_mean)
    VARI_mean <- as.data.frame(mean(VARI_mean, na.rm = T))
    colnames(VARI_mean) <- 'VARI_mean'
    VARI_sd <- (data$G - data$R) / (data$G + data$R - data$B)
    VARI_sd <- kaiserschmRn::remove_Inf_Nan(VARI_sd)
    VARI_sd <- as.data.frame(sd(VARI_sd, na.rm = T))
    colnames(VARI_sd) <- "VARI_sd"
    GRRI_mean <- (data$G / data$R)
    GRRI_mean <- kaiserschmRn::remove_Inf_Nan(GRRI_mean)
    GRRI_mean <- as.data.frame(mean(GRRI_mean, na.rm = T))
    colnames(GRRI_mean) <- "GRRI_mean"
    GRRI_sd <- (data$G / data$R)
    GRRI_sd <- kaiserschmRn::remove_Inf_Nan(GRRI_sd)
    GRRI_sd <- as.data.frame(sd(GRRI_sd, na.rm = T))
    colnames(GRRI_sd) <- "GRRI_sd"
  }

  if(PC_type==1){
    return(data.frame(cbind(z_mean, z_sd, z_var, z_skew, z_kurt, z_q, z_iqr, z_VCI, z_Dens, z_LAI, Ch_to_Cw_X, Ch_to_Cw_Y, z_RI)))
  }
  if(PC_type==2){
    return(data.frame(cbind(z_mean, z_sd, z_var, z_skew, z_kurt, z_q, z_iqr, z_VCI, z_Dens, z_LAI,Ch_to_Cw_X, Ch_to_Cw_Y, z_RI,
                            i_mean, i_sd, i_var, i_skew, i_kurt, i_MAD, i_q, i_iqr, R_mean, R_sd, R_var, R_skew, R_kurt, R_MAD, R_q, R_iqr,
                            G_mean, G_sd, G_var, G_skew, G_kurt, G_MAD, G_q, G_iqr, B_mean, B_sd, B_var, B_skew, B_kurt, B_MAD, B_q, B_iqr,
                            VDVI_mean,  VDVI_sd, NGRDI_mean, NGRDI_sd, VARI_mean, VARI_sd, GRRI_mean, GRRI_sd)))
  }
}

if(type == 'lidar'){


  if(is.null(out_path)) out_path <- file.path(dirname(data_path), '/export')
  pizzR::setcreate.wd(out_path)


  for (i1 in seq_along(data_path)){
    segments <- list.files(data_path[i1], pattern = '.las', full.names = F)
    segment_names <- segments[order(nchar(segments), segments)]
    seg_ID <- pizzR::file_path_sans_ext(segment_names)
    pc_files <- paste0(data_path, segment_names)

    result <- data.frame()

    for (i2 in seq_along(pc_files)){
      pizzR::loop_progress(i2, 3)
      data <- readLAS(pc_files[i2])
      lidar_metrics <- get.lidar.metrics.RGB(data, stepsize=0.1, data_type = 2, date = date[i1])
      result <- rbind(result, lidar_metrics)
    }

    result <- cbind(seg_ID, result)
    colnames(result)[1] <- 'Segment_ID'
    write.csv2(result, paste0(date[i1],'_LiDAR_metrics.csv'))
  }
  if(anyNA(result)) print('NAs found in the final extract!')
}
if(type == 'ortho'){
  stopifnot("'polygons' must be specified in order to compute metrics for the areas inside the polygons. Please provide a full path to a shp-layer. First column in attribute table needs to be ID" = !is.null(polygons))
  stopifnot('must be the same number of paths for raster and shape layers. Use shape paths multiple times if it is the same layer'= length(data_path)!=length(polygons))

  if(is.null(out_path)) out_path <- file.path(dirname(data_path), '/export')
  pizzR::setcreate.wd(out_path)

  for (i1 in seq_along(data_path)){

    rst <- terra::rast(data_path[i1])
    rst <- rst[[1:5]]
    shp <- terra::vect(polygons[i1])
    shp <- shp[order(shp[[1]])]
    checkID <- shp[[1]]
    aufloesung <- round(res(rst)[1], 3)

    extr <- terra::extract(rst, shp)
    cat(paste(date[i1], aufloesung,  ': extract successful. Continuing with texture extract...'))

    result <- data.frame()
    for (i2 in seq_along(checkID)){
      cat(paste0('/r', Sys.time(),' Loops remainings: ', length(checkID) - i2 + 1), '    ')
      data <- subset(extr, extr$ID==i2)
      result.tmp <- get.metrics(data, date[i1])
      result <- rbind(result, result.tmp)
    }
    result <- cbind(checkID,  result)
    colnames(result)[1] <- 'Segment_ID'

    feather::write_feather(result, paste(date[i1], aufloesung, 'extract.feather', sep = "_"))
    write.csv2(result, paste(date[i1], aufloesung, 'extract.csv', sep = "_"))
  }
  if(anyNA(result)) print('NAs found in the final extract!')

}
if(type == 'texture'){
  stopifnot("'polygons' must be specified in order to compute metrics for the areas inside the polygons. Please provide a full path to a shp-layer. First column in attribute table needs to be ID" = !is.null(polygons))
  stopifnot('must be the same number of paths for raster and shape layers. Use shape paths multiple times if it is the same layer'= length(data_path)!=length(polygons))

  if(is.null(out_path)) out_path <- file.path(dirname(data_path), '/export')
  pizzR::setcreate.wd(out_path)

  for (i1 in seq_along(data_path)){

    rst <- terra::rast(data_path[i1])
    rst <- rst[[1:5]]
    shp <- terra::vect(polygons[i1])
    shp <- shp[order(shp[[1]])]
    checkID <- shp[[1]]
    aufloesung <- round(res(rst)[1], 3)

    extr <- terra::extract(rst, shp)
    cat(paste(date[i1], aufloesung,  ': extract successful. Continuing with texture extract...'))

    result <- data.frame()
    for (i2 in seq_along(checkID)){
      cat(paste0('/r', Sys.time(),' Loops remainings: ', length(checkID) - i2 + 1), '    ')
      data <- subset(extr, extr$ID==i2)
      result.tmp <- get.texture(data_texture = data, date[i1])
      result <- rbind(result, result.tmp)
    }
    result <- cbind(checkID,  result)
    colnames(result)[1] <- 'Segment_ID'

    feather::write_feather(result, paste(date[i1], aufloesung, 'extract.feather', sep = "_"))
    write.csv2(result, paste(date[i1], aufloesung, 'texture_extract.csv', sep = "_"))
  }
  if(anyNA(result)) print('NAs found in the final extract!')
  }

}

