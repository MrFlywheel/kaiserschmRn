red_Corr <- function (data, class_column = 1, cutoff = 0.8){
  cor_matrix <- cor(data[-class_column])
  high_corr <- findCorrelation(cor_matrix, cutoff = cutoff)
  data_littlecorr <- cbind(data[class_column], data[-high_corr])
  return(data_littlecorr)
}


