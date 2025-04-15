remove_Inf_Nan <- function(metric){
  metric[is.infinite(metric)] <- NA
  metric[is.nan(metric)] <- NA
  metric <- c(metric)
}
