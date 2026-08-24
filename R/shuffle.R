shuffle <- function(x){
  stopifnot('x must by a data.frame' = is.data.frame(x))
  return(x[sample(1:nrow(x), nrow(x)),])
}



