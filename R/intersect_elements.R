intersect_elements <- function(vectors){
  stopifnot('Input must be a vector of integer elements' = typeof(vectors)=='integer')
  vectors <- list(vectors)
  uniqueID <- Reduce(intersect, vectors)
  return(uniqueID)
}
