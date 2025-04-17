search_strings <- function(data,  expr=NULL){
  n_expr <- length(expr)
  rel <- Filter(function(x) {
    all(sapply(expr, function(e) grepl(e, x)))
  }, data)
  return(rel)
}
