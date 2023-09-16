#' Implemement storage efficient indexing
#'
#' @param x Records in X_1 that share an agreement pattern
#' @param R Maximum number of labels to keep
#' @return Reduced vector of record labels
#' @export
#'

sei <- function(x, R){
  if(length(x) <= R){
    return(x)
  } else {
    sample(x, R, replace = FALSE)
  }
}
