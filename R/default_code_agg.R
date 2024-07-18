#' Sum or count distinct
#'
#' a function that either returns the sum of an input, if numeric, or the length
#' of an input if not numeric 
#'
#' @param index an index
#'
#' @return the sum of an index or a length of an index

default_code_agg <- function(index) {
  #If it's a number, sum, otherwise count distinct
  if(is.numeric(index)==TRUE) {
    sum(index)
  } else {
    length(unique(index))
  }
}
