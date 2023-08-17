#' Sum or count distinct
#'
#' FILL IN LATER
#'
#' @param index an index
#'
#' @return an index or a length of an index
#' @export
#'
#' @examples default_code_agg()
default_code_agg <- function(index) {
  #If it's a number, sum, otherwise count distinct
  if(is.numeric(index)==TRUE) {
    sum(index)
  } else {
    length(unique(index))
  }
}
