#' Fill in later
#'
#' @param x N/A
#'
#' @return N/A
#' @export
#'
#' @examples contr.phewas(n/a)
contr.phewas=function(x){
  y=contrasts(x)
  colnames(y)=paste0('-',colnames(y))
  y
}
