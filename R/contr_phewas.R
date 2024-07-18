#' Fill in later
#'
#' @param x A contrasts function used to contrast a column and rename the column
#'
#' @return A data frame with the contrasts and a renamed column name
#' 
#'
#' 
contr.phewas=function(x){
  y=contrasts(x)
  colnames(y)=paste0('-',colnames(y))
  y
}