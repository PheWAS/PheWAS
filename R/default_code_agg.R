default_code_agg <- function(index) {
  #If it's a number, sum, otherwise count distinct
  if(is.numeric(index)==TRUE) {
    sum(index)
  } else {
    length(unique(index))
  }
}