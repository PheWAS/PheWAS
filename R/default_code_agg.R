default_code_agg <- function(index) {
  #If it's a number, sum, otherwise count distinct
  if(class(index)=="numeric") {
    sum(index)
  } else {
    length(unique(index))
  }
}