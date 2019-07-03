default_code_agg <- function(index) {
  #If it's a number, sum, otherwise count distinct
  if(max(class(index)=="numeric")==TRUE) {
    sum(index)
  } else {
    length(unique(index))
  }
}