num_distinct_dates <- function(index) {
  length(distinct(as.Date(index)))
}