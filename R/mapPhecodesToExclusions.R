#' Map phecodes to their exclusions
#'
#'This function maps phecodes (optionally with ids for individuals) to a set of
#' phecode exclusions. It has replaced \code{mapPheWAStoExclusions}.
#'
#' @param phecodes A vector of phecodes.
#' @param ids An optional vector of ids to pair with the provided phecodes.
#'
#' @return A data frame containing phecodes and their exclusions. IDs for those
#' codes and exclusions are included if they were supplied.
#' \item{id}{If ids were provided, the individual ids are included as the first
#' column}
#' \item{exclusion_criteria}{Input phecodes}
#' \item{exclusion}{The exclusion phecodes for the codes provided}

mapPhecodesToExclusions <-
  function(phecodes, ids) {
    if(missing(ids)) {
      input = tbl_df(data.frame(id=0,exclusion_criteria=phecodes,stringsAsFactors = FALSE))
    }
    else {
      input = tbl_df(data.frame(id=ids, exclusion_criteria=phecodes,stringsAsFactors = FALSE))
    }
    output = inner_join(input,phecode_exclude)
    output = output %>% transmute(id, exclusion=code) %>% distinct()
    output
  }
