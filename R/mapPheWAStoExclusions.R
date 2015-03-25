mapPheWAStoExclusions <-
  function(phewas.codes, ids) {
    if(missing(ids)) {
      input = tbl_df(data.frame(id=0,exclusion_criteria=phewas.codes,stringsAsFactors = FALSE))
    }
    else {
      input = tbl_df(data.frame(id=ids, exclusion_criteria=phewas.codes,stringsAsFactors = FALSE))
    }
    output = inner_join(input,phewas_exclude)
    output = output %>% transmute(id, exclusion=code) %>% distinct()
    output
  }
