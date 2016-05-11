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
