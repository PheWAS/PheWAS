mapPheWAStoExclusions <-
  function(phewas.codes, ids) {
    if(missing(ids)) {
      input = data.frame(id=0,exclusion_criteria=phewas.codes,stringsAsFactors = FALSE)
    }
    else {
      input = data.frame(id=ids, exclusion_criteria=phewas.codes,stringsAsFactors = FALSE)
    }
#    output = rbind_all(by(input,input$id,function(by_id){
#      codes=unique(unlist(lapply(ex_split,function(by_range) {
#        match=max(by_id$exclusion_criteria>=by_range$low&
#              #by_id$exclusion_criteria>=by_range$lownum&
#              by_id$exclusion_criteria<=by_range$high#&
#              #by_id$exclusion_criteria<=by_range$highnum
#              )
#        if(match) {by_range$codes}
#        else {c()}
#      })))
#      data.frame(id=by_id$id[1],phe=unique(codes))
#    }))
    #output = output %>% transmute(id,phe=code) %>% distinct()
    exc_groups = rbind_all(by(input,input$id,function(by_id){
       merge(by_id,ex_df,by=c()) %>%
        filter(exclusion_criteria>=low,exclusion_criteria<=high) %>%
        dplyr::select(id,low,high) %>% distinct()
    }))

    output = inner_join(exc_groups,ex_reps) %>%
      dplyr::select(id,phe) %>% distinct()
    output
  }
