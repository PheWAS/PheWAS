mapPheWAStoExclusions <-
  function(phewas.codes, ids=NA) {
    if(is.na(ids)[1]) {
      input = data.frame(exclusion_criteria=phewas.codes)
    }
    else {
      input = data.frame(id=ids, exclusion_criteria=phewas.codes)
    }
    output = inner_join(input,phewas_exclude,by="exclusion_criteria")
    output = output[,!(names(output)=="exclusion_criteria")]
    names(output)[names(output)=="code"]="phe"
    output=unique(output)
    output
  }
