mapPheWAStoExclusions <-
  function(phewas.codes, ids=NA) {
    if(is.na(ids)[1]) {
      input = data.frame(exclusion_criteria=phewas.codes)
    }
    else {
      input = data.frame(id=ids, exclusion_criteria=phewas.codes)
    }
    output = merge(input,phewas_exclude)
    output = output[,!(names(output)=="exclusion_criteria")]
    names(output)[names(output)=="code"]="exclusion"
    output=unique(output)
    output
  }
