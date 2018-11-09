mapCodesToPhecodes <-
  function(input, 
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS::phecode_rollup_map,
           make.distinct=TRUE) {
    if(sum(names(input) %in% c("vocabulary_id","code"))!=2) {
      stop("Must supply a data frame with 'vocabulary_id' and 'code' columns")
    }
    if(!class(input[["code"]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
    
    if(!is.null(vocabulary.map)){
      #Perform the direct map
      output = inner_join(input,vocabulary.map,by=c("vocabulary_id","code"))
      #Remove old columns
      output = output %>% select(-code,-vocabulary_id) %>% rename(code=phecode) 
    } else {
      #Warn if the vocabulary IDs are not phecodes
      if(sum(input$vocabulary_id!="phecode")!=0) {warning("Phecode mapping was not requested, but the vocabulary_id of all codes is not 'phecode'")}
      #Prepare for just the phecode expansion
      output=input %>% filter(vocabulary_id=="phecode") %>% select(-vocabulary_id)
    }
    #Make distinct
    if(make.distinct) {output = distinct(output)}
    #Perform the rollup
    if(!is.null(rollup.map)) {
      output = inner_join(output ,rollup.map,by="code") %>% select(-code) %>% rename(phecode=phecode_unrolled)
      #Make distinct
      if(make.distinct) {output = distinct(output)}
    }
    else{
      #Rename output column to phecode
      output = rename(phecode=code)
    }
    
    #Return the output
    output 
  }
