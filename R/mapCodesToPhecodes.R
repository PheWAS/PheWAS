mapCodesToPhecodes <-
  function(input, 
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS:::phecode_rollup_map,
           make.distinct=T) {
    if(sum(names(input) %in% c("vocabulary_id","code"))!=2) {
      stop("Must supply a data frame with 'vocabulary_id' and 'code' columns")
    }
    if(!class(input[,"code"]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
    #Perform the direct map
    output = inner_join(input,vocabulary.map,by=c("vocabulary_id","code"))
    #Remove old columns
    output = output %>% mutate(vocabulary_id="phecode") %>% select(-code) %>% rename(code=phecode)
    #Make distinct
    if(make.distinct) {output = distinct(output)}
    #Perform the rollup
    output = inner_join(output,phecode_rollup_map,by="code") %>% select(-code) %>% rename(code=phecode_unrolled)
    #Make distinct
    if(make.distinct) {output = distinct(output)}
    
    #Return the output
    output 
  }
