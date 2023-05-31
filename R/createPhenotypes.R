createPhenotypes <-
  function(id.vocab.code.index, min.code.count=2, add.phecode.exclusions=T, translate=T, id.sex, 
           full.population.ids=unique(id.vocab.code.index[[1]]),
           aggregate.fun=PheWAS:::default_code_agg,
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS::phecode_rollup_map,
           exclusion.map=PheWAS::phecode_exclude,
           sex.restriction=PheWAS::sex_restriction,
           map.codes.make.distinct = FALSE)
  {
    id.name=names(id.vocab.code.index)[1]
    
    #Warn if id.sex information is not provided.
    if(missing(id.sex)) { warning("It is recommended to provide id.sex information to help address spurious sex-specific associations.") }
    
    if(!translate) {
      #Warn about exclusions if input is not translated and not phecodes. Same with id.sex
      if(add.phecode.exclusions & sum(tolower(id.vocab.code.index[[2]])=='phecode')!=nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")}
      if(!missing(id.sex) & sum(tolower(id.vocab.code.index[[2]])=='phecode')!=nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but id.sex is supplied for sex-based exclusions. Ensure that the code column has only phecodes or omit id.sex for accurate results.")}
      phemapped=tbl_df(data.frame(id=id.vocab.code.index[[1]],code=id.vocab.code.index[[3]],index=id.vocab.code.index[[4]],stringsAsFactors = F))
    } else {
      #check to make sure numeric codes were not passed in
      if(!class(id.vocab.code.index[[3]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      names(id.vocab.code.index)=c("id","vocabulary_id","code","index")
      message("Mapping codes to phecodes...")
      phemapped=mapCodesToPhecodes(id.vocab.code.index, make.distinct=map.codes.make.distinct, vocabulary.map=vocabulary.map, rollup.map=rollup.map) %>% transmute(id, code=phecode, index)
    }
    
    #Warn if there are ICD9CM and ICD10CM codes and numeric indexes, but default aggregation and code count settings
    if(missing(aggregate.fun) & map.codes.make.distinct == FALSE & min.code.count==2 &
       is.numeric(id.vocab.code.index[[4]]) & 
       "ICD9CM" %in% id.vocab.code.index[[2]] & "ICD10CM" %in% id.vocab.code.index[[2]] ) {
         warning("You are using ICD9CM and ICD10CM codes and numeric counts as an index. With default settings, individuals dual-coded on a single day with ICD9CM and ICD10CM codes that indicate a single phecode will meet a minimum code count of 2, as each code counts separately. Consider providing dates as an index instead. If that's not possible, be sure to evaluate the scientific impacts and consider adjusting `min.code.count`, `map.codes.make.distinct`, or `aggregate.fun`.")
       }
    
    message("Aggregating codes...")
    phecode=ungroup(summarize(group_by(phemapped,id,code),count=aggregate.fun(index)))
    phecode=phecode[phecode$count>0,]
    
    #Check exclusions, and add them to the list
    if(add.phecode.exclusions) {
      message("Mapping exclusions...")
      exclusions = inner_join(phecode %>% rename(exclusion_criteria=code), exclusion.map, by = "exclusion_criteria")
      exclusions = exclusions %>%  transmute(id, code, count=-1) %>% distinct()
      phecode=rbind(phecode,exclusions)
    }
    
    #If there is request for a min code count, adjust counts to -1 if needed
    if(!is.na(min.code.count)&(max(!is.na(phecode$count)&phecode$count<min.code.count))) {
      phecode[!is.na(phecode$count)&phecode$count<min.code.count,]$count=-1
    } 
    
    if(!is.na(min.code.count)|add.phecode.exclusions) {
      message("Coalescing exclusions and min.code.count as applicable...")
      phecode=ungroup(summarize(group_by(phecode,id,code),count=max(count)))
    }

    message("Reshaping data...")
    phens=spread(phecode,code,count,fill=0)
    
    #Set exclusions to NA, preserving IDs just in case one is -1
    tmp_id=phens[,1]
    phens[phens==-1]=NA
    phens[,1]=tmp_id
    
    #Add in inds present in input or the full population list, but without mapped phecodes
    missing_ids=setdiff(full.population.ids,phens[["id"]])
    if(length(missing_ids)>0) {
      empty_record=phens[1,-1]
      empty_record[]=0
      phens=rbind(phens,data.frame(id=missing_ids,empty_record,check.names=F))
    }
    
    #Change to logical if there is a min code count
    if(!is.na(min.code.count)) {phens[,-1]=phens[,-1]>0}
    
    
    #If there are sex restrictions, set them to NA
    if(!missing(id.sex)) {
      phens=restrictPhecodesBySex(phens,id.sex,sex.restriction)
    }
    
    #Limit to full population ids
    phens = filter(phens, id %in% full.population.ids)
    
    #Rename the ID column to the input ID column name
    names(phens)[1]=id.name
    
    phens
  }
