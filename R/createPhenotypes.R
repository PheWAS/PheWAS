createPhenotypes <-
  function(id.vocab.code.index, min.code.count=2, add.phecode.exclusions=T, translate=T, id.sex, 
           full.population.ids=unique(id.vocab.code.index[[1]]),
           aggregate.fun=PheWAS::num_distinct_dates, 
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS::phecode_rollup_map,
           exclusion.map=PheWAS::phecode_exclude)
  {
    id.name=names(id.vocab.code.index)[1]
    
    if(!translate) {
      #Warn about exclusions if input is not translated and not phecodes. Same with id.sex
      if(add.phecode.exclusions & sum(tolower(id.vocab.code.index[,2])=='phecode')==nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but exclusions are to be applied. Ensure that the code column has only phecodes or disable add.phecode.exclusions for accurate results.")}
      if(!missing(id.sex) & sum(tolower(id.vocab.code.index[,2])=='phecode')==nrow(id.vocab.code.index)){stop("Codes are not translated and vocab is not 'phecode' for every row, but id.sex is supplied for sex-based exclusions. Ensure that the code column has only phecodes or omit id.sex for accurate results.")}
      phemapped=tbl_df(data.frame(id=id.vocab.code.index[,1],code=id.vocab.code.index[,3],index=id.vocab.code.index[,4]))
    } else {
      #check to make sure numeric codes were not passed in
      if(!class(id.vocab.code.index[[3]]) %in% c("character","factor")) {stop("Please ensure character or factor code representation. Some vocabularies, eg ICD9CM, require strings to be represented accurately: E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      names(id.vocab.code.index)=c("id","vocabulary_id","code","index")
      message("Mapping codes to phecodes...")
      phemapped=mapCodesToPhecodes(id.vocab.code.index, vocabulary.map=vocabulary.map, rollup.map=rollup.map) %>% transmute(id, code=phecode, index)
    }
    
    message("Aggregating codes...")
    phecode=ungroup(summarize(group_by(phemapped,id,code),count=aggregate.fun(index)))
    phecode=phecode[phecode$count>0,]
    
    #Check exclusions, and add them to the list
    if(add.phecode.exclusions) {
      message("Mapping exclusions...")
      exclusions = inner_join(phecode %>% rename(exclusion_criteria=code), phecode_exclude, by = "exclusion_criteria")
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
      #Set exclusions to NA
      phecode$count=ifelse(phecode$count==-1,NA,phecode$count)
    }

    #For min.code.count, use logical fill
    if(!is.na(min.code.count)) { 
      fill=FALSE
    } else {
      #Fill for no min.code.count
      fill=0
    }
    
    message("Reshaping data...")
    phens=spread(phecode,code,count,fill=fill)
    
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
      phens=restrictPhecodesBySex(phens,id.sex)
    }
    
    #Limit to full population ids
    filter(phens, id %in% full.population.ids)
    
    #Rename the ID column to the input ID column name
    names(phens)[1]=id.name
    
    phens
  }