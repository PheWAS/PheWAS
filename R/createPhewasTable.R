createPhewasTable <-
  function(id.icd9.count, min.code.count=2, add.exclusions=T, translate=T, aggregate.fun=sum, id.gender)
  {
    .Deprecated("createPhenotypes", package="PheWAS", msg="Please see PheWAS::createPhenotypes() for more features. This function will be removed in the future.")
    id.name=names(id.icd9.count)[1]
    ids <- unique(id.icd9.count[,1])
    if(!translate) {
      #Warn about exclusions if input is not translated.
      if(add.exclusions){warning("Codes are not translated, but exclusions are to be applied. Ensure that the icd9 column has phecodes or disable add.exclusions for accurate results.")}
      phemapped=tbl_df(data.frame(id=id.icd9.count[,1],phe=id.icd9.count[,2],count=id.icd9.count[,3]))
    } else {
      #check to make sure numeric ICD9 codes were not passed in
      if(class(id.icd9.count[,2]) %in% c("integer","numeric")) {stop("Numeric ICD-9 codes passed in, so an accurate mapping is not possible. E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      names(id.icd9.count)=c("id","icd9","count")
      message("Mapping ICD-9 codes to phecodes...")
      phemapped=mapICD9ToPhecodes(id.icd9.count)
      phemapped=phemapped %>% transmute(id,phe=phecode,count) 
    }
    
    message("Aggregating codes...")
    phecode=ungroup(summarize(group_by(phemapped,id,phe),count=aggregate.fun(count)))
    phecode=phecode[phecode$count>0,]
    
    #Check exclusions, and add them to the list
    if(add.exclusions) {
      message("Mapping exclusions...")
      exclusions=mapPhecodesToExclusions(phecode$phe,phecode$id)
      exclusions$count=-1
      phecode=rbind(phecode,exclusions %>% transmute(id,phe=exclusion,count))
    }
    
    #If there is request for a min code count, adjust counts to -1 if needed
    if(!is.na(min.code.count)&(max(!is.na(phecode$count)&phecode$count<min.code.count))) {
      phecode[!is.na(phecode$count)&phecode$count<min.code.count,]$count=-1
    } 
    
    if(!is.na(min.code.count)|add.exclusions) {
      message("Coalescing exclusions and min.code.count as applicable...")
      phecode=ungroup(summarize(group_by(phecode,id,phe),count=max(count)))
    }
    
    
    #For min.code.count, use logical fill
    if(!is.na(min.code.count)) { 
      fill=FALSE
    } else {
      #Fill for no min.code.count
      fill=0
    }

    message("Reshaping data...")
    phens=spread(phecode,phe,count,fill=fill)
    
    #Set exclusions to NA
    phens[phens==-1]=NA
    
    #Add in inds present in input, but without mapped phecodes
    missing_ids=setdiff(ids,phens[["id"]])
    if(length(missing_ids)>0) {
      empty_record=phens[1,-1]
      empty_record[]=0
      phens=rbind(phens,data.frame(id=missing_ids,empty_record,check.names=F))
    }
    
    #Change to logical if there is a min code count
    if(!is.na(min.code.count)) {phens[,-1]=phens[,-1]>0}


    #If there are gender restrictions, set them to NA
    if(!missing(id.gender)) {
      phens=restrictPhecodesByGender(phens,id.gender)
    }
    
    #Rename the ID column to the input ID column name
    names(phens)[1]=id.name
      
    phens
  }