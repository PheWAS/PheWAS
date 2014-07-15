createPhewasTable <-
  function(id.icd9.count, min.code.count=2, add.exclusions=T, translate=T, aggregate.fun=sum, fast=T, id.gender)
  {

    if(fast) {
      if(!require(SparseM)) {
        stop("fast option requires 'SparseM' package")
      }else {
        phens=createPhewasTable_fast(id.icd9.count, min.code.count, add.exclusions, translate, aggregate.fun)
      }
    } else {
    if(!translate) {
      #Warn about exclusions if input is not translated.
      if(add.exclusions){warning("Codes are not translated, but exclusions are to be applied. Ensure that the icd9 column is phewas codes or disable add.exclusions for accurate results.")}
      phecode=data.frame(id=id.icd9.count[,1],phe=id.icd9.count[,2],count=id.icd9.count[,3])
    } else {
      #check to make sure numeric ICD9 codes were not passed in
      if(class(id.icd9.count[,2]) %in% c("integer","numeric")) {stop("Numeric ICD-9 codes passed in, so an accurate mapping is not possible. E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      names(id.icd9.count)=c("id","icd9","count")
      message("Mapping ICD-9 codes to PheWAS codes...")
      phemapped=mapICD9toPheWAS(id.icd9.count)
      phecode=data.frame(id=phemapped[,"id"],phe=phemapped[,"phewas_code"],count=phemapped[,"count"])
    }
    message("Aggregating PheWAS codes...")
    phecode=aggregate(phecode$count, list(phecode$id,phecode$phe), aggregate.fun)
    names(phecode)=c("id","phe","count")
    message("Reshaping matrix...")

    phens=reshape(phecode, direction="wide", idvar="id",timevar="phe")

    phens[is.na(phens)]=0
    
    if(add.exclusions) {
      message("Mapping and reshaping exclusions...")
      present_phecode=phecode[phecode$count>0,]
      exclusions=mapPheWAStoExclusions(present_phecode$phe,present_phecode$id)
      exclusions$count=1
      exclusions=reshape(exclusions, direction="wide", idvar="id",timevar="exclusion")
      exclusions[is.na(exclusions)]=0
      message("Merging exclusions and codes...")
      #Need to map exclusions to phens and then apply exclusions
      both=intersect(names(phens),names(exclusions))[-1]
      both_ids=intersect(phens$id,exclusions$id)
      #Slow call
      phens[phens$id %in% both_ids,both][(phens[phens$id %in% both_ids,both]==0)
                                         &!(exclusions[exclusions$id %in% both_ids,both]==0)]=NA
    }
    
    if(is.na(min.code.count)) {
      #Do nothing
    } else {

      phens=data.frame(phens$id,
                ifelse(phens[,-1]<min.code.count & phens[,-1]>0,
                       NA,
                       phens[,-1]>=min.code.count))
      names(phens)=c("id",substring(names(phens)[-1],7))
    }
    
    
    }
    #If there are gender restrictions, set them to NA
    if(!missing(id.gender)) {
      phens=restrictPhewasByGender(phens,id.gender)
    }
    phens
  }