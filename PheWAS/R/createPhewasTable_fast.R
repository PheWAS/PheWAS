createPhewasTable_fast <-
  function(id_icd9_count, min_code_count=2, add_exclusions=T, translate=T, aggregate_fun=sum)
  {
    if(!translate) {
      #Warn about exclusions if input is not translated.
      if(add_exclusions){warning("Codes are not translated, but exclusions are to be applied. Ensure that the icd9 column is phewas codes or disable add_exclusions for accurate results.")}
      phecode=data.frame(id=id_icd9_count[,1],phe=id_icd9_count[,2],count=id_icd9_count[,3])
    } else {
      #check to make sure numeric ICD9 codes were not passed in
      if(class(id_icd9_count[,2]) %in% c("integer","numeric")) {stop("Numeric ICD-9 codes passed in, so an accurate mapping is not possible. E.G.: 250, 250.0, and 250.00 are different codes and necessitate string representation")}
      names(id_icd9_count)=c("id","icd9","count")
      message("Mapping ICD-9 codes to PheWAS codes...")
      phemapped=mapICD9toPheWAS(id_icd9_count)
      phecode=data.frame(id=phemapped[,"id"],phe=phemapped[,"phewas_code"],count=phemapped[,"count"])
    }
    message("Aggregating PheWAS codes...")
    phecode=aggregate(phecode$count, list(phecode$id,phecode$phe), aggregate_fun)
    message("Reshaping matrix...")
    names(phecode)=c("id","phe","count")
    jtab=unique(phecode$phe)
    jtab=data.frame(ja=1:length(jtab),phe=jtab)
    jdim=nrow(jtab)
    itab=unique(phecode$id)
    itab=data.frame(ia=1:length(itab),id=itab)
    idim=nrow(itab)
    phecode=merge(merge(phecode,jtab),itab)
    phens=as.matrix.csr(new("matrix.coo",ra=phecode$count,ja=phecode$ja,ia=phecode$ia,dimension=c(idim,jdim)))
    
    if(add_exclusions) {
      message("Mapping and reshaping exclusions...")
      present_phecode=phecode[phecode$count>0,]
      exclusions=mapPheWAStoExclusions(present_phecode$phe,present_phecode$id)
      exclusions$count=1
      exclusions=exclusions[!duplicated(exclusions[,c("id", "exclusion")]),]
      jtab_ex=jtab
      names(jtab_ex)=c("ja","exclusion")
      exclusions=merge(merge(exclusions,jtab_ex),itab)
      exclusions=as.matrix.csr(new("matrix.coo",ra=exclusions$count,ja=exclusions$ja,ia=exclusions$ia,dimension=c(idim,jdim)))
      exclusion_inds=which(as.logical(SparseM::as.matrix((phens<1)&(exclusions==1))))
      #phens[(phens<1)&(exclusions==1)]=NA
    } else {
      exclusion_inds=c()
    }
    message("Merging exclusions and codes...")
    phens=SparseM::as.matrix(phens)
    phens[exclusion_inds]=NA
    phens=as.data.frame(phens)
    phens=cbind(itab$id,phens)
    names(phens)=c("id",as.character(jtab$phe))
    
    
    if(is.na(min_code_count)) {
      #Do nothing
    } else {
      
      phens=data.frame(phens$id,
                       ifelse(phens[,-1]<min_code_count & phens[,-1]>0,
                              NA,
                              phens[,-1]>=min_code_count))
      names(phens)=c("id",substring(names(phens)[-1],2))
    }
    
      
    phens
  }