phe_as_clogit <-
  function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, factor.contrasts=NA, strata, my.data, ...) {
    if(!missing(my.data)) data=my.data
    #Retrieve the targets for this loop
    phe=phe.gen[[1]]
    gen=phe.gen[[2]]
    gens=gen
    cov=phe.gen[[3]]
    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation
    
    #Subset the data
    d=data %>% select(one_of(na.omit(unlist(c(phe,gen,cov,strata)))))
    
    #Exclude the exclusions for the target phenotype
    d=d[!is.na(d[[phe]]),]
    n_no_snp=sapply(d %>% select(one_of(gen)),
                    FUN=function(x){sum(is.na(x))})
    #Exclude rows with missing data
    d=na.omit(d)
    #Exclude strata with only one predictor class (or value) represented
    strata.keep=d %>% group_by_at(.vars=strata) %>% summarize_at(.funs="n_distinct",na.rm=TRUE,.vars=gen) %>% 
      filter_at(.vars=gen,.vars_predicate=all_vars(.>1)) %>% select(one_of(strata))
    d = inner_join(d,strata.keep,by=strata)
    
    n_total=nrow(d)
    n_cases=NA_integer_
    n_controls=NA_integer_
    n_strata=NA_integer_
    allele_freq=NA_real_
    HWE_pval=NA_real_
    or=NA_real_
    se=NA_real_
    p=NA_real_
    beta=NA_real_
    type=NA_character_
    note=""
    model=NA
    formula.string=NA_character_
    expanded_formula=NA_character_
    gen_expansion=1:length(gen)
    
    #Drop columns with no variability
    drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
    if(length(drop.cols>0)) {
      note=paste(note,"[Note: Column(s) dropped due to lack of variability: ",paste0(drop.cols,collapse=", "),"]")
      d=select(d, -one_of(drop.cols))
      #Remove dropped columns from covs- sticks around in the listed "covariates"
      cov=setdiff(cov,drop.cols)
    }

    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation
    
    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(sum(c(phe,gen) %in% names(d))!=length(c(phe,gen))) {
      note=paste(note,"[Error: non-varying phenotype or genotype]")
    } else {
      if(additive.genotypes) {
        snp.details=lapply(d %>% select(one_of(gen)),
                           FUN=function(x){
                             if(class(x) %in% c("numeric","integer")){
                               a.f=sum(x)/(2*n_total)
                               if(sum(!(na.omit(x) %in% 0:2))==0) {
                                 P=a.f
                                 Q=1-a.f
                                 AA=sum(x==2)
                                 xAA=P^2*n_total
                                 Aa=sum(x==1)
                                 xAa=2*P*Q*n_total
                                 aa=sum(x==0)
                                 xaa=Q^2*n_total
                                 hwe=pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
                               } else {hwe=NA_real_}
                             } else {
                               a.f=NA_real_
                               hwe=NA_real_
                             }
                             c(a.f,hwe)
                           })
        
        
        allele_freq=sapply(snp.details,FUN=`[`,1)
        HWE_pval=sapply(snp.details,FUN=`[`,2)
        
        #Report a warning as needed.
        if(min(sapply(d %>% select(one_of(gen)),
                      FUN=function(x){class(x) %in% c("numeric","integer") & sum(!(na.omit(x) %in% 0:2))==0}))!=TRUE){
          note=paste(note,"[Warning: At least one genotype was not coded 0,1,2, but additive.genotypes was TRUE.]")}
      } 
      #Alter factors to use special contrasts
      if(suppressWarnings(!is.na(factor.contrasts))) {
        d=data.frame(lapply(d,function(x){
          if("factor" %in% class(x)){
            x=droplevels(x)
            contrasts(x)=factor.contrasts(x)
          }
          x}),check.names=F)
      }
      #Create the formula:
      formula.string=paste0("`",phe,"` ~ `",paste(c(gen,cov),collapse = "` + `"),'`'," + strata(`",strata,"`)")
      my.formula = as.formula(formula.string)
      
      #Check if phenotype is logical (boolean)
      if(class(d[[phe]]) %in% c("logical")) {
        type = "conditional logistic"
        #Create the logistic model
        n_cases=sum(d[[phe]])
        n_controls=n_total-n_cases
        n_strata=n_distinct(d[[strata]])
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {
          
          model = tryCatch(clogit(my.formula, data=d), warning = function(w) {w$message}, error = function(e) {e$message})
          #If the models did not converge, report NA values instead.
          if(class(model)[1]!="character") {
            #Find the observed genotype columns
            gen_expansion=attr(model.matrix(my.formula, data=d),"assign")
            gen_list=which(gen_expansion %in% 1:length(gen))
            gen_expansion=gen_expansion[gen_list]
            gen_list=gen_list-1
            #Create the model
            modsum= summary(model)
            #Find the rows with results that gets merged across all loops
            gens=row.names(modsum$coefficients)[gen_list]
            or=modsum$coefficients[gen_list,2]
            beta=modsum$coefficients[gen_list,1]
            se=modsum$coefficients[gen_list,3]
            p=modsum$coefficients[gen_list,5]
            expanded_formula=paste0(c(names(model$coefficients),paste0("strata(`","strata","`)")),collapse=" + ")
          } else {
            note=paste0(note,"[Error: Potential fitting problem (try changing covariates). clogit error: ",model,"]")
            model=NA
          }
        }
      } else {
        type = "linear"
        note=paste(note,"[Error: clogit requires a logical/Boolean outcome]")
      }
    }
    
    output=data.frame(phenotype=phe,snp=gens,
                      covariates=covariates,
                      beta=beta, SE=se,
                      OR=or,
                      p=p, type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,n_strata=n_strata,
                      HWE_p=HWE_pval[gen_expansion],allele_freq=allele_freq[gen_expansion],n_no_snp=n_no_snp[gen_expansion], 
                      formula=formula.string,
                      expanded_formula = expanded_formula,
                      note=note, stringsAsFactors=F)
    
    #Add confidence intervals if requested.
    if(!is.na(confint.level)) {
      if(!is.na(model)[1]){
        modsum= summary(model,conf.int=confint.level)
        lower=modsum$conf.int[gen_list,3]
        upper=modsum$conf.int[gen_list,4]
      } else {
        lower=NA_real_
        upper=NA_real_
      }
      output$lower=lower
      output$upper=upper
      
      output=output[,c("phenotype","snp", "covariates","beta","SE",
                       "lower","upper","OR","p","type",
                       "n_total","n_cases","n_controls","n_strata",
                       "HWE_p","allele_freq","n_no_snp","formula","expanded_formula","note")]
    }
    
    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$model=model}
    attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe)
    attributes(output)$successful.genotype=ifelse(is.na(p),NA,gen)
    #Return this to the loop to be merged.
    output
  }
