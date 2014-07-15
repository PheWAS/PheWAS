phewas <-
function(phenotypes,genotypes,data,covariates=c(NA),adjustments=list(NA), outcomes, predictors, cores=1, additive.genotypes=T, 
         significance.threshold, alpha=0.05, unadjusted=F, return.models=F, min.records=20) {
  if(missing(phenotypes)) {
    if(!missing(outcomes)) phenotypes=outcomes
    else stop("Either phenotypes or outcomes must be passed in.")
  }
  if(missing(genotypes)) {
    if(!missing(predictors)) genotypes=predictors
    else stop("Either genotypes or predictors must be passed in.")
  }
  association_method=phe_as
  if(unadjusted) {
    association_method=phe_as_unadjusted
    if(!is.na(covariates) | !is.na(adjustments)) warning("Covariates and adjustments are ignored in unadjusted mode.")
  }
  #If input is formatted as a set of data frames, reshape it
  if(missing(data)) {
    phe=phenotypes
    gen=genotypes
    cov=covariates
    adjustment=adjustments
    id=intersect(names(phenotypes),names(genotypes))
    if(length(id)==0) {stop("There is no shared column to merge phenotypes and genotypes!")}
    message(paste("Merging data using these shared columns: ",id))
    phenotypes=names(phenotypes)
    phenotypes=phenotypes[!(phenotypes %in% id)]
    genotypes=names(genotypes)
    genotypes=genotypes[!(genotypes %in% id)]
    data=merge(phe,gen,by=id)
    if(!is.null(names(covariates)[-1])) {
      covariates=names(covariates)
      if(sum(id %in% covariates)!=length(id)){stop(paste("The shared ID column(s) do not all exist in covariates: ",id))}
      covariates=covariates[!(covariates %in% id)]
      data=merge(data,cov,by=id)
    }
    if(!is.null(names(adjustments)[-1])) {
      adjustments=names(adjustments)
      if(sum(id %in% adjustments)!=length(id)){stop(paste("The shared ID column(s) do not all exist in adjustments: ",id))}
      adjustments=as.list(c(NA,adjustments[!(adjustments %in% id)]))
      data=merge(data,adjustment,by=id)
    }
  }
  para=(cores>1)
  #Create the list of combinations to iterate over
  full_list=data.frame(t(expand.grid(phenotypes,genotypes,adjustments,stringsAsFactors=F)),stringsAsFactors=F)
  message("Finding associations...")
  #If parallel, run the snowfall version.
  if(para) {
    require(snowfall)
    sfInit(parallel=para, cpus=cores)
    sfExport("data", "covariates")
    #Loop across every phenotype- iterate in parallel
    result <- sfClusterApplyLB(full_list, phe_as, additive.genotypes, min.records,return.models)
    sfStop()
  } else {
    #Otherwise, just use lapply.
    result=lapply(full_list,FUN=phe_as, additive.genotypes, min.records,return.models, data, covariates)
  }
  if(return.models) {
    models=lapply(result,function(x){attributes(x)$model})
    names(models)=sapply(models,function(x){paste0(as.character(terms(x))[c(2,1,3)],collapse=" ")})
  }
  successful.phenotypes=na.omit(sapply(result,function(x){attributes(x)$successful.phenotype}))
  n.tests=length(successful.phenotypes)
  successful.phenotypes=unique(successful.phenotypes)
  successful.genotypes=unique(na.omit(sapply(result,function(x){attributes(x)$successful.genotype})))
  sig=rbind.fill(result)
  
  
  #Add significance thresholds
  attributes(sig)$alpha=alpha
  attributes(sig)$n.tests=n.tests
  if(!missing(significance.threshold)) {
    message("Finding significance thresholds...")
    thresh=match(c("p-value","bonferroni","fdr","simplem-genotype","simplem-phenotype","simplem-product"),significance.threshold)
    sm.g=1
    sm.p=1
    #p.value
    if(!is.na(thresh[1])) {
      sig$p.value=sig$p<=alpha
    }
    #bonferroni
    if(!is.na(thresh[2])) {
      sig$bonferroni=sig$p<=alpha/n.tests
      attributes(sig)$bonferroni=alpha/n.tests
    }
    #fdr
    if(!is.na(thresh[3])) {
      sig$fdr=p.adjust(sig$p,method="fdr")<=alpha
    }
    #simplem-genotype
    if(!is.na(thresh[4])|!is.na(thresh[6])) {
      if(length(successful.genotypes)>1){
        eigs=eigen(cor(data[,genotypes],use="pairwise.complete.obs",method="spearman"))[[1]]
        max.eig=sum(eigs)
        sm.g=which.max(cumsum(eigs)>.995*max.eig)
      } else {
        sm.g=1
      }
      sig$simplem.genotype=sig$p<=alpha/sm.g
      attributes(sig)$simplem.genotype=alpha/sm.g
      attributes(sig)$simplem.genotype.meff=sm.g
    }
    #simplem-phenotype
    if(!is.na(thresh[5])|!is.na(thresh[6])) {
      if(length(successful.phenotypes>1)) {
        eigs=try(cor(data[,successful.phenotypes],use="pairwise.complete.obs",method="spearman"),silent=T)
        if(class(eigs)!="try-error"){
          eigs[is.na(eigs)]=0
          eigs=eigen(eigs)[[1]]
          max.eig=sum(eigs)
          sm.p=which.max(cumsum(eigs)>.995*max.eig)
        } else {
          warning("Phentoype correlation generation failed; this is typically due to sparse phenotypes.")
          sm.p=NA
        }
      } else {
        sm.p=1
      }
      sig$simplem.phenotype=sig$p<=alpha/sm.p
      attributes(sig)$simplem.phenotype=alpha/sm.p    
      attributes(sig)$simplem.phenotype.meff=sm.p
    }
    #simplem-product
    if(!is.na(thresh[6])) {
      sm=sm.g*sm.p
      sig$simplem.product=sig$p<=alpha/sm
      attributes(sig)$simplem.product=alpha/sm
      attributes(sig)$simplem.product.meff=sm
    }
  }
  if(!missing(outcomes)) names(sig)[names(sig)=="phenotype"]="outcome"
  if(!missing(predictors)) names(sig)[names(sig)=="snp"]="predictor"
  if(return.models){sig=list(results=sig,models=models)}
  return(sig)
}
