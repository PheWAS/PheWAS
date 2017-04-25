phewas_ext <-
function(phenotypes, genotypes, data, covariates=NA, outcomes, predictors, cores=1, additive.genotypes=T,
         method="glm", strata=NA, factor.contrasts=contr.phewas,
         return.models=F, min.records=20, MASS.confint.level=NA, quick.confint.level) {
  #Require an input data frame
  if(missing(data)) {
    stop("A data frame must be supplied in 'data'")
  }
  #Rename outcomes and predictors parameters if used.
  if(missing(phenotypes)) {
    if(!missing(outcomes)) phenotypes=outcomes
    else stop("Either phenotypes or outcomes must be passed in.")
  }
  if(missing(genotypes)) {
    if(!missing(predictors)) genotypes=predictors
    else stop("Either genotypes or predictors must be passed in.")
  }
  #Convert covariates to a list if it is not one
  if(class(covariates)!="list") { covariates=list(covariates)}
  
  #Checks for each of the PheWAS methods
  if(method=="glm") {
    association_method=phe_as_ext
  } else if (method=="clogit") {
    association_method=phe_as_clogit
    #Check for a strata parameter
    if(is.na(strata)) { stop("clogit requires groups- please provide the column name in the 'strata' parameter.")}
  } else if (method=="lrt") {
    association_method=phe_as_lrt
    #Setup the genotypes as a single complete list if not done already
    if(class(genotypes)!="list") {genotypes=list(genotypes)}
  } else if (method == "logistf") {
    association_method=phe_as_logistf
  } else if (method == "logistf_multips") {
    association_method=phe_as_logistf_multips
  } else {
    stop("Method must be one of: 'glm', 'clogit', 'lrt', or 'logistf'.")
  }
  
  para=(cores>1)
  #Create the list of combinations to iterate over
  full_list=data.frame(t(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F)),stringsAsFactors=F)

  #If parallel, run the parallel version.
  if(para) {
    #Check to make sure there is no existing phewas cluster.
    if(exists("phewas.cluster.handle")) {
      #If there is, kill it and remove it
      message("Old cluster detected (phewas.cluster.handle), removing...")
      try(stopCluster(phewas.cluster.handle), silent=T)
      rm(phewas.cluster.handle, envir=.GlobalEnv)
    }
    message("Starting cluster...")
    assign("phewas.cluster.handle", makeCluster(cores), envir = .GlobalEnv)
    message("Cluster created, finding associations...")
    clusterExport(phewas.cluster.handle,c("data"), envir=environment())
    clusterCall(phewas.cluster.handle,library,package="dplyr",character.only=T)
    #Loop across every phenotype- iterate in parallel
    result <-parLapplyLB(phewas.cluster.handle, full_list, association_method, additive.genotypes=additive.genotypes, 
                         confint.level=MASS.confint.level, min.records=min.records,
                         return.models=return.models,factor.contrasts=factor.contrasts,strata=strata)
    #Once we have succeeded, stop the cluster and remove it.
    stopCluster(phewas.cluster.handle)
    rm(phewas.cluster.handle, envir=.GlobalEnv)
  } else {
    #Otherwise, just use lapply.
    message("Finding associations...")
    result=lapply(full_list,FUN=association_method, additive.genotypes=additive.genotypes, 
                  confint.level=MASS.confint.level, my.data=data, min.records=min.records,
                  return.models=return.models,factor.contrasts=factor.contrasts,strata=strata)
  }

  if(return.models) {
    message("Collecting models...")
    models=lapply(result,function(x){attributes(x)$model})
    names(models)=sapply(models,function(x){paste0(as.character(terms(x))[c(2,1,3)],collapse=" ")})
  }
  
  message("Compiling results...")
  successful.phenotypes=na.omit(sapply(result,function(x){attributes(x)$successful.phenotype}))
  n.tests=length(successful.phenotypes)
  successful.phenotypes=unique(successful.phenotypes)
  successful.genotypes=unique(na.omit(sapply(result,function(x){attributes(x)$successful.genotype})))
  sig=bind_rows(result)
  
  #Report warning if any convergence errors
  if(max(grepl(pattern = "[Error: The model did not converge]", sig$note, fixed=TRUE))){
    warning("Not all models converged, check the notes column for details.")
  }
  
  message("Cleaning up...")
  
  if(!missing(outcomes)) names(sig)[names(sig)=="phenotype"]="outcome"
  if(!missing(predictors)) names(sig)[names(sig)=="snp"]="predictor"
  if(return.models){sig=list(results=sig,models=models)}
  if(!missing(quick.confint.level)) {
    if(quick.confint.level>=1|quick.confint.level<=0) {warning("Quick confidence interval requested, but a value in the range (0,1) was not supplied")}
    else {
      sig.names=names(sig)
      two.sided=(1-quick.confint.level)/2
      sig=sig %>% mutate(lower.q=beta+se*qnorm(two.sided),upper.q=beta+se*qnorm(two.sided,lower.tail=F))
      sig=sig  %>% mutate(lower.q=ifelse(sig$type=="logistic",exp(lower.q),lower.q),
                      upper.q=ifelse(sig$type=="logistic",exp(upper.q),upper.q))
      sig=sig[,c(sig.names[1:5],"lower.q","upper.q",sig.names[6:length(sig.names)])]
    }
  }
  return(sig)
}
