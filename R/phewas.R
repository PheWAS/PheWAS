#'  Function to perform a PheWAS analysis with multiple methods
#'
#' This function will perform a PheWAS analysis, optionally adjusting for other
#' variables. It is parallelized using the base package parallel.
#'
#'  These results can be directly plotted using the
#'  \code{\link[PheWAS:phewasManhattan]{phewasManhattan}} function, assuming
#'  that models are not returned. If they are, the \code{results} item of the
#'  returned list needs to be used.
#' @param phenotypes  The names of the outcome variables in data under study.
#' These can be logical (for logistic regression) or continuous
#'  (for linear regression) columns.
#' @param genotypes The names of the prediction variables in data under study.
#' @param data Data frame containing all variables for the anaylsis.
#' @param covariates The names of the covariates to appear in every analysis.
#' @param outcomes An alternate to \code{phenotypes}. It will be ignored if
#' \code{phenotypes} exists.
#' @param predictors An alternate to \code{genotypes}. It will be ignored if
#' \code{genotypes} exists.
#' @param cores The number of cores to use in the parallel socket cluster
#' implementation. If \code{cores=1}, \code{lapply} will be used instead.
#' @param additive.genotypes Are additive genotypes being supplied? If so,
#' it will attempt to calculate allele frequencies and HWE values. Default
#' is TRUE.
#' @param method  Determines the statistical method to check associations.
#' One of: 'glm', 'clogit', 'lrt', or 'logistf'.
#'
#' If clogit, requires the \code{strata} parameter to be defined.
#'
#' If lrt, an atomic vector of \code{genotypes} will test all at once. A list
#' of vectors in \code{genotypes} will perform each vector as a test (EG,
#' provide a list of single items to see glm with LRT p-values).
#' @param strata Name of the grouping / strat column necessary for clogit.
#' @param factor.contrasts Contrasts used for factors to generate names used
#' in clogit.
#' @param return.models Return a list the complete models, with the names equal
#' to the string formula used to create them, as well as the results. Default is
#'  FALSE.
#' @param min.records The minimum number of records to perform a test. For
#' logistic regression, there must be at least this number of each cases and
#'  controls, for linear regression this total number of records. Default is 20.
#' @param MASS.confint.level Uses the \code{MASS} package and the \code{confint}
#'  function to calculate a confidence interval at the specified level.
#'   \code{confint} uses a profile likelihood method, which takes some time to
#'    compute. Output is stored in the \code{lower} and \code{upper} columns.
#'    Logistic models will report OR CIs and linear models will report beta CIs.
#'     Default is NA, which does not calculate confidence intervals.
#' @param quick.confint.level Calculate a confidence interval based on
#' \code{beta + or - qnorm * SE}. Output is stored in the \code{lower.q} and
#' \code{upper.q} columns. Logistic models will return have the exponentiated OR
#'  confidence intervals.
#'
#' @return The following are the default rows included in the returned data
#' frame. The attributes of the returned data frame contain additional
#'  information about the anaylsis. If a model did not have sufficient cases
#'  or controls for analysis or failed to converge, NAs will be reported and a
#'  note will be added in the note field.
#' \item{phenotype}{The outcome under study}
#' \item{snp}{The predictor under study}
#' \item{adjustment}{The one off adjustment used}
#' \item{beta}{The beta coefficient for the predictor}
#' \item{SE}{The standard error for the beta coefficient}
#' \item{lower.p}{The lower bound of the quick confidence interval, if requested}
#' \item{upper.p}{The upper bound of the quick confidence interval, if requested}
#' \item{lower}{The lower bound of the \code{confint} confidence interval, if
#' requested}
#' \item{upper}{The upper bound of the \code{confint} confidence interval, if
#' requested}
#' \item{OR}{For logistic regression, the odds ratio for the predictor}
#' \item{p}{The p-value for the predictor}
#' \item{type}{The type of regression model used}
#' \item{n_total}{The total number of records in the analysis}
#' \item{n_cases}{The number of cases in the analysis (logical outcome only)}
#' \item{n_controls}{The number of controls in the analysis (logical outcome only)}
#' \item{HWE_p}{The Hardy-Weinberg equilibrium p-value for the predictor,
#' assuming 0,1,2 allele coding}
#' \item{allele_freq}{The allele frequency in the predictor for the coded allele}
#' \item{n_no_snp}{The number of records with a missing predictor}
#' \item{note}{Additional warning or error information}
#' If there are any requested significance thresholds, boolean variables will be
#' included reporting significance.
#' If \code{return.models=T}, a list is returned. The named item \code{results}
#'  contains the above data frame. The named item \code{models} contains a list
#'  of the models generated in the analysis. To distinguish models, the list is
#'  named by the full formula used in generation.
#' @export
#'
#' @examples
#'  \donttest{
#' #Generate some example data
#'data <- sample_data
#' phenotype_data <- createPhenotypes(data$id.vocab.code.count,
#'  id.sex = data$id.sex) 
#' final_data <- dplyr::inner_join(dplyr::inner_join(data$id.sex, 
#' data$genotypes),  phenotype_data)
#' test_phewas <- phewas_ext(names(phenotype_data)[-1], 
#'                          genotypes = c('rsEXAMPLE'), covariates = 'sex', 
#'                          data = final_data)
#' }
phewas_ext <-
  function(phenotypes, genotypes, data, covariates=NA, outcomes, predictors, cores=1, additive.genotypes=T,
           method="glm", strata=NA, factor.contrasts=contr.phewas,
           return.models=F, min.records=20, MASS.confint.level=NA, quick.confint.level) {
       # print(phenotypes)
    #print(genotypes)
    #print(head(data))
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
    if(!is(covariates, 'list')) { covariates=list(covariates)}

    #Checks for each of the PheWAS methods
    if(method=="glm") {
      association_method=phe_as_ext
     # print('ext')
    } else if (method=="clogit") {
      association_method=phe_as_clogit
      #Check for a strata parameter
      if(is.na(strata)) { stop("clogit requires groups- please provide the column name in the 'strata' parameter.")}
    } else if (method=="lrt") {
      association_method=phe_as_lrt
      #Setup the genotypes as a single complete list if not done already
      if(!is(genotypes, 'list')) {genotypes=list(genotypes)}
    } else if (method == "logistf") {
      association_method=phe_as_logistf
    } else {
      stop("Method must be one of: 'glm', 'clogit', 'lrt', or 'logistf'.")
    }

    para=(cores>1)
    #Create the list of combinations to iterate over
    #print(nrow(phenotypes))
    #print(nrow(genotypes))
   # print(nrow(covariates))
    #print('expand.grid')
    #print(nrow(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F)))
    #print('t')
    #print(nrow(t(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F))))
    #print(nrow(data.frame(t(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F)),stringsAsFactors=F)))
    full_list=data.frame(t(expand.grid(phenotypes,genotypes,covariates,stringsAsFactors=F)),stringsAsFactors=F)
    #print(full_list)
    #If parallel, run the parallel version.
    if(para) {
      #Check to make sure there is no existing phewas cluster.
      if(exists("phewas.cluster.handle")) {
        #If there is, kill it and remove it
        message("Old cluster detected (phewas.cluster.handle), removing...")
        try(parallel::stopCluster(phewas.cluster.handle), silent=T)
        rm(phewas.cluster.handle, envir=.GlobalEnv)
      }
      message("Starting cluster...")
      assign("phewas.cluster.handle", parallel::makeCluster(cores), envir = .GlobalEnv)
      message("Cluster created, finding associations...")
      parallel::clusterExport(phewas.cluster.handle,c("data"), envir=environment())
      parallel::clusterCall(phewas.cluster.handle,library,package="dplyr",character.only=T)
      #Loop across every phenotype- iterate in parallel
      result <-parallel::parLapplyLB(phewas.cluster.handle, full_list, association_method, additive.genotypes=additive.genotypes,
                           confint.level=MASS.confint.level, min.records=min.records,
                           return.models=return.models,factor.contrasts=factor.contrasts,strata=strata)
      #Once we have succeeded, stop the cluster and remove it.
      parallel::stopCluster(phewas.cluster.handle)
      rm(phewas.cluster.handle, envir=.GlobalEnv)
    } else {
      #Otherwise, just use lapply.
      message("Finding associations...")
      #print(association_method)
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
