#' Fill in later
#'
#' @param phe.gen n/a
#' @param min.records n/a
#' @param return.models n/a
#' @param my.data n/a
#' @param ... n/a
#'
#' @return n/a
#' @export
#' @importFrom lmtest lrtest
#'
#' @examples phe_as_lrt(n/a)
phe_as_lrt <-
  function(phe.gen, min.records=20,return.models=T, my.data, ...) {
    if(!missing(my.data)) data=my.data

    #Retrieve the targets for this loop
    phenotype=phe.gen[[1]]
    gen=phe.gen[[2]]
    gens=paste0(gen,collapse = ", ")
    cov=phe.gen[[3]]
    #Turn covariates into a string, if not NA
    if(!is.na(cov[1])) {covariates=paste(cov,collapse=",")}
    else {covariates=NA_character_} #Make sure it is a character NA for aggregation

    #Subset the data
    d=data[,na.omit(unlist(c(phenotype,gen,cov)))]
    d=na.omit(d)
    #Exclude rows with missing data
    n_total=nrow(d)
    n_cases=NA_integer_
    n_controls=NA_integer_
    type=NA_character_
    note=""
    lrt=NA
    p=NA

    #Drop columns with no variability
    drop.cols = names(d)[sapply(d, function(col) length(unique(col)))<=1]
    if(length(drop.cols>0)) {
      note=paste(note,"[Note: Column(s) dropped due to lack of variability: ",paste0(drop.cols,collapse=", "),"]")
      d=select(d, -one_of(drop.cols))
      #Remove dropped columns from covs- sticks around in the listed "covariates"
      cov=setdiff(cov,drop.cols)
    }

    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(!(phenotype %in% names(d))) {
      note=paste(note,"[Error: non-varying phenotype]")
    } else {

      #Check if phenotype is logical (boolean)
      if(class(d[,phenotype]) %in% c("logical")) {
        type = "logistic - LRT"
        #Create the logistic model
        n_cases=sum(d[,phenotype])
        n_controls=n_total-n_cases
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {
          model = glm(as.formula(paste("`",phenotype,"`"," ~ .", sep="", collapse="")), data=d, family=binomial)
          model.base = glm(as.formula(paste("`",phenotype,"`"," ~ .", sep="", collapse="")), data=d %>% select(-one_of(gen)), family=binomial)
          lrt=lrtest(model, model.base)
          p=lrt$Pr[2]
        }
      } else {
        type = "gaussian"
        if(n_total<min.records) {
          note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
        } else {
          model = lm(as.formula(paste(phenotype," ~ .", sep="", collapse="")), data=d)
          model.base = lm(as.formula(paste(phenotype," ~ .", sep="", collapse="")), data=d %>% select(-one_of(gen)))
          lrt=lrtest(model, model.base)
          p=lrt$Pr[2]
        }
      }
    }
    output=data.frame(phenotype=phenotype,genotype=gens,covariates=covariates,p=p,type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                      note=note, stringsAsFactors=F)



    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$lrt=lrt}
    attributes(output)$successful.phenotype=ifelse(is.na(lrt),NA,phenotype)
    #Return this to the loop to be merged.
    output
  }
