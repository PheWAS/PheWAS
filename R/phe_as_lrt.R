phe_as_lrt <-
  function(phenotype, min.records=20,return.models=T, my.data, my.covariates, ...) {
    if(!missing(my.data)) data=my.data
    if(!missing(my.covariates)) covariates=my.covariates
    
    #Subset the data
    d=data[,c(phenotype,covariates)]
    d=na.omit(d)
    #Exclude rows with missing data
    n_total=nrow(d)
    n_cases=NA_integer_
    n_controls=NA_integer_
    type=NA_character_
    note=""
    lrt=NA
    p=NA
    if(n_total<min.records) {
      note=paste(note,"[Error: <",min.records," complete records]")
    } else if(length(unique(na.omit(d[,phenotype])))<=1 ) {
      note=paste(note,"[Error: non-varying phenotype]")
    } else {
      
      #Check if phenotype is logical (boolean)
      if(class(d[,phenotype]) %in% c("logical")) {
        type = "logistic"
        #Create the logistic model
        n_cases=sum(d[,phenotype])
        n_controls=n_total-n_cases
        if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
        else {
          model = glm(as.formula(paste("`",phenotype,"`"," ~ .", sep="", collapse="")), data=d, family=binomial)
          model.base = glm(as.formula(paste("`",phenotype,"`"," ~ .", sep="", collapse="")), data=d %>% select(-starts_with("rs")), family=binomial)
          lrt=lrtest(model, model.base)
          p=lrt$Pr[2]
        }
      } else {
        type = "gaussian"
        if(n_total<min.records) {
          note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
        } else {
          model = lm(as.formula(paste(phenotype," ~ .", sep="", collapse="")), data=d)
          model.base = lm(as.formula(paste(phenotype," ~ .", sep="", collapse="")), data=d %>% select(-starts_with("rs")))
          lrt=lrtest(model, model.base)
          p=lrt$Pr[2]
        }
      }
    }
    output=data.frame(phenotype=phenotype,p=p,type=type,
                      n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                      note=note, stringsAsFactors=F)
    
    
    
    #If the complete models were requested, add them as well.
    if(return.models) {attributes(output)$lrt=lrt}
    attributes(output)$successful.phenotype=ifelse(is.na(lrt),NA,phenotype)
    #Return this to the loop to be merged.
    output
  }
