phe_as_ext <-
function(phe.gen, additive.genotypes=T,min.records=20,return.models=F,confint.level=NA, my.data) {
  if(!missing(my.data)) data=my.data
  #Retrieve the targets for this loop
  phe_o=phe.gen[[1]]
  phe=phe_o
  gen=phe.gen[[2]]
  gens=gen
  covariates=phe.gen[[3]]
  #Subset the data
  d=data %>% select(one_of(na.omit(unlist(c(phe,gen,covariates)))))
  #Turn covariates into a string, if not NA
  if(!is.na(covariates[1])) {covariates=paste(covariates,collapse=",")}
  else {covariates=NA_character_} #Make sure it is a character NA for aggregation
  #Alter phe_o if necessary for the regression formula
  if(suppressWarnings(!is.na(as.numeric(phe_o)))) {
    phe=paste0("pheno_",phe_o)
    names(d)[1]=phe
  }
  #Exclude the exclusions for the target phenotype
  d=d[!is.na(d[[phe]]),]
  n_no_snp=sum(is.na(d[[gen]]))
  #Exclude rows with missing data
  d=na.omit(d)
  n_total=nrow(d)
  n_cases=NA_integer_
  n_controls=NA_integer_
  allele_freq=NA_real_
  HWE_pval=NA_real_
  or=NA_real_
  se=NA_real_
  p=NA_real_
  beta=NA_real_
  type=NA_character_
  note=""
  model=NA
  if(n_total<min.records) {
    note=paste(note,"[Error: <",min.records," complete records]")
  } else if(min(sapply(d[,c(phe,gen)],FUN=function(x){length(unique(x))}))<=1) {
    note=paste(note,"[Error: non-varying phenotype or genotype]")
  } else {
    if(additive.genotypes) {
      allele_freq=sapply(d %>% select(one_of(gen))),
                    FUN=function(x){
                      if(class(x) %in% c("numeric","integer")){
                        sum(x)/(2*n_total)
                      } else {NA_real_}
                    })
      HWE_pval=sapply(d %>% select(one_of(gen))),
            FUN=function(x){
              if(class(x) %in% c("numeric","integer") & sum(!(na.omit(x) %in% 0:2))==0) {
                P=allele_freq
                Q=1-allele_freq
                AA=sum(x==2)
                xAA=P^2*n_total
                Aa=sum(x==1)
                xAa=2*P*Q*n_total
                aa=sum(x==0)
                xaa=Q^2*n_total
                pchisq((AA-xAA)^2/(xAA)+(Aa-xAa)^2/(xAa)+(aa-xaa)^2/(xaa),1)
              } else {NA_real_}
            }
      )
      #Report a warning as needed.
      if() {note=paste(note,"[Warning: At least one genotype was not coded 0,1,2, but additive.genotypes was TRUE.]")}
    } 
    #Check if genotype was available
    #Check if phenotype is logical (boolean)
    if(class(d[[phe]]) %in% c("logical")) {
      type = "logistic"
      #Create the logistic model
      n_cases=sum(d[[phe]])
      n_controls=n_total-n_cases
      if(n_cases<min.records|n_controls<min.records) {note=paste(note,"[Error: <",min.records," cases or controls]")}
      else {
        model = glm(as.formula(paste(phe," ~ .")), data=d, family=binomial)
        modsum= summary(model)
        #If the models did not converge, report NA values instead.
        if(model$converged) {
          #Find the rows with results that gets merged across all loops
          gen_list=grep(gen,row.names(modsum$coef))
          gens=row.names(modsum$coef)[gen_list]
          or=exp(modsum$coef[gen_list,1])
          beta=modsum$coef[gen_list,1]
          se=modsum$coef[gen_list,2]
          p=modsum$coef[gen_list,4]
        } else {
          note=paste(note,"[Error: The model did not converge]")
        }
      }
    } else {
      type = "linear"
      if(n_total<min.records) {
        note=paste(note,"[Error: <",min.records," records with phenotype and genotype]")
      } else {
        model = lm(as.formula(paste(phe," ~ .", sep="", collapse="")), data=d)
  
        modsum= summary(model)
        #If the models did not converge, report NA values instead.
        if(model$converged) {
          #Find the rows with results that gets merged across all loops
          gen_list=grep(gen,row.names(modsum$coef))
          gens=row.names(modsum$coef)[gen_list]
          beta=modsum$coef[gen_list,1]
          se=modsum$coef[gen_list,2]
          p=modsum$coef[gen_list,4]
        } else {
          note=paste(note,"[Error: The model did not converge]")
        }
      }
    }
  }
  
  #Check to see there were numbers (ie phewas codes) as the predictor, and clean up.
  if(suppressWarnings(!is.na(as.numeric(gen)))){
    gens=substring(gens,2)
    #Clean up if there was just one TRUE entry
    if(length(gens)==1 & substring(gens,length(gens)-4)=="TRUE") {
      gens=substring(gens,1,length(gens)-4)
    }
  }
  output=data.frame(phenotype=phe_o,snp=gens,
                    covariates=covariates,
                    beta=beta, SE=se,
                    OR=or,
                    p=p, type=type,
                    n_total=n_total, n_cases=n_cases, n_controls=n_controls,
                    HWE_p=HWE_pval,allele_freq=allele_freq,n_no_snp=n_no_snp, 
                    note=note, stringsAsFactors=F)

  #Add confidence intervals if requested.
  if(!is.na(confint.level)) {
    if(!is.na(model)[1]){
      suppressMessages(conf<-confint(model,c("(Intercept)",gens),level=confint.level))
      lower=conf[-1,1]
      upper=conf[-1,2]
      if(type=="logistic") {
        lower=exp(lower)
        upper=exp(upper)
      }
    } else {
      lower=NA_real_
      upper=NA_real_
    }
    output$lower=lower
    output$upper=upper
    
    output=output[,c("phenotype","snp","covariates","beta","SE",
                                "lower","upper","OR","p","type",
                                "n_total","n_cases","n_controls",
                                "HWE_p","allele_freq","n_no_snp","note")]
  }
  
  #If the complete models were requested, add them as well.
  if(return.models) {attributes(output)$model=model}
  attributes(output)$successful.phenotype=ifelse(is.na(p),NA,phe_o)
  attributes(output)$successful.genotype=ifelse(is.na(p),NA,gen)
  #Return this to the loop to be merged.
  output
}
