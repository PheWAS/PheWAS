phewasMeta <- function(results, fixed=T, keep.both=T, ...) {
  #Stop the function if the meta package cannot be loaded
  if(!require(meta)) {stop("PheWAS meta functions require the 'meta' package.")}
  #Replace NA adjustment values with "NA"- necessary for by
  results$adjustment=ifelse(is.na(results$adjustment),"NA",results$adjustment)
  #Iterate across all phenotype, snp, adjustment combinations.
  out=by(results,INDICES=list(results$phenotype,results$snp,results$adjustment),
         FUN=function(a){
           #Get information about the current step
           phenotype=a$phenotype[1]
           snp=a$snp[1]
           adjustment=a$adjustment[1]
           #Remove NA results from the input
           a=a[!is.na(a$p),]
           #If there are no records left, return an NA row (so one still knows it was there)
           if(nrow(a)==0) {
             #create an NA row to return.
             c=data.frame(phenotype=phenotype,snp=snp,adjustment=adjustment,
                          beta=NA_real_, OR=NA_real_,SE=NA_real_,p=NA_real_,type=NA_character_,
                          n_total=sum(a$n_total,na.rm=TRUE),
                          n_cases=sum(a$n_cases,na.rm=TRUE),
                          n_controls=sum(a$n_controls,na.rm=TRUE),
                          HWE_p.min=NA_real_,allele_freq=NA_real_,
                          n_no_snp=NA_integer_,k_studies=0,tau2=NA_real_,I2.percent=NA_real_,
                          Q=NA_real_,Q.df=NA_real_,Q.p=NA_real_,beta.fixed=NA_real_,OR.fixed=NA_real_,
                          SE.fixed=NA_real_,p.fixed=NA_real_,beta.random=NA_real_,
                          OR.random=NA_real_,SE.random=NA_real_,p.random=NA_real_)
           } else {#If there was at least one good analysis, calculate the meta-analysis
             #Define the Odds Ratio summary measure if appropriate
             sm=""
             type=a$type[1]
             if(type=="logistic") {
               sm="OR"
             } 
             #Warn if the type is not linear or logistic
             if(!type %in% c("logistic","linear")) {
               warning(paste0("No match for study type ",type,". Assuming no summary measure for metagen."))
             }
             #Throw an error if study types, e.g. logistic and linear, don't match
             if(length(na.omit(unique(a$type)))>1) stop(paste0("Study types do not match for ",phenotype,", ",snp,", ",adjustment,"."))
             #Perform the meta analysis.
             b=metagen(TE=a$beta,seTE=a$SE, sm=sm,n.e=a$n_cases,n.c=a$n_controls,
                       studlab=a$study, title=paste0(phenotype," ",snp," ",adjustment), ...)
             #Calculate I2, setting to 0 if below 0
             I2=(b$Q-b$df.Q)/b$Q*100
             I2=ifelse(I2>=0,I2,0)
             #Calculate the total N
             n_total=sum(b$n.e,b$n.c)
             #Create the data frame with the meta analysis results
             c=data.frame(phenotype=phenotype,snp=snp,adjustment=adjustment,
                          beta=NA_real_,
                          OR=NA_real_,
                          SE=NA_real_,
                          p=NA_real_,
                          type=type,
                          n_total=n_total,
                          n_cases=sum(b$n.e),
                          n_controls=sum(b$n.c),
                          HWE_p.min=min(a$HWE_p),
                          allele_freq=sum(a$allele_freq*a$n_total)/n_total,
                          n_no_snp=sum(a$n_no_snp),
                          k_studies=b$k,
                          tau2=b$tau^2,
                          I2.percent=I2,
                          Q=b$Q,
                          Q.df=b$df.Q,
                          Q.p=pchisq(b$Q,b$df.Q,lower.tail=FALSE),
                          beta.fixed=b$TE.fixed,
                          OR.fixed=exp(b$TE.fixed),
                          SE.fixed=b$seTE.fixed,
                          p.fixed=b$pval.fixed,
                          beta.random=b$TE.random,
                          OR.random=exp(b$TE.random),
                          SE.random=b$seTE.random,
                          p.random=b$pval.random
             )
             
             #Set OR measures to NA if the summary measure is not Odds Ratio
             if(sm!="OR") {
               c$OR.fixed=NA_real_
               c$OR.random=NA_real_
             }
           }
           #Set the default attributes based on user preference
           if(fixed) {
             c$beta=c$beta.fixed
             c$OR=c$OR.fixed
             c$SE=c$SE.fixed
             c$p=c$p.fixed
           } else {
             c$beta=c$beta.random
             c$OR=c$OR.random
             c$SE=c$SE.random
             c$p=c$p.random
           }
           #remove the fixed/random specific attributes if requested.
           if(!keep.both) {
             c=c[,-21:-28]
           }
           c
         })
  rbind_all(out)
}