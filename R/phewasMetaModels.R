phewasMetaModels <- function(results, ...) {
  #Stop the function if the meta package cannot be loaded
  if(!require(meta)) {stop("PheWAS meta functions require the 'meta' package.")}
  results$adjustment=ifelse(is.na(results$adjustment),"NA",results$adjustment)
  out=by(results,INDICES=list(results$phenotype,results$snp,results$adjustment),
         FUN=function(a){
           #Get information about the current step
           phenotype=a$phenotype[1]
           snp=a$snp[1]
           adjustment=a$adjustment[1]
           title=paste0(phenotype," ",snp," ",adjustment)
           #Remove NA results from the input
           a=a[!is.na(a$p),]
           #If there are no records left, return a data frame with a note (so one still knows it was there)
           if(nrow(a)==0) {
             data.frame(note="Error: No succesful analysis in this subset",title=title) 
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
             metagen(TE=a$beta,seTE=a$SE, sm=sm,n.e=a$n_cases,n.c=a$n_controls,
                       studlab=a$study, title=paste0(phenotype," ",snp," ",adjustment), ...)
           }})
  names(out)=sapply(out,function(x){x$title})
  out
}