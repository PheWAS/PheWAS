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
           type=suppressWarnings(max(as.character(a$type),na.rm=TRUE))
           #If type is -Inf, no succesful anaylses exist in this subset
           if(type==-Inf) { data.frame(note="Error: No succesful analysis in this subset",title=title) }
           else { #Otherwise, perform the analysis
           #Throw an error if study types, e.g. logistic and linear, don't match
           if(length(unique(a$type))!=1) stop(paste0("Study types do not match for ",phenotype,", ",snp,", ",adjustment,"."))
           sm=""
           #Define the Odds Ratio summary measure if appropriate
           if(type=="logistic") {
             sm="OR"
           } 
           #Warn if the type is not linear or logistic
           if(!type %in% c("logistic","linear")) {
             warning(paste0("No match for study type ",type,". Assuming no summary measure."))
           }
           metagen(TE=a$beta,seTE=a$SE, sm=sm,n.e=a$n_cases,n.c=a$n_controls,
                   studlab=a$study, title=title,...)}})
  names(out)=sapply(out,function(x){x$title})
  out
}