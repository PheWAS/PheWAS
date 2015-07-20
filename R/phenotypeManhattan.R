phenotypeManhattan <-
  function(d, suggestive.line=0.05, significant.line, 
           OR.size=F,OR.direction=F,
           annotate.level,
           y.axis.interval=5,
           y.axis.label=expression(-log[10](italic(p))),
           max.y,
           ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) stop("Data input must contain columns phenotype and p.")
    if((OR.size|OR.direction)&!length(d$OR)) stop("OR size or direction requested, but d$OR not provided.")
    #Check for size or direction parameter in addition to OR.size or OR.direction
    if(OR.size & !is.null(list(...)$sizes)) stop("You cannot use OR.size=T for OR point sizes and use the second point sizing parameter 'sizes'.")
    if(OR.direction & !is.null(list(...)$direction)) stop("You cannot use OR.direction=T for OR direction shapes and use the second point shape parameter 'direction'.")
    #Check for renamed parameters
    if(max(names(list(...)) %in% c("genomewide.line","sort.by.p","sort.by.category.p"))>0){
      stop("'genomewide.line','sort.by.p', and 'sort.by.category.p' parameters have been replaced with 'significant.line','sort.by.value', and 'sort.by.category.value'")
    }
    #Remove records with NA p values
    d=d[!is.na(d$p),]
    #Transform the significance thresholds
    if(missing(significant.line)) significant.line=suggestive.line/nrow(d)
    significant.line=-log10(significant.line)
    suggestive.line=-log10(suggestive.line)
    if(missing(annotate.level)) {
      if(is.na(significant.line)) {annotate.level=-log10(min(d$p,na.rm=T))}
      else{annotate.level=significant.line}
    }
    else {annotate.level=-log10(annotate.level)}
    #Nudge all of the p=0 results to the smallest double
    if(sum(d$p==0)>0)d[d$p==0,]$p=.Machine$double.xmin
    #Restrict to only those with appropriate p-values
    d=d[d$p>0 & d$p<=1,]
    
    #Create the - log p value for plotting. No errors given restriction to p>0
    d$value = -log10(d$p)
    
    max.y=ifelse(missing(max.y),max(ceiling(max(d$value)),4.4),max.y)
    
    #If OR sizes are requested, normalize them to magnitude only
    #Commented: restrict to those reaching the annotation significance level
    if(OR.size){
      d$size = d$OR
      d[d$size<1,]$size = 1/d[d$size<1,]$size
      #d[d$value>=annotate.level,]$new.OR = 1
    }
    
    #If the OR direction is requested, create it
    if(OR.direction) d$direction = d$OR>=1
    plot=phenotypePlot(d,suggestive.line=suggestive.line,significant.line=significant.line,
                        sizes=OR.size,direction=OR.direction,
                        annotate.level=annotate.level,
                        y.axis.interval=y.axis.interval,
                        y.axis.label=y.axis.label, max.y=max.y,
                        ...)
    if(OR.size) plot=suppressWarnings(plot+scale_size("Odds Ratio", range = c(2, 4), breaks=c(1, 1.2, 1.6)))
    plot
  }