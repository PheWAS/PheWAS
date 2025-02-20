#' Fill in later
#'
#' @param d Data frame containing phenotype and p, or in the case of 
#' phenotypePlot, phenotype and value.
#' @param suggestive.line For phenotypeManhattan, Draw a blue line at this 
#' p-value to show a suggestive significance threshold. Defaults to 0.05. A 
#' missing or NA value as applicable will result in no line. For phenotypePlot, 
#' Blue line drawn at value to emphasize interesting data points. Missing or an
#'  NA value will result in no line.
#' @param significant.line For phenotypeManhattan, draw a red line at this 
#' p-value to show an adjusted significance threshold. Defaults to 
#' suggestive.line divided by the number of non-NA p-values. An NA value will 
#' result in no line. For phenotypePlot, the red line is drawn at the specified 
#' value to emphasize interesting data points. Missing or an NA value will 
#' result in no line.
#' @param OR.size Adjust point size based on odds ratios? Requires d to contain
#'  the OR column. Default is FALSE. Mutually exclusive with sizes.
#' @param OR.direction Adjust point shape based on odds ratio direction? 
#' Requires d to contain the OR column. Default is FALSE.
#' @param sizes n/a
#' @param annotate.level For phenotypeManhattan or phenotypeManhattan, annotate
#'  points at or below this p with a default of significant.line. Set to NA to 
#'  disable adding these annotations. For phenotypePlot, minimum value to 
#'  annotate points. Default of no annotation with a missing arugment, or set 
#'  to NA to disable.
#' @param y.axis.interval The interval for y axis labeling. Defaults to 5. 
#' Scale is in -log10 units for the Manhattan plots and is the same as value for
#'  phenotypePlot.
#' @param y.axis.label The label for the y axis. Defaults to a stylized 
#' "-log10(p)" for the Manhattan plots and "Values" for phenotypePlot.
#' @param max.y For the Manhattan plots, the maximum -log10(p) on the axis.
#'  Defaults to the larger of 4.4 or the maximum -log10(p) in the supplied data. 
#'  For phenotypePlot, the maximum value for the plot. If missing, defaults to 
#'  the maximum value in the supplied data.
#' @param ... n/a
#' @return A ggplot2 plot. It contains PheWAS codes against their -log10 
#' transformed p-values for phewasManhattan, generic phenotypes versus -log10 
#' transformed p-values for phenotypeManhattan, and generic phenotypes versus 
#' generic values for phenotypePlot.

phenotypeManhattan <-
  function(d, suggestive.line=0.05, significant.line,
           OR.size=F,OR.direction=F, sizes,
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

    #Set default of sizing to FALSE
    sizing=FALSE
    #Check sizing parameters
    if(!missing(sizes)) {
      if(sizes==T){
        #If sizes is provided and true, check OR.size and d$size to ensure validity
        if(OR.size) { stop("Both sizes and OR.size are TRUE. Only one size can be used at a time. Please ensure only one parameter is TRUE.") }
        if(!length(d$size)) { stop("sizes were requested, but no d$size was not provided. Please provide the point size scaling variable in d$size") }
        sizing=TRUE
      }
    }
    #If OR sizes are requested, normalize them to magnitude only
    if(OR.size){
      sizing=TRUE
      d$size = d$OR
      d[d$size<1,]$size = 1/d[d$size<1,]$size
    }
    #If the OR direction is requested, create it
    if(OR.direction) d$direction = d$OR>=1
    plot=phenotypePlot(d,suggestive.line=suggestive.line,significant.line=significant.line,
                       sizes=sizing,direction=OR.direction,
                       annotate.level=annotate.level,
                       y.axis.interval=y.axis.interval,
                       y.axis.label=y.axis.label, max.y=max.y,
                       ...)
    if(OR.size) plot=suppressMessages(plot+scale_size("Odds Ratio", range = c(2, 4), breaks=c(1, 1.2, 1.6)))
    plot
  }