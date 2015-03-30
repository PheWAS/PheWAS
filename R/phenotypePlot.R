phenotypePlot <-
  function(d, max.y,max.x, suggestive.line, significant.line, 
           size.x.labels=9, size.y.labels=9, switch.axis=F, sort.by.value=F, sort.by.category.value=F,
           #annotation
           annotate.phenotype.description,
           annotate.angle=30, annotate.size=5, annotate.level,
           annotate.phenotype=F, 
           annotate.snp.w.phenotype=F,
           annotate.snp=F, annotate.snp.angle=0,
           annotate.list, annotate.only.largest=T,
           #labels
           lc.labels=F,
           x.group.labels=T, x.phenotype.labels=F,
           sizes=F, direction=F, point.size=3, 
           #plot characteristics
           use.color=T,
           color.palette,
           title= paste0("Phenotype Plot ", date()),
           x.axis.label="Phenotypes",
           y.axis.label="Values",
           y.axis.interval=10) {
    #Check to ensure the input contains columns phenotype and value
    if ( sum(c("phenotype","value") %in% names(d))<2 ) stop("Data input must contain columns phenotype and value.")
    
    #Check for annotation information
    if (!missing(annotate.phenotype.description)) {
      if(class(annotate.phenotype.description)=="data.frame" & sum(c("phenotype","description") %in% names(annotate.phenotype.description))==2) {
        #Add annotation
        d=merge(d,annotate.phenotype.description)
        annotate.phenotype.description=T
      } else if(annotate.phenotype.description==T & length(d$description)) {
        #Do nothing, as it is ready.
      } else {
        annotate.phenotype.description=F
        stop("Annotate.phenotype must contain columns phenotype and description, or be TRUE with provided d$description.")
      }
    } else {
      annotate.phenotype.description=F
    }
    if ((annotate.snp|annotate.snp.w.phenotype) & !("snp" %in% names(d))) stop("You requested SNP annotation but d$snp is not defined.")
    if(annotate.snp.w.phenotype&annotate.snp) warning("You requested SNP annotation twice")
    
    #Is there any annotation?
    annotate=annotate.phenotype.description|annotate.phenotype|annotate.snp|annotate.snp.w.phenotype
    #If so, but nothing listing to annotate, warn and remove annotate flag
    if(annotate&missing(annotate.level)&missing(annotate.list)) {
      warning("You requested annotation, but did not specify annotate.level or annotate.list.")
      annotate=F
    }
    #Check for conflicting color commands
    if(!use.color & !missing(color.palette)) stop("You requested no color, but provided a color palette.")
    #Check for color
    if(use.color) {
      if(missing(color.palette)&!("color" %in% names(d))) stop("You requested color, but did not provide a color attribute in d or color.palette")
      #Set up color if using a palette
      if(!missing(color.palette)&!length(d$color)) {
        if(length(d$groupnum)) d$color=d$groupnum
        else stop("You requested use.color, but d$color or d$groupnum were not provided to distinguish groups.")
      }
    } else {
      #Set the colors to all black if no color is requested
      d$color="#000000"
    }
    
    #Check for point sizing/direction information if requested
    if(sizes&!length(d$size)) stop("You requested size information, but did not provide d$size")
    if(direction&!length(d$direction)) stop("You requested direction information, but did not provide d$direction")
    #Set point size if sizes are not requested
    if(!sizes) d$size=point.size
    
    #Remove lc.labels flag if no annotations
    if(!annotate.phenotype.description&!annotate.phenotype&!annotate.snp&!annotate.snp.w.phenotype) lc.labels=F
    
    #One cannot sort by value and have x.group.labels
    if (sort.by.value & x.group.labels) stop("One cannot sort universally by value and have x group labels. Try sort.by.category.value or not labeling the groups.")
    #Check for group information if requested
    if((sort.by.category.value|x.group.labels)&(sum(c("groupnum","group") %in% names(d))<2)) {
      stop("Requested group information, but did not provide d$groupnum and d$group.")
    }
    
    #If no group info, just assign all to group 0
    if(!("groupnum"%in%names(d))) { 
      d$groupnum=0
    }
    
    #Remove items with NA values or phenotypes
    d=d[!(is.na(d$phenotype)|is.na(d$value)),]
    
    #Sort by the phenotype
    d=d[order(d$phenotype),]
        
    #Set the maximum x value to fit all phenotypes if not specified
    if(missing(max.x)) max.x = length(unique(d$phenotype))
    
    
    #Create the list of phenotypes, finding the best values for each phenotype
    phenotypes=aggregate(value ~ phenotype + groupnum, d,FUN=max)
    #Remove the least significant phenotypes; only has an effect if max.x was specified.
    phenotypes=phenotypes[order(phenotypes$value, decreasing=T),][1:max.x,]
    
    #If the user requested sorting by values
    if (sort.by.value) {
      #Sort by values
      phenotypes=phenotypes[order(phenotypes$value, decreasing=T),]
    }  else if (sort.by.category.value) {
      #If the user requested soring by values within each group.
      #Sort by group and then value
      phenotypes=phenotypes[order(-phenotypes$groupnum,phenotypes$value, decreasing=T),]
    } else {
      #Restore the phenotype sorting order if no other sorting was requested
      #TODO: What about non-numeric phenotypes? Need this to sort 008 etc.
      phenotypes=phenotypes[order(phenotypes$phenotype),]
    }
    phenotypes$seq = 1:nrow(phenotypes)
    
    #Limit to phenotype and seq, as they are the only relevant columns
    #Include value as min.value for annotation purposes
    phenotypes=phenotypes[,c("phenotype","seq","value")]
    names(phenotypes)[3]="min.value"
    
    #Add sequence information
    d=merge(d,phenotypes,by="phenotype")
    
    #If labeling the phenotype groups
    if (x.group.labels) ticks=sort(aggregate(seq ~ groupnum, d, mean)$seq)
    #Define the max y axis value if not provided
    if(missing(max.y)) max.y=ceiling(max(d$value))
    

    #Generate the inital plot
    plot=ggplot(d,ylab=y.axis.label,xlab=x.axis.label)

    #Include lines for significance thresholds
    if (!missing(suggestive.line)&!is.na(suggestive.line)) plot=plot+geom_hline(yintercept=suggestive.line,colour="blue", alpha=I(1/3),size=1)
    if (!missing(significant.line)&!is.na(significant.line)) plot=plot+geom_hline(yintercept=significant.line,colour="red",alpha=I(1/3),size=1)

    plot=plot+aes(seq,value,size=size,colour=color)
    if(!sizes) plot=plot+scale_size(range=c(point.size,point.size),guide="none")
    #Add points
    plot=plot+geom_point() 
    
    if(missing(color.palette)) {
      color.palette = unique(d[order(d$seq),]$color)
      names(color.palette)=color.palette
    } else {
      names(color.palette) = unique(d[order(d$seq),]$color)
    }
    
    #Color as defined 
    plot = plot + scale_colour_manual(values= color.palette, guide="none") 

    #TODO: X phenotype labels
    #If label the X axis with the groups if requested
    if (x.group.labels) {
      plot=plot+scale_x_continuous(name=x.axis.label, limits=c(1,max.x), breaks=ticks, labels=(unique(d$group)), expand=c(.01,0))
      
    } else {
      plot=plot+scale_x_continuous(name=x.axis.label, limits=c(1,max.x), breaks=c(-100), labels=c(""), expand=c(.015,0))
    }
    
    #Set the Y scale and labels
    plot=plot+scale_y_continuous(y.axis.label, limits=c(0,max.y), breaks=seq(0,max.y,y.axis.interval), expand=c(0,.2))
    
    #Set the default theme
    plot=plot+theme(
      panel.background=element_blank(), 
      panel.grid.minor=element_blank(),
      axis.text.x=element_text(size=size.x.labels, colour="black", angle=-40, hjust=0, vjust=1), 
      axis.text.y=element_text(size=size.y.labels, colour="black"), 
      axis.line =element_line(colour="black"),
      axis.ticks=element_line(colour="black")
    ) 
 
    #Hide the legend by default
    plot = plot+theme(legend.position = "none")
    
    #Add OR information
    if(sizes){			
      plot= suppressWarnings(plot + theme(legend.position="right") +
        scale_size("Size", range = c(point.size, 2*point.size)))
    }
    if(direction){
      plot=plot+aes(shape = factor(direction), fill=color) + 
        scale_shape("Direction", solid = TRUE) + scale_shape_manual(values=c(25,24)) + 
        scale_fill_manual(values= color.palette, guide="none") 
    }		
    #If annotation present, start the definitions, otherwise skip it
    if (annotate) 	{
      d$annotate=F
      #If provided with a list of phenotypes to annotate, select those.
      if(!missing(annotate.list)) d[d$phenotype %in% annotate.list, ]$annotate=T
      #Include those above the given threshold
      if((!missing(annotate.level))&(sum(d$value>=annotate.level)>0)) d[d$value>=annotate.level, ]$annotate=T
      #Select only the largest value for each phenotype to label, unless requested otherwise
      if(sum(d$value!=d$min.value)>0&annotate.only.largest) d[d$value!=d$min.value,]$annotate=F
      #If no base descriptions were included, add empty ones
      if(!annotate.phenotype.description) {
        d$description=""
      }     
      #If requested the phenotype annotation, add it
      if(annotate.phenotype) d$description = d$phenotype %+% ifelse(d$description=="","",":" %+% d$description )
      
      #Add snp annotation to the phenotype label
      if (annotate.snp.w.phenotype) d$description = d$snp %+% ifelse(d$description=="","",":" %+% d$description )
      
      #Cap annotation length
      d$description = substr(d$description,1,60)
      #Add leading space
      d$description = paste0("  ",d$description)
      #If lower case labels are requested, lower case them.
      if(lc.labels) d$description=tolower(d$description)
      if(sum(d$annotate)==0) {
        warning("Annotation requested, but no points met criteria")
      } else {
        #Add annotations
        if(annotate.phenotype.description|annotate.phenotype|annotate.snp.w.phenotype) plot = plot + geom_text(aes(label=description),colour="black",data=d[d$annotate,],hjust=0,size=annotate.size,angle=annotate.angle)
        #Add SNP annotations
        if(annotate.snp) plot = plot + geom_text(aes(label=snp),colour="black",data=d[d$annotate,],hjust=0,size=annotate.size,angle=annotate.snp.angle)
      }
    } 
    
    #Add the title
    plot=plot+labs(title=title) + theme(title=element_text(size=8))
    
    plot
  }
