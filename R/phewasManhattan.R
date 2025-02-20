#' Plotting methods for PheWAS results, phenotype association studies, or 
#' generic phenotype values
#'
#' This functions permit the plotting of a variety of data types found when
#' working with phenotypes. The function is comprised of sequentially nested 
#' calls as outlined below. Supplying lower level parameters will permit tweaks 
#' to the higher level functions as well.
#'
#' These functions are nested: phewasManhattan calls phenotypeManhattan which 
#' calls phenotypePlot. 
#' 
#' phewasManhattan generates a PheWAS Manhattan plot and can be used with a 
#' result data frame from phewas using phecodes. This function only adds PheWAS
#'  code groupings and descriptions as applicable before calling phenotypeManhattan.
#' phenotypeManhattan performs the transformations to the data and paramters for
#'  a -log10 scale. It is useful in cases of phewas results using other 
#'  phenotypes besides phecodes
#'  
#' phenotypePlot does the actual plotting and can be used to visualize non-p-value data.
#' All of these functions return ggplot2 objects. Aspects of the plots can be 
#' modified post hoc, for instance to adjust the point size scale one can use 
#' plot + scale_size(<parameters>)
#'
#' @param d 	Data frame containing phenotype and p, or in the case of 
#' phenotypePlot, phenotype and value
#' @param annotate.phenotype.description Contains either TRUE or a data frame 
#' that contains description annotations for each phenotype. Should contain 
#' columns "phenotype" and "description". Missing or FALSE yields no description
#'  annotation.
#' @param pheinfo.map a map of phecode descriptions for visualization
#' @param ... n/a
#'
#' @return A ggplot2 plot. It contains PheWAS codes against their -log10
#'  transformed p-values for phewasManhattan, generic phenotypes versus -log10 
#'  transformed p-values for phenotypeManhattan, and generic phenotypes versus 
#'  generic values for phenotypePlot.
#' @export
#'
#' @examples 
#' 
#' phewasManhattan(PheWAS:::test_phewas, pheinfo.map = PheWASmaps::pheinfo_X)
phewasManhattan <-
  function(d, annotate.phenotype.description=T
           , pheinfo.map = PheWASmaps::pheinfo
           ,  ...) {
    if(sum(c("phenotype","p") %in% names(d))<2 ) 
      stop("Data input must contain columns phenotype and p.")
    if(!is(d$phenotype, 'character')) {

      
        stop("Non-character phenotypes passed in, so an accurate
             phecode mapping is not possible.")
      
    }
    #Check to see if it looks 0-padded
    if(min(nchar(d$phenotype))<3) stop('Phenotypes with length <3 observed,ensure they are are 0-padded')

    #Add the groups and phecode descriptions as requested
    d=addPhecodeInfo(d,groupnums =T, groupcolors = T
                     ,pheinfo = pheinfo.map
                     )

    phenotypeManhattan(d, annotate.phenotype.description=annotate.phenotype.description, ...)
  }