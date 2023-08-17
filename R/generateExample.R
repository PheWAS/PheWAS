#'  Generate example data for PheWAS
#'
#'  \code{generateExample} creates simulated ICD9CM/ICD10CM, id.sex, and genotype data to test PheWAS methods. Requires the package \code{MASS}.
#'
#' The output from this function will be random results with a significant association found between the genotype data and the selected phecode. It does not represent actual treatment data, though some highly connected phenotypes can also have strong associations due to the data creation methods.
#'
#' @param n The number of individuals in the study.
#' @param phenotypes.per The base number of phenotypes per individual.
#' @param hit The phecode for which an artificial hit will be generated.
#'
#' @return List with three named items:
#' \item{id.vocab.code.count}{A data frame of columns id, vocabulary_id, ICD code, and count suitable for transformation by \code{createPhenotypes}.}
#' \item{genotypes}{A data frame of columns id and rsEXAMPLE, suitable for use in \code{phewas} after joining with phenotypic data.}
#' \item{id.sex}{A data frame of columns id and sex, suitable for use in \code{phewas} and \code{createPhenotypes}.}
#' @export
#' @import dplyr
#' @examples generateExample()
generateExample <- function(n=5000,phenotypes.per=10,hit="335") {
  phecode=unique(PheWAS::pheinfo$phecode)
  #Exclude the code to add
  phecode=phecode[phecode!=hit]
  #Assign individuals random phenotypes
  random=data.frame(id=rep.int(1:n,phenotypes.per),phecode="",count=0)
  random$phecode=sample(phecode,nrow(random),replace=TRUE)
  #Create the signal
  signal=as.data.frame(mvrnorm(n=n,mu=c(0,0),Sigma=rbind(c(.5,.1),c(.1,1))))
  names(signal)=c("phenotype","rsEXAMPLE")
  #Normalize the phenotype to logical and the genotype to 0,1,2.
  signal$phenotype=signal$phenotype>.25
  signal$rsEXAMPLE=signal$rsEXAMPLE+abs(min(signal$rsEXAMPLE))
  signal$rsEXAMPLE=floor(signal$rsEXAMPLE*2.99/max(signal$rsEXAMPLE))
  signal$id=1:n
  signal$count=0
  signal$phecode=hit
  random=rbind(random,signal[signal$phenotype,c("id","phecode","count")])
  #Map to ICD9/10CM codes and pick one per phecode
  random=inner_join(random,PheWAS::phecode_rollup_map, by=c("phecode"="code")) %>%
    inner_join(PheWAS::phecode_map,by=c("phecode_unrolled"="phecode")) %>%
    group_by(id, phecode_unrolled) %>% sample_n(1) %>% ungroup()
  random$count= rpois(nrow(random),4)
  random[random$phecode==hit,]$count=random[random$phecode==hit,]$count+2

  random=random %>% filter(count>0) %>% select(id,vocabulary_id,code, count) %>% unique()
  return(list(id.vocab.code.count=random,
              genotypes=signal[,c("id","rsEXAMPLE")],
              id.sex=data.frame(id=signal$id,sex=sample(c("F","M"), nrow(signal), replace=TRUE))
  ))
}

