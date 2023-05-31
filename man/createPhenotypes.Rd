\name{createPhenotypes}
\alias{createPhenotypes}
\title{
Creates a phenotype table from id, ICD9CM, ICD10CM (or phecode, etc), data.
}
\description{
This function takes a data frame with four columns: id, vocabulary_id, code, and index. It returns a wide table with phecodes as TRUE/FALSE/NA. It can optionally use the PheWAS exclusion criteria.
}
\usage{
createPhenotypes(id.vocab.code.index, min.code.count=2, 
           add.phecode.exclusions=T, translate=T, id.sex, 
           full.population.ids=unique(id.vocab.code.index[[1]]),
           aggregate.fun=PheWAS:::default_code_agg, 
           vocabulary.map=PheWAS::phecode_map,
           rollup.map=PheWAS::phecode_rollup_map,
           exclusion.map=PheWAS::phecode_exclude,
           sex.restriction=PheWAS::sex_restriction)
}

\arguments{
  \item{id.vocab.code.index}{
Data frame with four columns of information: id, vocabulary_id, code, and index. The id and index columns can have other names, but id must be consistent among input files. The vocabulary_id and code must match up with the vocabulary.map file. The default supports the vocabularies "ICD9CM" and "ICD10CM". Code contains the raw code value.

Note the beta ICD10 map is provided in \code{\link[PheWAS:phecode_map_icd10]{PheWAS::phecode_map_icd10}}, which can be provided as a parameter to \code{vocabulary.map}.
}
  \item{min.code.count}{
The minimum code count to be considered a case. NA results in a continuous output.
}
  \item{add.phecode.exclusions}{
Apply PheWAS exclusions to phecodes.
}
  \item{translate}{
Should the input be translated to phecodes? Defaults to TRUE. Generally recommended, though can be skipped if phecodes are provided.
}
  \item{aggregate.fun}{
Aggregate function for duplicated phenotypes (phecodes, etc) in an individual. The default supports will use \code{sum} for numeric values, otherwise it will count the distinct values, eg, for dates.
}
  \item{id.sex}{
If supplied, restrict the phecodes by sex. This should be a data frame with the first column being the id and the second the sex, "M" or "F", of the individual. Individuals with any other specification will have all sex specific phenotypes set to NA.
}
  \item{full.population.ids}{
List of IDs in the "complete" population. This allows for individuals with no observed codes to have appropriate "control" status, eg 0s or FALSE in every field.
}
  \item{aggregate.fun}{
Aggregate function for duplicated phenotypes (phecodes, etc) in an individual. The default supports a naive "distinct date" approach. Use \code{sum} to support count data.
}
  \item{vocabulary.map}{
Map between supplied vocabularies and phecodes. Allows for custom phecode maps. By default uses \code{\link[PheWAS:phecode_map]{PheWAS::phecode_map}}, which supports ICD9CM (v1.2) and ICD10CM (beta-2018). The package also includes the ICD10 beta map (\code{\link[PheWAS:phecode_map_icd10]{PheWAS::phecode_map_icd10}}), which can be used in this parameter.
}
  \item{rollup.map}{
Map between phecodes and all codes that they expand to, eg parent codes. By default uses the PheWAS::phecode_rollup_map.
}
  \item{exclusion.map}{
Map between phecodes and their exclusions. By default uses the PheWAS::phecode_exclude.
}
  \item{sex.restriction}{
List that specifies sex-specific phecodes. By default uses the PheWAS::sex_restriction.
}
}
\details{
By default, this function returns a wide format data frame with boolean phenotypes suitable for PheWAS analysis. Specifying a \code{min.code.count=NA} will permit continuous code count phenotypes.

The default exclusions can be skipped with \code{add.exclusions=F}. In conjuntion with \code{translate=F} (and optionally adjusting \code{min.code.count} and \code{aggregate.fun}), one can use this function as a simple reshaping wrapper.
}
\value{
A data frame. The first column contains the supplied id for each individual (preserving the name of the original column). The following columns are all present phewas codes. They contain T/F/NA for case/control/exclude or continuous/NA if min.code.count was NA.
}
\author{
Robert Carroll
Laura Wiley
}
\examples{
#Simple example
id_icd9_count=data.frame(id=c(1,2,2,2,3),vocabulary_id="icd9",code=c("714","250.11","714.1","714","250.11"),
  count=c(1,5,1,1,0))
createPhenotypes(id_icd9_count)
\donttest{
#Complex example
ex=generateExample(n=500,hit="335")
#Extract the two relevant parts from the returned list
id.vocab.code.count=ex$id.vocab.code.count
id.sex=ex$id.sex
#Create the phecode table- translates the codes, adds 
#exclusions, and reshapes to a wide format.
#Sum up the counts in the data where applicable.
phenotypes=createPhenotypes(id.vocab.code.count, 
  aggregate.fun=sum, id.sex=id.sex)
#Create the phecode table for a PheWAS
phenotypes=createPhewasTable(id.icd9.count)
phenotypes[1:10,1:10]
}
}
\keyword{ utilities }
