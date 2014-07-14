source("http://knowledgemap.mc.vanderbilt.edu/research/sites/default/files/R_PheWAS_Install.R")
#Load the PheWAS package
suppressWarnings(library(PheWAS))
#Set the random seed so it is replicable
set.seed(1)
#Generate some example data
ex=generateExample()
#Extract the two parts from the returned list
id.icd9.count=ex$id.icd9.count
genotypes=ex$genotypes
#Create the PheWAS code table- translates the icd9s, adds exclusions, and reshapes to a wide format
phenotypes=createPhewasTable(id.icd9.count,fast=TRUE)
#Run the PheWAS
results=phewas(phenotypes,genotypes,cores=4,significance.threshold=c("bonferroni"))
#Plot the results
phewasManhattan(results, annotate.angle=0)
#Add PheWAS descriptions
results_d=addPhewasDescription(results)
#List the significant results
results_d[results_d$bonferroni&!is.na(results_d$p),]
#List the top 10 results
results_d[order(results_d$p)[1:10],]
