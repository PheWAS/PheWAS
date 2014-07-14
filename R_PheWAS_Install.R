#Install the recommended packages
if(!require("plyr")) {install.packages("plyr")}
if(!require("snowfall")) {install.packages("snowfall")}
if(!require("shiny")) {install.packages("shiny")}
if(!require("MASS")) {install.packages("MASS")}
if(!require("ggplot2")) {install.packages("ggplot2")}
if(!require("SparseM")) {install.packages("SparseM")}
if(!require("meta")) {install.packages("meta")}
#If necessary, install the R PheWAS package
cur_ver="0.9.3.1"
if(!(require("PheWAS")&&(packageVersion("PheWAS")>=package_version(cur_ver)))) {
  #Retrieve the package file
  download.file(paste0("http://knowledgemap.mc.vanderbilt.edu/research/sites/default/files/PheWAS_",cur_ver,".tar.gz"),paste0("~/PheWAS_",cur_ver,".tar.gz"))
  #Install the PheWAS package
  install.packages(paste0("~/PheWAS_",cur_ver,".tar.gz"),repos=NULL, type="source")
}