#' Startup Message
#'
#' @param libname N/A
#' @param pkgname N/A
#'
#' @return A start up message
#' @export
#'
#' @examples n/a
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the new version of PheWAS. This version has many updates; please see https://github.com/PheWAS/PheWAS/tree/legacy for the legacy release if needed. Check ?PheWAS for more documentation")
}
