shiny=function(results){
if(!require(shiny)) {stop("Requires the package shiny to function.")}
runApp(list(ui=createUI(results),server=plotServer(results)))
}