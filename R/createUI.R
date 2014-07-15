createUI <- function(results) {
pageWithSidebar(
  
  # Application title
  headerPanel("PheWAS results"),
  
  # Sidebar with controls to select the variable to plot 
  sidebarPanel(
    selectInput("snp", "Plot for which SNP:", 
                choices = c("All",as.character(unique(results$snp)))
                ),
    selectInput("mtcor", "Significance Threshold:", 
                choices = c("Alpha","Bonferroni","SimpleM-product"), selected="Bonferroni"),
    textInput("alpha", "Unadjusted Alpha:", attributes(results)$alpha),
    
    checkboxInput("annotate", "Show PheWAS code Annotation", T),
    conditionalPanel(
      condition = "input.annotate == true",
      textInput("annotate_level", "Level for annotation:", "4e-3"),
      textInput("annotate_angle", "Angle for annotation:", "8")
    )
    
  ),
  
  # Show the caption and plot 
  mainPanel(
    h3(textOutput("caption")),
    
    plotOutput("plot")
  )
)
}