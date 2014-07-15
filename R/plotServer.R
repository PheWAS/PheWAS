plotServer <-function(results) { function(input, output) {
  
  # Get the reactive value
  annotate <- reactive({
    input$annotate
  })
  annotate_level <- reactive({
    as.numeric(input$annotate_level)
  })
  annotate_angle <- reactive({
    as.numeric(input$annotate_angle)
  })
  my_snp <- reactive({
    input$snp
  })
  alpha <- reactive({
    as.numeric(input$alpha)
  })
  a.corrected <- reactive({
    if(input$mtcor=="Alpha") {
      return(as.numeric(input$alpha))
    }
    if(input$mtcor=="Bonferroni") {
      return(as.numeric(input$alpha)/nrow(results))
    }
    if(input$mtcor=="SimpleM") {
      if(is.null(attributes(results)$simplem.product.meff)) warning("Warning: Must calculate SimpleM product in the phewas call to use this threshold")
      return(as.numeric(input$alpha))
    }
  })
  
  
  # Return the formula text for printing as a caption
  output$caption <- renderText({
    my_snp()
  })
  
  # Generate a plot of the requested snp 
  output$plot <- renderPlot({
    if(my_snp()=="All") { my_sig=results
    } else { my_sig = results[results$snp==my_snp(),] }
    plot=phewasManhattan(my_sig, 
                     annotate.level=annotate_level(), 
                     annotate.phenotype=annotate(), 
                     annotate.angle=annotate_angle(),
                     suggestive.line=alpha(),
                     genomewide.line=a_corrected())
    print(plot)
  })
}}