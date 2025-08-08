function(input, output, session) {
  
  # Render dynamic tabs based on input
  output$dynamicTabs <- renderUI({
    generateTabs(input$data_input)
  })
  
  ##store some data
  volcanodata <- reactiveVal()
  
  # Reactive database connection
  db_con <- reactive({
    getDbConnection(input$data_input)
  })
  
  # Render gene barcode table
  output$gene_barcode <- DT::renderDataTable({
    req(input$data_input)
    renderGeneBarcodeTable(input$data_input)
  })
  
  # Render UI for XY plot
  output$UIXYplot <- renderUI({
    withSpinner(plotOutput("XYplot", height = input$plotHeight, width = input$plotWidth,
               hover = hoverOpts("plot_xy", delay = 10, delayType = "debounce")))
  })
  
  output$UIVolcanoplot <- renderUI({
    withSpinner(plotOutput("Volcanoplot", height = input$plotHeight, width = input$plotWidth,
                           hover = hoverOpts("plot_volcano", delay = 10, delayType = "debounce"),
                           brush = "plot_volcano_brush"))
  })
  
  
  
  # Reactive theme object
  theme_object <- reactive({
    generateTheme(input)
  })
  
  # Update picker input for gene highlighting
  observe({
    req(db_con())
    con <- db_con()
    genes <- tbl(con, "barcodeCount") %>% select(Gene) %>% distinct() %>% collect() %>% pull()
    updatePickerInput(session, inputId = "highlightGene", choices = genes)
    
    ##make this data static
    volcano_data <- prepareVolcanoData(con, input$VolcanoType)
    volcanodata(volcano_data)  
  })
  
  # Render XY plot
  output$XYplot <- renderPlot({
    req(db_con())
    con <- db_con()
    plotX <- input$XYX
    plotY <- input$XYY
    
    df <- prepareXYData(con, plotX, plotY)
    base_plot <- buildXYPlot(df, input, plotX, plotY)
    
    highlightPart <- getHighlightedGenes(df, input)
    topGenes <- getTopGenes(df, input, plotX)
    
    final_plot <- addHighlightsAndMarginals(
      base_plot, df, input, plotX, plotY,
      highlightPart, topGenes$up, topGenes$down
    )
    
    final_plot +
      facet_wrap(Alias ~ ., ncol = 2, labeller = as_labeller(GENOME_WIDE_LABEL_MAP)) +
      theme_object()
  })
  
  # Render Volcano plotXYplot
  ##TODO: check WT
  output$Volcanoplot <- renderPlot({
    req(volcanodata())
    
    plotType = input$VolcanoType
    
    df <- volcanodata()
    highlightPart <- getHighlightedGenes(df, input)
    topGenes <- getTopGenesVolcano(df, input, plotType)
    
    base_plot <- ggplot(data = df, aes(y = Pvalue, x = log2fraction)) +
      geom_point() +
      NULL
    
    final_plot <- addHighlightsAndMarginalsVolcano(
      base_plot, input, highlightPart, topGenes$up, topGenes$down)
    
    final_plot = final_plot + facet_wrap(Alias ~ ., ncol = 2, labeller = as_labeller(GENOME_WIDE_LABEL_MAP)) +
      theme_object() +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
      NULL
    
    ##TODO: check brush
    #if(!is.null(input$plot_volcano_brush)){
    #  brush = input$plot_volcano_brush
    #  print(brush)
    #}
    final_plot
  })
  
  ##for the hover
  output$hover_volcano <- renderUI({
    hover <- input$plot_volcano
    if(is.null(hover)){
      return()
    }
    df <- volcanodata()
    createHoverTooltip(hover, df)
  })
  
  # Optional: observe tab selection
  observeEvent(input$tabs, {
    if (input$tabs == "XY") {
      # Add logic if needed when XY tab is selected
    }
  })
}
