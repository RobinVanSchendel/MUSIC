function(input, output, session) {
  
  # Render dynamic tabs based on input
  output$dynamicTabs <- renderUI({
    generateTabs(input$data_input)
  })
  
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
    plotOutput("XYplot", height = input$plotHeight, width = input$plotWidth,
               hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce"))
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
  
  # Optional: observe tab selection
  observeEvent(input$tabs, {
    if (input$tabs == "XY") {
      # Add logic if needed when XY tab is selected
    }
  })
}
