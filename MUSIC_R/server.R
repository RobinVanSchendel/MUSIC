function(input, output, session) {
  
  # Render dynamic tabs based on input
  output$dynamicTabs <- renderUI({
    generateTabs(input$data_input)
  })
  
  ##store some data
  volcanodata <- reactiveVal()
  volcanofocuseddata <- reactiveVal()
  
  # Reactive database connection
  db_con_genome_wide <- reactive({
    return(genome_wide_conn())
  })
  db_con_candidate <- reactive({
    return(candidate_conn())
  })
  
  # Render gene barcode table
  output$gene_barcode <- DT::renderDataTable({
    req(input$data_input)
    renderGeneBarcodeTable(input$data_input)
  })
  
  
  # Render UI for Waterfall plot
  output$UIWaterfall <- renderUI({
    withSpinner(plotOutput("Waterfallplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_waterfall", delay = 10, delayType = "debounce")))
  })
  
  # Render UI for XY plot
  output$UIXYplot <- renderUI({
    withSpinner(plotOutput("XYplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
               hover = hoverOpts("plot_xy", delay = 10, delayType = "debounce")))
  })
  
  output$UIVolcanoplot <- renderUI({
    withSpinner(plotOutput("Volcanoplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_volcano", delay = 10, delayType = "debounce"),
                           brush = "plot_volcano_brush"))
  })
  
  
  output$UIVolcanoFocusedplot <- renderUI({
    withSpinner(plotOutput("VolcanoFocusedplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_volcano_focused", delay = 10, delayType = "debounce"),
                           brush = "plot_volcano_focused_brush"))
  })
  
  
  
  
  
  # Reactive theme object
  theme_object <- reactive({
    generateTheme(input)
  })
  
  # Update picker input for gene highlighting
  observeEvent(input$data_input, {
    if(input$data_input == 1){
      req(db_con_genome_wide())
      print("update highlightGene")
      con <- db_con_genome_wide()
      genes <- tbl(con, "barcodeCount") %>% select(Gene) %>% distinct() %>% collect() %>% pull()
      updatePickerInput(session, inputId = "highlightGene", choices = genes)
      print("update highlightGene finished")
    }
    else if(input$data_input == 2){
      req(db_con_candidate())
      print("update highlightGene")
      con <- db_con_candidate()
      genes <- tbl(con, "barcodeCount") %>% select(Gene) %>% distinct() %>% collect() %>% pull()
      updatePickerInput(session, inputId = "highlightGeneFocused", choices = genes)
      print("update highlightGeneFocused finished")
    }
  })
  
  observeEvent(input$VolcanoType, {
      req(db_con_genome_wide())
      con <- db_con_genome_wide()
      ##make this data static
      print("prepareVolcanoData")
      volcano_data <- prepareVolcanoData(con, input$VolcanoType)
      print("got prepareVolcanoData")
      volcanodata(volcano_data)
  })
  
  observeEvent(input$VolcanoFocusedType, {
    req(db_con_candidate())
    con <- db_con_candidate()
    print("prepareVolcanoFocusedData")
    ##make this data staic
    volcano_focused_data <- prepareVolcanoFocusedData(con, input$VolcanoFocusedType)
    volcanofocuseddata(volcano_focused_data)
  })
  
  output$Waterfallplot <- renderPlot({
    req(db_con_genome_wide())
    con <- db_con_genome_wide()
    df <- prepareWaterfallData(con)
    
    df <- df %>% group_by(Alias, Type) %>%
      arrange(fraction) %>%
      mutate(rowNr = row_number())
    
    sds = df %>% group_by(Alias, Type) %>%
      summarise(mean = mean(fraction), sd = sd(fraction))
    
    ggplot(df, aes(x = rowNr, y = fraction)) +
      geom_hline(data = sds, aes(yintercept = mean+3*sd), linetype = "dashed", color = "grey30") +
      geom_hline(data = sds, aes(yintercept = mean-3*sd), linetype = "dashed", color = "grey30") +
      geom_point(size = .5, alpha = .6) +
      facet_wrap(~Type, scales = "free") +
      theme_object() +
      NULL
  })
  
  
  # Render XY plot
  output$XYplot <- renderPlot({
    req(db_con_genome_wide())
    con <- db_con_genome_wide()
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
    if(nrow(df) == 0){
      return()
    }
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
  
  # Render Volcano plotXYplot
  ##TODO: check WT
  output$VolcanoFocusedplot <- renderPlot({
    req(volcanofocuseddata())
    
    plotType = input$VolcanoType
    
    df <- volcanofocuseddata()
    highlightPart <- getHighlightedGenes(df, input)
    topGenes <- getTopGenesVolcano(df, input, plotType)
    
    base_plot <- ggplot(data = df, aes(y = mutEvents, x = log2fraction)) +
      geom_point() +
      NULL
    
    final_plot <- addHighlightsAndMarginalsVolcano(
      base_plot, input, highlightPart, topGenes$up, topGenes$down)
    
    final_plot = final_plot + facet_wrap(Alias ~ ., ncol = 3, labeller = as_labeller(FOCUSED_MAP)) +
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
  
  output$hover_volcano_focused <- renderUI({
    hover <- input$plot_volcano_focused
    if(is.null(hover)){
      return()
    }
    df <- volcanofocuseddata()
    createHoverTooltip(hover, df)
  })
  
  
  
  
  # Optional: observe tab selection
  observeEvent(input$tabs, {
    if (input$tabs == "XY") {
      # Add logic if needed when XY tab is selected
    }
  })
}
