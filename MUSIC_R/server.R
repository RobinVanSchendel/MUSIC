function(input, output, session) {
  
  # Render dynamic tabs based on input
  output$dynamicTabs <- renderUI({
    generateTabs(input$data_input)
  })
  
  ##store some data
  volcanodata <- reactiveVal()
  volcanofocuseddata <- reactiveVal()
  ##TODO: implement this to increase speed!
  xydata <- reactiveVal()
  umap_gene_data <- reactiveVal()
  
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
                           hover = hoverOpts("plot_waterfall", delay = 10, delayType = "debounce"),
                           brush = "plot_waterfallplot_brush"))
  })
  
  # Render UI for XY plot
  output$UIXYplot <- renderUI({
    withSpinner(plotOutput("XYplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
               hover = hoverOpts("plot_xy", delay = 10, delayType = "debounce"),
               brush = "plot_xy_brush",
               dblclick = "plot_xy_dblclick"))
  })
  
  output$UIVolcanoplot <- renderUI({
    withSpinner(plotOutput("Volcanoplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_volcano", delay = 10, delayType = "debounce"),
                           brush = "plot_volcano_brush",
                           dblclick = "plot_volcano_dblclick"))
  })
  
  
  output$UIVolcanoFocusedplot <- renderUI({
    withSpinner(plotOutput("VolcanoFocusedplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_volcano_focused", delay = 10, delayType = "debounce"),
                           brush = "plot_volcano_focused_brush"))
  })
  
  ##umap-gene plot
  output$UIUmapPlot <- renderUI({
    withSpinner(plotOutput("UmapPlot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_umap_focused", delay = 10, delayType = "debounce"),
                           brush = "plot_umap_focused_brush"))
  })
  
  
  ##volcano brush
  observeEvent(input$plot_volcano_brush, {
    data = volcanodata()
    brushed <- brushedPoints(data, input$plot_volcano_brush)
    output$volcano_table <- renderDT(makeInteractiveTable(brushed))
  })
  
  ##xy brush
  observeEvent(input$plot_xy_brush, {
    req(xydata())
    data = xydata()
    brushed <- brushedPoints(data, input$plot_xy_brush)
    output$xy_table <- renderDT(makeInteractiveTable(brushed))
  })
  
  ##umap_gene brush
  observeEvent(input$plot_umap_focused_brush, {
    req(umap_gene_data())
    data = umap_gene_data()
    brushed <- brushedPoints(data, input$plot_umap_focused_brush)
    ##no need for some of these columns
    brushed <- brushed %>% select(-size, -alpha)
    output$umap_gene_table <- renderDT(makeInteractiveTable(brushed))
  })
  
  
  
    ##for zooming purposes
  ranges_volcano_brush <- reactiveValues(x = NULL, y = NULL)
  ranges_xy_brush <- reactiveValues(x = NULL, y = NULL)
  
  
  # Zoom on double-click
  observeEvent(input$plot_volcano_dblclick, {
    updateZoomRanges(input$plot_volcano_brush, ranges_volcano_brush)
  })
  # Zoom on double-click
  observeEvent(input$plot_xy_dblclick, {
    updateZoomRanges(input$plot_xy_brush, ranges_xy_brush)
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
  
  ##TODO change to only read in data when umap is needed
  observe({
    message("Reading in umap_gene_data")
    umap_data = read_in_umap_gene_data()
    umap_gene_data(umap_data)
  })
  
  output$Waterfallplot <- renderPlot({
    req(input$waterfall_type, db_con_genome_wide())
    con <- db_con_genome_wide()
    df <- prepareWaterfallData(con, input)
    
    df <- df %>% group_by(Alias, Type) %>%
      arrange(desc(fraction)) %>%
      mutate(rowNr = row_number()) %>%
      ungroup()
    
    sds = df %>% group_by(Alias, Type) %>%
      summarise(mean = mean(fraction), sd = sd(fraction))
    
    ##ensure order is correct
    #levels = c("")
    #df$Type = as.factor(df$Type, levels = c())
    ##ensure the waterfall types are coloured
    WATERFALL_TYPES = setNames(c(
      "#FF8C00", "#124E8B","#05878C","#EC71A7", "#B03060", "#E6332A"),
      GENOME_WIDE_TYPES
    )
    #only keep selected types
    WATERFALL_TYPES <- WATERFALL_TYPES[names(WATERFALL_TYPES) %in% input$waterfall_type]
    
    df$Type = factor(df$Type, levels = names(WATERFALL_TYPES))
    sds$Type = factor(sds$Type, levels = names(WATERFALL_TYPES))
    
    
    ggplot(df, aes(x = rowNr, y = fraction)) +
      geom_hline(data = sds, aes(yintercept = mean+3*sd), linetype = "dashed", color = "grey30") +
      geom_hline(data = sds, aes(yintercept = mean-3*sd), linetype = "dashed", color = "grey30") +
      geom_point(size = .5, alpha = .6) +
      facet_wrap2(~Type, scales = "free",
                  strip = strip_themed(
                    background_x = elem_list_rect(
                      fill = WATERFALL_TYPES,
                      color = "black"
                    ),
                    text_x = elem_list_text(
                      color = "white",
                      face  = "bold"
                    ))
                 ) +
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
    ##store xydata in a reactive
    xydata(df)
    base_plot <- buildXYPlot(df, input, plotX, plotY)
    
    highlightPart <- getHighlightedGenes(df, input)
    topGenes <- getTopGenes(df, input, plotX)
    
    final_plot <- addHighlightsAndMarginals(
      base_plot, df, input, plotX, plotY,
      highlightPart, topGenes$up, topGenes$down
    )
    
    final_plot +
      facet_wrap(Alias ~ ., ncol = 2, labeller = as_labeller(GENOME_WIDE_LABEL_MAP)) +
      theme_object() +
      coord_cartesian(xlim = ranges_xy_brush$x, ylim = ranges_xy_brush$y)
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
    
    final_plot <- final_plot + coord_cartesian(xlim = ranges_volcano_brush$x, ylim = ranges_volcano_brush$y)
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
  
  output$UmapPlot <- renderPlot({
    req(umap_gene_data())
    df = umap_gene_data()
    
    base_plot = ggplot(df, aes(x = X2, y = X1, col = Pathway1, alpha = alpha, size = size)) +
      geom_point()
    
    final_plot = addLabelsUMAP(base_plot, df, input)
    
    #  geom_text_repel(data=subset(df, Gene != "NonTargeting"), aes(label = Gene), max.overlaps = 100)+
    final_plot = final_plot +
      theme_object() +
      scale_color_manual(values = UMAP_GENE_COLORS)+
      ##remove legends for size and alpha
      guides(size = "none", alpha = "none")+
      theme(legend.position = "top") +
      theme_object() +
      NULL
    final_plot
  })
  
  ##for the hover
  output$hover_xy <- renderUI({
    req(xydata())
    hover <- input$plot_xy
    if(is.null(hover)){
      return()
    }
    ##get the data from the reactive
    df <- xydata()
    createHoverTooltip(hover, df)
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
  
  ##hover for umap
  output$hover_umap_gene <- renderUI({
    hover <- input$plot_umap_focused
    if(is.null(hover)){
      return()
    }
    df <- umap_gene_data()
    createHoverTooltip(hover, df)
  })
    
    
  
  
  
  # Optional: observe tab selection
  observeEvent(input$tabs, {
    if (input$tabs == "XY") {
      # Add logic if needed when XY tab is selected
    }
  })
}
