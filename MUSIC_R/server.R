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
  heatmap_data <- reactiveVal()
  
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
  
  ##heat map plot
  output$UIHeatmapPlot <- renderUI({
    withSpinner(plotOutput("HeatmapPlot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_heatmap_hover", delay = 10, delayType = "debounce"),
                           brush = brushOpts(id = "plot_heatmap_brush", resetOnNew = T),
                           dblclick = "plot_heatmap_dblclick"))
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
  ranges_heatmap_brush <- reactiveValues(x = NULL, y = NULL)
  
  
  # Zoom on double-click
  observeEvent(input$plot_volcano_dblclick, {
    updateZoomRanges(input$plot_volcano_brush, ranges_volcano_brush)
  })
  # Zoom on double-click
  observeEvent(input$plot_xy_dblclick, {
    updateZoomRanges(input$plot_xy_brush, ranges_xy_brush)
  })
  # Zoom on double-click
  observeEvent(input$plot_heatmap_dblclick, {
    updateZoomRanges(input$plot_heatmap_brush, ranges_heatmap_brush)
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
      message("prepareVolcanoData Type",input$VolcanoType)
      volcano_data <- prepareVolcanoData(con, input$VolcanoType)
      message("got prepareVolcanoData")
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
    if(is.null(umap_gene_data())){
      message("Reading in umap_gene_data")
      umap_data = read_in_umap_gene_data()
      umap_gene_data(umap_data)
    }
    if(is.null(heatmap_data())){
      message("Reading in heatmap_data")
      heatmap_df = read_in_heatmap_data()
      heatmap_data(heatmap_df)
    }
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
  
  ##heatmap plot
  output$HeatmapPlot <- renderPlot({
    req(heatmap_data())
    heatmap = heatmap_data()
    
    ##convert to a matrix first, takes care of NA and limits (-2, 2)
    matrix = convert_data_frame_to_matrix(heatmap)
    
    highlight_barcodes = c("Polq-2", "Polq-3", "Polq-5")
    
    ##retrieve a df in long format, but uses heatmap.2 for its ordering
    ##currently uses integers to plot as that allows for additional components to be added to the plot
    ##e.g. in (-10 - -1)
    long_df = retrieve_long_format_sorted_by_heatmap(matrix)
    
    ggplot(long_df, aes(x=x_num, y = y_num, fill = log2fraction)) +
      geom_tile() +
      # Overlay highlighted barcodes
      geom_tile(
        data = long_df %>% filter(Barcode %in% highlight_barcodes),
        color = "red", size = 0.25, fill = NA
      ) +
      scale_fill_gradientn(
        colours = rev(brewer.pal(11, "RdBu")),
        name = "log2 fraction"
      ) +
      ##ensure labels are there
      scale_y_continuous(expand = c(0,0), breaks = 1:length(levels(long_df$Outcome_Alias)), labels = levels(long_df$Outcome_Alias)) +
      scale_x_continuous(expand = c(0,0), breaks = 1:length(levels(long_df$Barcode)), labels = levels(long_df$Barcode)) +
      theme_object() +
      coord_cartesian(xlim = ranges_heatmap_brush$x, ylim = ranges_heatmap_brush$y) +
      NULL
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
    message(input$tabs," selected")
  })
}
