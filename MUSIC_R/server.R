function(input, output, session) {
  
  # Render dynamic tabs based on input
  output$dynamicTabs <- renderUI({
    generateTabs(input$data_input)
  })
  
  ##store some data
  volcanodata <- reactiveVal()
  volcanofocuseddata <- reactiveVal()
  xydata <- reactiveVal()
  annotated_waterfall <- reactiveVal()
  umap_gene_data <- reactiveVal()
  umap_outcome_data <- reactiveVal()
  heatmap_data <- reactiveVal()
  barcode_data_focused <- reactiveVal()
  
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
                           brush = brushOpts(id = "plot_waterfall_brush", resetOnNew = T),
                           dblclick = "plot_waterfall_dblclick"))
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
                           brush = brushOpts(id = "plot_volcano_brush", resetOnNew = T),
                           dblclick = "plot_volcano_dblclick"))
  })
  
  
  output$UIVolcanoFocusedplot <- renderUI({
    withSpinner(plotOutput("VolcanoFocusedplot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_volcano_focused", delay = 10, delayType = "debounce"),
                           brush = brushOpts(id = "plot_volcano_focused_brush", resetOnNew = T),
                           dblclick = "plot_volcano_focused_dblclick"
                           ))
  })
  
  ##umap-gene plot
  output$UIUmapPlot <- renderUI({
    withSpinner(plotOutput("UmapPlot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_umap_focused", delay = 10, delayType = "debounce"),
                           brush = "plot_umap_focused_brush"))
  })
  
  ##umap-outcome plot
  output$UIUmapOutcomePlot <- renderUI({
    withSpinner(plotOutput("UmapOutcomePlot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_umap_outcome_focused", delay = 10, delayType = "debounce"),
                           brush = "plot_umap_outcome_brush"))
  })
  
  
  output$UIUmapOutcomeGeneSpecificPlot <- renderUI({
    withSpinner(plotOutput("UmapOutcomeGeneSpecificPlot", height = input$plotHeight, width = getPlotWidth(input$plotWidth)))
  })
  
  ##heat map plot
  output$UIHeatmapPlot <- renderUI({
    withSpinner(plotOutput("HeatmapPlot", height = input$plotHeight, width = getPlotWidth(input$plotWidth),
                           hover = hoverOpts("plot_heatmap_hover", delay = 10, delayType = "debounce"),
                           brush = brushOpts(id = "plot_heatmap_brush", resetOnNew = T),
                           dblclick = "plot_heatmap_dblclick"))
  })
  
  
  
  ##dna plot
  output$UIDNAPlot <- renderUI({
    withSpinner(plotOutput("DNAPlot", height = 150, width = getPlotWidth(input$plotWidth)))
  })
  
  output$DNAPlot <- renderPlot({
    req(input$VolcanoFocusedType)
    outcome = input$VolcanoFocusedType
    plot_event_dna(outcome)
  })
  
  ##waterfall data
  observeEvent(input$waterfall_type, {
    df <- prepareWaterfallData() %>%
      filter(Type %in% input$waterfall_type)
    df_annotated <- annotate_waterfall_data(
      df,
      sdLimit = WATERFALL_SD_LIMIT   # make this adjustable?
    )
    
    annotated_waterfall(df_annotated)
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
  ranges_volcano_focused_brush <- reactiveValues(x = NULL, y = NULL)
  ranges_waterfall_brush <- reactiveValues(x = NULL, y = NULL)
  
  
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
  #zoom for plot_volcano_focused_brush
  observeEvent(input$plot_volcano_focused_dblclick, {
    updateZoomRanges(input$plot_volcano_focused_brush, ranges_volcano_focused_brush)
  })
  
  # Zoom on double-click waterfall
  observeEvent(input$plot_waterfall_dblclick, {
    updateZoomRanges(input$plot_waterfall_brush, ranges_waterfall_brush)
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
      ##get all barcode info in barcode_data_focused
      barcode_genes = tbl(con, "barcodeCount") %>% select(Gene,Barcode) %>% distinct() %>% collect()
      barcode_data_focused(barcode_genes)
      print("update highlightGeneFocused finished")
    }
  })
  
  observeEvent(input$VolcanoType, {
      req(db_con_genome_wide(), input$tabs)
      message("tabs found is:", input$tabs)
      con <- db_con_genome_wide()
      ##make this data static
      message("prepareVolcanoData Type:",input$VolcanoType)
      volcano_data <- prepareVolcanoData(con, input$VolcanoType)
      message("got prepareVolcanoData")
      volcanodata(volcano_data)
  })
  
  ##Update the selectatble genes based on the tab selected
  observeEvent(input$tabs, {
    message("input tab|",input$tabs,"|")
    ##TODO: this is probably not only for heatmap
    if(input$tabs == "Heatmap"){
      genes = heatmap_data() %>% select(Gene) %>% distinct() %>% pull()
      updatePickerInput(session, inputId = "highlightGeneFocused", choices = genes)
    } else{
      req(barcode_data_focused())
      genes = barcode_data_focused() %>% select(Gene) %>% distinct() %>% pull()
      updatePickerInput(session, inputId = "highlightGeneFocused", choices = genes)
    }
  })
  
  ##for when the genome wide volcano is chosen and the data is not yet loaded
  ##Only happens the first time
  ##Likely need to index the DB table, because it is very slow
  observeEvent(input$tabs, {
    req(db_con_genome_wide(),input$VolcanoType)
    if(input$tabs != "Volcano"){
      return()
    }
    ##do we need to load volcanodata?
    if(is.null(volcanodata())){
      con <- db_con_genome_wide()
      message("prepareVolcanoData genome wide first attempt:",input$VolcanoType)
      volcano_data <- prepareVolcanoData(con, input$VolcanoType)
      message("got prepareVolcanoData")
      volcanodata(volcano_data)
    }
    message("here we go:",input$tabs)
  })
  
  observeEvent(input$VolcanoFocusedType, {
    req(db_con_candidate())
    con <- db_con_candidate()
    message("prepareVolcanoFocusedData: ", input$VolcanoFocusedType)
    ##make this data staic
    volcano_focused_data <- prepareVolcanoFocusedData(con, input$VolcanoFocusedType)
    volcanofocuseddata(volcano_focused_data)
    message("prepareVolcanoFocusedData done")
  })
  
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
    if(is.null(umap_outcome_data())){
      message("Reading in umap_outcome_data")
      outcome_df = read_in_umap_outcome_data()
      umap_outcome_data(outcome_df)
    }
  })
  
  output$Waterfallplot <- renderPlot({
    req(input$waterfall_type)
    df_plot <- annotated_waterfall()
    req(df_plot)
    
    sd_lines <- waterfall_sd_lines(df_plot)
    
    WATERFALL_COLORS <- WATERFALL_COLORS[names(WATERFALL_COLORS) %in% input$waterfall_type]
    
    highlight = df_plot %>% filter(Gene %in% input$highlightGene)
    
    p <- ggplot(df_plot, aes(x = RankUp, y = meanLFC)) +
      
      geom_hline(
        data = sd_lines,
        aes(yintercept = mean_lfc + 2.5 * sd_lfc),
        linetype = "dashed",
        color = "grey30"
      ) +
      geom_hline(
        data = sd_lines,
        aes(yintercept = mean_lfc - 2.5 * sd_lfc),
        linetype = "dashed",
        color = "grey30"
      ) +
      
      geom_point(
        data = df_plot %>% filter(!hit),
        color = "grey70",
        size = 1,
        alpha = 0.4
      ) +
      
      geom_point(
        data = df_plot %>% filter(hit, !assigned),
        aes(color = Pathway1),
        size = 2,
        alpha = 0.4
      ) +
      
      geom_point(
        data = df_plot %>% filter(hit, assigned),
        aes(color = Pathway1),
        size = 2,
        alpha = 0.8
      ) +
      
      geom_text_repel(
        data = df_plot %>% filter(hit, assigned),
        aes(label = Gene, color = Pathway1),
        show.legend = FALSE,
        size = 8,
        max.overlaps = Inf,
        seed = 123,
        force = 0.75,
        force_pull = 5,
        nudge_x = 0.2,
        nudge_y = 0.2
      ) +
      
      scale_color_manual(values = WATERFALL_PATHWAY_COLORS) +
      facet_wrap2(~Type, strip = strip_themed(
        background_x = elem_list_rect(fill = WATERFALL_COLORS, color = "black"),
        text_x = elem_list_text(color = "white", face = "bold")
      )) +
      theme_object()
    
    if(nrow(highlight) > 0){
      p <- p + geom_point(data = highlight, color = "blue", size = 4) +
        geom_text_repel(data = highlight, aes(label = Gene), color = "blue", size = 8)
    }
    
    p + coord_cartesian(xlim = ranges_waterfall_brush$x, ylim = ranges_waterfall_brush$y)
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
    
    df <- volcanofocuseddata()
    highlightPart <- getHighlightedFocusedGenes(df, input)
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
    
    final_plot <- final_plot + coord_cartesian(xlim = ranges_volcano_focused_brush$x, ylim = ranges_volcano_focused_brush$y)
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
      scale_color_manual(values = UMAP_GENE_COLORS)+
      ##remove legends for size and alpha
      guides(size = "none", alpha = "none")+
      theme(legend.position = "top") +
      theme_object() +
      NULL
    final_plot
  })
  
  ##UMap outcome plot
  output$UmapOutcomePlot <- renderPlot({
    req(umap_outcome_data())
    message("UmapOutcomePlot")
    df = umap_outcome_data()
    
    left_plot <- ggplot(df) +
      geom_point(aes(x = X2, y = X1, fill = Type), size = 5, alpha = 0.3, shape = 21) +
      scale_fill_manual(values = UMAP_OUTCOME_COLORS) +
      theme_object() +
      NULL
    
    top_r_plot <- ggplot(df) +
      geom_point(aes(x = X2, y = X1, color = DelSize), size = 5, alpha = 0.3) +
      scale_color_gradient(low = "grey95", high = "steelblue4", limits=c(0, 40), oob = scales::squish) +
      theme_object() +
      NULL
    
    bottom_r_plot <- ggplot(df) +
      geom_point(aes(x = X2, y = X1, color = MHDel), size = 5, alpha = 0.3) +
      scale_color_gradient(low = "grey95", high = "maroon4", limits=c(0, 4), oob = scales::squish) +
      theme_object() +
      NULL
    
    right_col <- grid.arrange(
      top_r_plot,
      bottom_r_plot,
      ncol = 1
    )
    ##make the complete arrangement
    grid.arrange(
      left_plot,
      right_col,
      ncol = 2,
      widths = c(2, 1)  # left plot wider
    )
  })
  
  ##reactive for umap outcome specific plot
  final_umap_gene_specific_data <- reactive({
    ##only when this tab is selected
    req(input$tabs == "UMAP - Outcome")
    req(
      umap_outcome_data(),
      db_con_candidate(),
      input$highlightGeneFocused
    )
    
    message("retrieving data for umap_gene_specific")
    hm_data <- read_in_umap_gene_specific_data(
      db_con_candidate(),
      input$highlightGeneFocused
    )
    message("retrieved data for umap_gene_specific ", nrow(hm_data))
    
    build_umap_gene_specific_data(
      umap_df = umap_outcome_data(),
      hm_data = hm_data,
      genes    = input$highlightGeneFocused
    )
  })
  
  ##outcome specific plot
  output$UmapOutcomeGeneSpecificPlot <- renderPlot({
    req(final_umap_gene_specific_data())
   
    final_data = final_umap_gene_specific_data()
    
    ggplot(final_data, aes(x = X2, y = X1, color = Pvalue)) +
      geom_point(data = final_data %>% filter(is.na(log2fraction)), aes(shape="NA"), shape = 21, size = 1.5, fill = "white", alpha = 0.3) +
      geom_point(alpha = 0.3, shape=21, aes_string(col = "log2fraction", size = "Pvalue", fill = "log2fraction")) +
      scale_fill_gradient2('GeneLog2fc', na.value="white", low = "navy", mid = "grey95", high = "red", midpoint = 0, limits=c(-2, 2), oob = scales::squish)+ 
      scale_colour_gradient2('GeneLog2fc', na.value="black", low = "navy", mid = "grey95", high = "red", midpoint = 0, limits=c(-2, 2), oob = scales::squish)+ #log2fc
      scale_size_area(
        max_size = 12,
        breaks = c(1,2.5,7.5,25,75),
        labels = c("1","2.5","7.5","25","75+"),
        guide = "legend",
        limits = c(1,75),
        oob = scales::squish
      ) +
      facet_wrap(~Gene, ncol = 3) +
      theme_object() +
      ggtitle("Umap Outcome Specific Gene(s)") +
      NULL
  })
  
  ##heatmap plot
  output$HeatmapPlot <- renderPlot({
    req(heatmap_data())
    heatmap = heatmap_data()
    
    ##convert to a matrix first, takes care of NA and limits (-2, 2)
    matrix = convert_data_frame_to_matrix(heatmap)
    
    if(!is.null(input$highlightGeneFocused)){
      ##now get the Barcodes
      ##get the selected barcodes based on the gene
      highlight_barcodes = heatmap %>% filter(Gene %in% input$highlightGeneFocused) %>% 
        select(Barcode) %>% distinct() %>% pull()
    }
    
    ##retrieve a df in long format, but uses heatmap.2 for its ordering
    ##currently uses integers to plot as that allows for additional components to be added to the plot
    ##e.g. in (-10 - -1)
    long_df = retrieve_long_format_sorted_by_heatmap(matrix)
    
    ##get the replicat information in as well
    replicate_info = heatmap %>% select(Outcome, Alias) %>% distinct() %>%
      mutate(Outcome_Alias = paste(Outcome, Alias))
    replicate_info = dplyr::left_join(long_df %>% select(Outcome_Alias, y_num) %>% distinct(), replicate_info, by = "Outcome_Alias") %>%
      select(Alias, y_num) %>%
      ##set at position -1
      mutate(x_num = -1) %>%
      ##name it Replicate, while leaving Alias intact
      mutate(Replicate = Alias)
    
    nr_x_labels = 60
    if(input$plotWidth > 0){
      nr_x_labels = round(input$plotWidth/15)
    }
    nr_y_labels = round(input$plotHeight/15)
    
    plot <- ggplot(long_df, aes(x=x_num, y = y_num, fill = log2fraction)) +
      geom_tile() +
      
      scale_fill_gradientn(
        colours = rev(brewer.pal(11, "RdBu")),
        name = "log2 fraction"
      ) +
      ##separate scale for these tiles
      ggnewscale::new_scale_fill() +
      geom_tile(data = replicate_info, aes(x=x_num, y = y_num, fill = Replicate), inherit.aes = F) +
      scale_fill_manual(values = HEATMAP_REPLICATE_COLORS) +
      ##ensure labels are there
      scale_y_continuous(
        expand = c(0, 0),
        
        breaks = function(lims) {
          n <- length(levels(long_df$Outcome_Alias))
          
          lo <- max(1, ceiling(lims[1]))
          hi <- min(n, floor(lims[2]))
          
          # target ~15 labels
          step <- max(1, ceiling((hi - lo + 1) / nr_y_labels))
          
          seq(lo, hi, by = step)
        },
        
        labels = function(y) {
          levels(long_df$Outcome_Alias)[y]
        }
      ) +
      #scale_y_continuous(expand = c(0,0), breaks = 1:length(levels(long_df$Outcome_Alias)), labels = levels(long_df$Outcome_Alias)) +
      #scale_x_continuous(expand = c(0,0), breaks = 1:length(levels(long_df$Barcode)), labels = levels(long_df$Barcode)) +
      #scale_x_continuous(expand = c(0,0), labels = levels(long_df$Barcode)) +
      scale_x_continuous(
        expand = c(0, 0),
        
        breaks = function(lims) {
          n <- length(levels(long_df$Barcode))
          
          # visible range
          lo <- max(1, ceiling(lims[1]))
          hi <- min(n, floor(lims[2]))
          
          # target ~15 labels max
          step <- max(1, ceiling((hi - lo + 1) / nr_x_labels))
          
          seq(lo, hi, by = step)
        },
        
        labels = function(x) {
          levels(long_df$Barcode)[x]
        }
      ) +
      theme_object() +
      coord_cartesian(xlim = ranges_heatmap_brush$x, ylim = ranges_heatmap_brush$y) +
      NULL
    
    if(!is.null(input$highlightGeneFocused) && length(highlight_barcodes) > 0){
      ##get the highlight boxes
      barcode_boxes <- long_df %>%
        filter(Barcode %in% highlight_barcodes) %>%
        group_by(Barcode) %>%
        summarise(
          xmin = min(x_num) - 0.5,
          xmax = max(x_num) + 0.5,
          ymin = min(y_num) - 0.5,
          ymax = max(y_num) + 0.5,
          .groups = "drop"
        )
      # Overlay highlighted barcodes
      plot <- plot + 
        geom_rect(
          data = barcode_boxes,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = NA,
          color = "red",
          linewidth = 0.25,
          inherit.aes = F
        )
    }
    
    plot
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
