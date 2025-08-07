# server.R
function(input, output, session) {
  output$dynamicTabs <- renderUI({
    if(input$data_input == "1"){
      tabsetPanel(id = "tabs",
        tabPanel("Genes & Barcodes",h1("genes and barcodes"),DT::dataTableOutput("gene_barcode")),
        tabPanel("Waterfall", uiOutput("test")),
        tabPanel("XY", p("XY plot is meant to quickly ascertain if a Gene deviates in a certain outcome type. You can highlight the top x genes (up and/or down)
                         and set the number of experiments have to show that deviation. You can also select you favourite gene to see if it deviates."),uiOutput("UIXYplot")),
        tabPanel("Volcano"),
      )
    } else{
      tabsetPanel(id = "tabs",
        tabPanel("Genes & Barcodes",DT::dataTableOutput("gene_barcode")),
        tabPanel("Heatmap"),
        tabPanel("UMAP - Gene"),
        tabPanel("UMAP - Outcome"),
        tabPanel("Volcano"),
      )
      
    }
  })
  
  db_con <- reactive({
    if(input$data_input == "1"){
      con = genome_wide_conn()
    } else if(input$data_input == "2"){
      con = candidate_conn()
    }
  })
  
  output$UIXYplot <- renderUI({
    plotOutput("XYplot", height = input$plotHeight, width = input$plotWidth,
               hover = hoverOpts("plot_hover", delay = 10, delayType = "debounce")
    )
  })
  output$gene_barcode <- DT::renderDataTable({
    req(db_con())
    con = db_con()
    count_df = tbl(con, "countTable") %>% collect()
    
    ##only display gene info:
    gene_info = count_df %>% select(Gene, Barcode) %>% 
      distinct() %>% 
      group_by(Gene) %>% 
      summarise(nrBarcodes = n())
    
    dt = DT::datatable(gene_info, caption = paste("Genes"),
                       rownames = FALSE,extensions = 'Buttons', options = list(
                         lengthMenu = list(c(10, 20, 100 -1), c('10', '20','100', 'All')),
                         pageLength = 20,
                         dom = 'ftlpiB',
                         buttons = c('copy', 'excelHtml5')
                       ))
    dt
    
  })
  
  observe({
    req(db_con())
    con = db_con()
    ##put the genes there
    genes = tbl(con, "barcodeCount") %>% select(Gene) %>% distinct() %>% collect() %>% pull()
    updatePickerInput(session, inputId = "highlightGene", choices = genes)
  })
  
  output$XYplot <- renderPlot({
    plotY = input$XYY
    plotX = input$XYX
    con = db_con()
    df = tbl(con, "geneAlt") %>% filter(Type %in% c(plotY, plotX)) %>% 
      select(Gene, Alias, Type, fraction) %>% collect() %>% 
      spread(key = "Type", value = "fraction")
    
    
    plots = list()  
    plot = ggplot(df,aes(x=!!as.name(plotX),y=!!as.name(plotY)))
    if(input$geom_which == "density_log10"){
      #scale_fill_gradient(low = "grey", high = "black") +
      plot = plot + geom_bin2d(bins = input$densityDots) +scale_fill_viridis_c(trans = "log10")
    } else if(input$geom_which == "density"){
      #scale_fill_gradient(low = "grey", high = "black") +
      plot = plot + geom_bin2d(bins = input$densityDots) +scale_fill_viridis_c()
    }
    else {
      plot = plot + geom_scattermore(alpha = 0.12, pointsize = 2.0)
    }
    
    highlight = input$highlightGene
    if(!is.null(highlight) & length(highlight) > 0){
      highlightPart = df %>% filter(Gene %in% highlight)
      plot = plot + geom_text_repel(data = highlightPart, aes(label=Gene), color = "blue", min.segment.length = unit(0, 'lines')) + 
        geom_point(data = highlightPart, aes(x=!!as.name(plotX),y=!!as.name(plotY)),color="blue", size = 4)
    }
    if(input$highlightTop>0){
      dfTop = df %>% group_by(Alias) %>% slice_max(order_by = !!as.name(plotX),n = input$highlightTop)
      dfTopLow = df %>% group_by(Alias) %>% slice_min(order_by = !!as.name(plotX), n = input$highlightTop)
      
      if(input$overlapping_genes_only){
        dfTopOverlap = dfTop %>% group_by(Gene) %>% count() %>% filter(n >= input$overlapping_genes_count)
        keepTop = dfTopOverlap$Gene
        dfTop = dfTop %>% filter(Gene %in% keepTop)
        
        dfTopLowOverlap = dfTopLow %>% group_by(Gene) %>% count() %>% filter(n>=input$overlapping_genes_count)
        keepTopLow = dfTopLowOverlap$Gene
        dfTopLow = dfTopLow %>% filter(Gene %in% keepTopLow)
      }
      if("up" %in% input$highlightShow){
        plot = plot + geom_point(data = dfTop, color = "red",size =2) + geom_text_repel(data = dfTop, aes(label=Gene), color = "red", min.segment.length = unit(0, 'lines'), max.overlaps = Inf)
      }
      if("down" %in% input$highlightShow){
        plot = plot + geom_point(data = dfTopLow, color = "red",size =2) + geom_text_repel(data = dfTopLow, aes(label=Gene), color = "red", min.segment.length = unit(0, 'lines'), max.overlaps = Inf)
      }
    }
    plot = plot + facet_wrap(Alias ~ ., ncol = 2, labeller = as_labeller(GENOME_WIDE_LABEL_MAP))
    if("horizontal" %in% input$marginalPlot){
      ymin = min(df[[plotY]], na.rm = T)
      ymax = max(df[[plotY]], na.rm = T)
      locY = ymin - (ymax-ymin)/8
      print(ymin)
      print(ymax)
      print(locY)
      plot = plot + geom_boxploth(aes(y = locY), width = (ymax-ymin)/8)
    }
    if("vertical" %in% input$marginalPlot){
      xmin = min(df[[plotX]], na.rm = T)
      xmax = max(df[[plotX]], na.rm = T)
      locX = xmin - (xmax-xmin)/8
      print(locX)
      plot = plot + geom_boxplot(aes(x = locX), width = (xmax-xmin)/8)
    }
    plot + theme_object()
  })
  
  
  ##create a theme_object for all plots
  theme_object <- reactive({
    themeObj = theme()
    themeObj$text = element_text(size = input$theme_font)
    
    if(input$theme_background){
      themeObj$panel.background = element_blank()
    }
    if(input$theme_grid == FALSE){
      themeObj$panel.grid.major = element_blank()
      themeObj$panel.grid.minor = element_blank()
    }
    
    themeObj$axis.line = element_line(colour = 'black', linewidth = input$theme_line_thickness)
    themeObj$axis.ticks = element_line(colour = 'black', linewidth = input$theme_line_thickness)
    
    if(input$theme_text_vertical){
      themeObj$axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
    }
    themeObj
  })
  
  
  observeEvent(input$tabs,{
    if(input$tabs == "XY"){
    }
  })
  
}