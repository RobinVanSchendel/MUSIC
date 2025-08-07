# server.R
function(input, output, session) {
  output$dynamicTabs <- renderUI({
    if(input$data_input == "1"){
      tabsetPanel(
        tabPanel("Genes & Barcodes",h1("genes and barcodes"),DT::dataTableOutput("gene_barcode")),
        tabPanel("Waterfall", uiOutput("test")),
        tabPanel("XY"),
        tabPanel("Volcano"),
      )
    } else{
      tabsetPanel(
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
  
  output$test <- renderUI({
    plotOutput("testPlot", height = input$plotHeight, width = input$plotWidth,
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
  
  output$testPlot <- renderPlot({
    req(db_con())
    print(db_con())
    ggplot(mtcars, aes(x = mpg, y = qsec)) + geom_point()
  })
  
  
  
}