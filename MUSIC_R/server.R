# server.R
function(input, output, session) {
  output$dynamicTabs <- renderUI({
    if(input$data_input == "1"){
      tabsetPanel(
        tabPanel("Waterfall"),
        tabPanel("XY"),
        tabPanel("Volcano"),
      )
    } else{
      tabsetPanel(
        tabPanel("Heatmap"),
        tabPanel("UMAP - Gene"),
        tabPanel("UMAP - Outcome"),
        tabPanel("Volcano"),
      )
      
    }
  })
}