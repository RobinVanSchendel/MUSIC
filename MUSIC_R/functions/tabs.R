generateTabs <- function(data_input) {
  if (data_input == "1") {
    tabsetPanel(id = "tabs",
                tabPanel("Genes & Barcodes", h1("genes and barcodes"), withSpinner(DT::dataTableOutput("gene_barcode"))),
                tabPanel("Waterfall", uiOutput("test")),
                tabPanel("XY", p("XY plot is meant to quickly ascertain if a Gene deviates..."), uiOutput("UIXYplot")),
                tabPanel("Volcano", uiOutput("UIVolcanoplot"), uiOutput("hover_volcano"))
    )
  } else {
    tabsetPanel(id = "tabs",
                tabPanel("Genes & Barcodes", DT::dataTableOutput("gene_barcode")),
                tabPanel("Heatmap"),
                tabPanel("UMAP - Gene"),
                tabPanel("UMAP - Outcome"),
                tabPanel("Volcano")
    )
  }
}
