generateTabs <- function(data_input) {
  if (data_input == "1") {
    tabsetPanel(id = "tabs",
                tabPanel(
                  "Welcome",
                  h1("Welcome"),
                  p("This application provides multiple interactive visualizations for exploring gene-level and barcode-level data."),
                  p("Use the tabs above to navigate between tables and plots, including waterfall, XY, and volcano plots."),
                  p("Hover over plots or select data points where available to view additional details.")
                ),
                tabPanel("Genes & Barcodes", h1("genes and barcodes"), withSpinner(DT::dataTableOutput("gene_barcode"))),
                tabPanel("Waterfall", uiOutput("UIWaterfall")),
                tabPanel("XY", p("XY plot is meant to quickly ascertain if a Gene deviates..."), uiOutput("UIXYplot"), uiOutput("hover_xy")),
                tabPanel("Volcano", uiOutput("UIVolcanoplot"), uiOutput("hover_volcano"), h3("Selected Data Points"), DTOutput("volcano_table"))
    )
  } else {
    tabsetPanel(id = "tabs",
                tabPanel("Genes & Barcodes", DT::dataTableOutput("gene_barcode")),
                tabPanel("Heatmap"),
                tabPanel("UMAP - Gene"),
                tabPanel("UMAP - Outcome"),
                tabPanel("Volcano ", uiOutput("UIVolcanoFocusedplot"), uiOutput("hover_volcano_focused"))
    )
  }
}
