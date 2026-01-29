generateTabs <- function(data_input) {
  if (data_input == "1") {
    tabsetPanel(id = "tabs",
                tabPanel(
                  "Welcome",
                  h3("Welcome to MUSIC"),
                  img(src = "figure_welcome.jpg", width = "100%"),
                  p("This application provides multiple interactive visualizations for exploring gene-level and barcode-level data."),
                  p("Use the tabs above to navigate between tables and plots, including waterfall, XY, and volcano plots."),
                  p("Hover over plots or select data points where available to view additional details.")
                ),
                tabPanel("Genes & Barcodes", h3("Genes and number of barcodes per gene (genome-wide)"),
                         withSpinner(DT::dataTableOutput("gene_barcode"))),
                tabPanel("Waterfall", h3("Aggregated Z-score for each gene for various outcomes"), uiOutput("UIWaterfall")),
                tabPanel("XY", h3("XY plot is meant to quickly plot if a Gene deviates based on two outcomes types"), uiOutput("UIXYplot"), uiOutput("hover_xy"), DTOutput("xy_table")),
                tabPanel("Volcano", h3("Volcano plot for each outcome type"), uiOutput("UIVolcanoplot"), uiOutput("hover_volcano"), h3("Selected Data Points"), DTOutput("volcano_table")),
                selected = "Welcome"
    )
  } else {
    tabsetPanel(id = "tabs",
                tabPanel("Genes & Barcodes", h3("Genes and the number of barcodes per gene (focused candidate)"), DT::dataTableOutput("gene_barcode")),
                tabPanel("Heatmap",h3("Heat map info..."),uiOutput("UIHeatmapPlot")),
                tabPanel("UMAP - Gene", h3("a UMAP representation of all hits, same as Figure x in our manuscript"),uiOutput("UIUmapPlot"), uiOutput("hover_umap_gene"),
                         DTOutput("umap_gene_table")),
                tabPanel("UMAP - Outcome",h3("Umap outcome plot"),uiOutput("UIUmapOutcomePlot"),
                         h3("Please select gene(s) to show here:"),
                         uiOutput("UIUmapOutcomeGeneSpecificPlot") ),
                tabPanel("Volcano ", uiOutput("UIVolcanoFocusedplot"), uiOutput("hover_volcano_focused"), uiOutput("UIDNAPlot")),
                selected = "UMAP - Outcome"
    )
  }
}
