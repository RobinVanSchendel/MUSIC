generateTabs <- function(data_input) {
  tags$head(
    tags$script("
    $(function () {
      $('[data-toggle=\"tooltip\"]').tooltip();
    });
  ")
  )
  if (data_input == "1") {
    tabsetPanel(
      id = "tabs",
      
      # Welcome tab for genome-wide dataset
      tabPanel(
        "Welcome",
        h3("Welcome to MUSIC"),
        p("Using a genome-wide CRISPR library in IB10 mES cells we established a comprehensive catalogue of MUtational Scars of Induced DNA Cleavage (MUSIC)"),
        p("For the genome-wide screen we used 3 different targets (sgTargeting-1, sgTargeting-2 and sgTargeting-3) each with 2 replicates"),
        br(), br(),
        img(src = "figure_welcome.jpg", width = "100%"),
        br(), br(),
        p("This application provides multiple interactive visualizations for exploring gene-level and sgRNA-level data."),
        p("You can start here and use the tabs above to navigate between plots, including Waterfall, XY, and Volcano plots."),
        p("Interactive features:"),
        tags$ul(
          tags$li("Hover over points in plots: display detailed information about the corresponding gene or barcode."),
          tags$li("Brush or select regions in plots: highlight subsets of data, which can be examined in the linked data tables below the plots."),
          tags$li("Tables are linked to the plots: selecting points in a plot will filter or highlight rows in the associated table, and vice versa when applicable.")
        ),
      ),
      
      tabPanel(
        "Waterfall",
        h3("Aggregated Z-score for each gene for various outcomes"),
        uiOutput("plot_header1"),
        uiOutput("UIWaterfall"),
        uiOutput("hover_waterfall"),
        h3("Selected Data points:"),
        DTOutput("waterfall_table")
      ),
      
      tabPanel(
        "XY",
        h3("XY plot is meant to quickly plot if a gene deviates based on two outcome types"),
        uiOutput("plot_header2"),
        uiOutput("UIXYplot"),
        uiOutput("hover_xy"),
        DTOutput("xy_table")
      ),
      
      tabPanel(
        "Volcano",
        h3("Volcano plot for each outcome type"),
        uiOutput("plot_header3"),
        uiOutput("UIVolcanoplot"),
        uiOutput("hover_volcano"),
        h3("Selected Data points:"),
        DTOutput("volcano_table")
      ),
      
      tabPanel(
        "Genes & Barcodes",
        h3("Genes and number of barcodes per gene (genome-wide)"),
        withSpinner(DT::dataTableOutput("gene_barcode"))
      ),
      
      selected = "Welcome"  # Make Welcome the default starting tab
    )
    
  } else {
    tabsetPanel(
      id = "tabs",
      
      # Welcome tab for focused dataset
      tabPanel(
        "Welcome",
        h3("Welcome to MUSIC - Focused Dataset"),
        #img(src = "figure_welcome_focused.jpg", width = "100%"),
        p("The focused candidate screeen was performed using a set of 2372 sgRNAs targeting 743 genes that were identified through the genome-wide screen."),
        p("The focused candidate screen also contained 170 Non-Targeting sgRNAs as controls"),
        p("Some of the plots: Heatmap and UMAP plots only include the 264 sgRNAs that we found to affect mutational outcomes."),
        p("You can start here and navigate using the tabs above, including Heatmap, UMAP, and Volcano plots."),
      ),
      
      tabPanel(
        "Heatmap",
        h3("Heat map showing only significant sgRNAs affecting selected outcomes"),
        uiOutput("plot_header4"),
        uiOutput("UIHeatmapPlot")
      ),
      
      tabPanel(
        "UMAP - Gene",
        h3("UMAP representation of all hits (sgRNAs significantly deviating from Non-Targeting)"),
        uiOutput("plot_header5"),
        uiOutput("UIUmapPlot"),
        uiOutput("hover_umap_gene"),
        DTOutput("umap_gene_table")
      ),
      
      tabPanel(
        "UMAP - Outcome",
        h3("UMAP outcome plot (only outcomes where sgRNAs affect log2 fraction)"),
        uiOutput("plot_header6"),
        uiOutput("UIUmapOutcomePlot"),
        uiOutput("hover_umap_outcome"),
        h3("Please select gene(s) to show here:"),
        uiOutput("UIUmapOutcomeGeneSpecificPlot")
      ),
      
      tabPanel(
        "Volcano ",
        h3("Volcano plots show per Target the deviation (log2, on x-axis) for a single mutational outcome. The y-axis show the number of (mutated) reads per sgRNA"),
        uiOutput("plot_header7"),
        uiOutput("UIVolcanoFocusedplot"),
        uiOutput("hover_volcano_focused"),
        p("Here the actual mutation is schematically drawn:"),
        uiOutput("UIDNAPlot")
      ),
      
      tabPanel(
        "Genes & sgRMAs",
        h3("Genes and the sgRNAs per gene"),
        DT::dataTableOutput("gene_barcode")
      ),
      
      selected = "Volcano "  # Start with the Welcome tab
    )
  }
}
