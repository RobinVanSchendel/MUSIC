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
        h3("Aggregated P-score for each gene for various outcomes"),
        uiOutput("plot_header1"),
        uiOutput("UIWaterfall"),
        uiOutput("hover_waterfall"),
        h3("Selected Data points:"),
        DTOutput("waterfall_table")
      ),
      
      tabPanel(
        "XY",
        h3("XY plot is meant to directly compare the effect of each gene on two outcome types"),
        p("each dot represents a fraction of mutagenic reads"),
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
        h3("Genes and number of barcodes per gene"),
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
        h3("Heatmap representation of a focused MUSIC"),
        p("For each of the 150 most abundant outcomes (y-axis) – chosen to balance complexity, statistical power, and dataset depth – sgRNAs that deviated by more than three standard deviations across all three biological replicates were selected. From these, we retained the top-100 sgRNAs per outcome showing the strongest positive or negative effects (by Manhattan distance). The mutation redistribution profiles of these sgRNAs (x-axis) are presented in a heatmap, with relative positions determined by unsupervised hierarchical clustering. Colours indicate the log2 fold-change of each outcome frequency relative to non-targeting sgRNAs: depletion is shown in shades of blue, indicating that the targeted gene facilitates the outcome, whereas enrichment is shown in shades of red, indicating a suppressive role. Outcomes are schematically depicted relative to the expected cut site, with each replicate shown separately. Both outcomes and genes are colour-coded according to mutational type and DNA repair pathway, respectively."),
        uiOutput("plot_header4"),
        uiOutput("UIHeatmapPlot")
      ),
      
      tabPanel(
        "UMAP - Gene",
        h3("UMAP representation of focused MUSIC data"),
        p("UMAP projection of the top outlier sgRNAs/genes. Each dot represents an individual sgRNA, positioned based on the similarity of its repair outcome distribution to that of other outlier or non-targeting sgRNAs. Genes are colour-coded by the DSB repair pathways in which they are involved. “53BP1” = 53BP1 sub-pathway of c-NHEJ, “BER/SSBR” = Base Excision Repair/Single-Strand Break Repair, “c-NHEJ” = canonical Non-Homologous End Joining, “FA/ICL” = Fanconi Anaemia Pathway / Interstrand Crosslink Repair, “Fork QC” = Replication Fork Quality Control, “HR” = Homologous Recombination, “MMR” = Mismatch Repair, “NER” = Nucleotide Excision Repair, “RER” = Ribosomal Excision Repair, “TMEJ” = Polymerase Theta-Mediated End Joining."),
        uiOutput("plot_header5"),
        uiOutput("UIUmapPlot"),
        uiOutput("hover_umap_gene"),
        DTOutput("umap_gene_table")
      ),
      
      tabPanel(
        "UMAP - Outcome",
        h3("UMAP-based separation of repair outcomes according to distinct genetic dependencies"),
        p("(A) Outcome-based UMAP projection based on the top 100 outlier sgRNAs, with outcomes plotted for each replicate. Relative positions are a function of the similarity in genetic dependency. Dots are colour coded by mutational type. “SNV”: Single-Nucleotide Variant within 2bp of the cut-site; “HDR”: Homology-Directed Repair; for deletions: “B” – bidirectional, spanning both sides of the break junction; “D” – directional, spanning only one side of the break junction). (B) Same data as in (A), with colour intensity indicating deletion size (in bp). (C) Same data as in (A), with colour intensity reflecting the length of microhomology (in bp). Non-deletion outcomes are represented as open circles."),
        uiOutput("plot_header6"),
        fluidRow(
          column(
            6,
            uiOutput("UIUmapOutcomePlot", height = "600px")
          ),
          column(
            6,
            column(12, uiOutput("UIUmapOutcomeDelSizePlot", height = "300px")),
            column(12, uiOutput("UIUmapOutcomeMHPlot", height = "300px"))
          )
        ),
        uiOutput("hover_umap_outcome"),
        h3("Please select gene(s) to show here:"),
        uiOutput("UIUmapOutcomeGeneSpecificPlot")
      ),
      
      tabPanel(
        "Volcano ",
        h3("Volcano plots for the indicated outcome, showing log2 fold-change (x-axis) versus total number of mutational events (y-axis)"),
        uiOutput("plot_header7"),
        uiOutput("UIVolcanoFocusedplot"),
        uiOutput("hover_volcano_focused"),
        p("Here the actual mutation is schematically drawn:"),
        uiOutput("UIDNAPlot")
      ),
      
      tabPanel(
        "Genes & sgRNAs",
        h3("Genes and the sgRNAs per gene"),
        DT::dataTableOutput("gene_barcode")
      ),
      
      selected = "Welcome"  # Start with the Welcome tab
    )
  }
}
