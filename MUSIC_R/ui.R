# ui.R
fluidPage(
  titlePanel("MUSIC app - Version 1.0"),
  sidebarLayout(
    sidebarPanel(
                 radioButtons("data_input", "Select MUSIC data set:",
                   choices = list("Genome Wide" = 1,"Focused Candidate" = 2),
                   selected = 1),
                 hr(),
                 sliderInput("plotHeight", "Plot height (# pixels): ",
                             value = 600, min = 0, max = 5000, step = 50
                 ),
                 sliderInput("plotWidth", "Plot width (# pixels, 0 = 100%):", 
                             value = 0, min = 0, max = 5000, step = 50
                 ),
                 ##Gene for Genome wide
                 conditionalPanel(
                   condition = "input.data_input == '1'",
                   pickerInput(
                     inputId = "highlightGene", 
                     label = "Select Gene(s) to highlight:",
                     choices = NULL,
                     multiple = T,
                     options = pickerOptions(
                       actionsBox = TRUE,
                       `live-search` = TRUE,
                       virtualScroll = 10,
                       size = 10
                     )
                     
                   )
                 ),
                 ##Gene for Genome wide
                 conditionalPanel(
                   condition = "input.data_input == '2'",
                   pickerInput(
                     inputId = "highlightGeneFocused", 
                     label = "Select Gene(s) to highlight:",
                     choices = NULL,
                     multiple = T,
                     options = pickerOptions(
                       actionsBox = TRUE,
                       `live-search` = TRUE,
                       virtualScroll = 10,
                       size = 10
                     )
                     
                   )
                 ),
                 
                 dropdownButton(
                   label = "Theme settings",
                   tooltip = T,
                   circle = F,
                   icon = icon("gears"),
                   checkboxInput(inputId = "theme_text_vertical",
                                 label = "X-axis text vertical",
                                 value = T),
                   checkboxInput(inputId = "theme_background",
                                 label = "theme background blank",
                                 value = F),
                   checkboxInput(inputId = "theme_grid",
                                 label = "add gridlines",
                                 value = T),
                   numericInput(inputId = "theme_font",
                                label = "font size",
                                value = 16),
                   numericInput(inputId = "theme_line_thickness",
                                label = "line thickness",
                                value = 0.25)
                 ),
                 hr(),
                 ###Waterfall plot
                 conditionalPanel(
                   condition = "input.tabs == 'Waterfall'",
                   pickerInput(
                     inputId = "waterfall_type", 
                     label = "Select outcome type(s):",
                     choices = GENOME_WIDE_TYPES,
                     selected = GENOME_WIDE_TYPES,
                     multiple = T,
                     options = pickerOptions(
                       actionsBox = TRUE,
                       `live-search` = TRUE,
                       virtualScroll = 10,
                       size = 10
                     )
                   )
                 ),
                 
                 ##XY plot
                 conditionalPanel(
                   condition = "input.tabs == 'XY'",
                   pickerInput(
                     inputId = "XYX", 
                     label = "Select outcome type (x-axis):",
                     choices = GENOME_WIDE_TYPES,
                     selected = "POLQ-deletion",
                     multiple = F,
                     options = pickerOptions(
                       actionsBox = TRUE,
                       `live-search` = TRUE,
                       virtualScroll = 10,
                       size = 10
                     )
                   ),
                   pickerInput(
                     inputId = "XYY", 
                     label = "Select outcome type (y-axis):",
                     selected = GENOME_WIDE_TYPES[4], ##HDR
                     choices = GENOME_WIDE_TYPES,
                     multiple = F,
                     options = pickerOptions(
                       actionsBox = TRUE,
                       `live-search` = TRUE,
                       virtualScroll = 10,
                       size = 10
                     )
                   ),
                   radioButtons(inputId = "geom_which",
                                label = "Visualize dots:",
                                choices = c("density", "density_log10" , "dots"),
                                selected = "density",
                                inline = T),
                   numericInput(inputId = "densityDots",
                                label = "Bin size for density plot",
                                value = 500),
                   checkboxGroupInput(inputId = "marginalPlot",
                                      label = "Add boxplot:",
                                      choices = c("horizontal","vertical"),
                                      inline = T
                   ),
                 ),
                 ##volcano genome-wide
                 conditionalPanel(
                   condition = "input.tabs == 'Volcano'",
                   pickerInput(
                     inputId = "VolcanoType",
                     label = "Select Type(s):",
                     choices = GENOME_WIDE_TYPES,
                     selected = GENOME_WIDE_TYPES[1],
                     options = list(
                       `actions-box` = TRUE, 
                       size = 10,
                       `live-search`=TRUE,
                       `selected-text-format` = "count > 3"
                     ), 
                     multiple = F
                   ),
                 ),
                 ##volcano genome-wide
                 conditionalPanel(
                   condition = "input.tabs == 'Volcano '",
                   pickerInput(
                     inputId = "VolcanoFocusedType",
                     label = "Select Outcome(s):",
                     choices = FOCUSED_OUTCOMES,
                     selected = FOCUSED_OUTCOMES[2],
                     options = list(
                       `actions-box` = TRUE, 
                       size = 10,
                       `live-search`=TRUE,
                       `selected-text-format` = "count > 3"
                     ), 
                     multiple = F
                   ),
                 ),
                 ##UMAP - gene
                 conditionalPanel(
                   condition = "input.tabs == 'UMAP - Gene'",
                   checkboxInput(inputId = "UMAP_gene_add_labels", label = "add Gene labels", value = T)
                 ),
                 wellPanel(
                   h4("Highlight Top Genes"),
                   numericInput("highlightTop", "Highlight top N genes: ", value = 10),
                   prettyCheckboxGroup(inputId = "highlightShow", label = "Direction of top genes to highlight:", inline = T, 
                                      choices = c("Reduced" = "down", "Increased" = "up"), selected =  c("down","up"), status = "primary"),
                 ),
                 
                 # Section 2: Overlap Filtering (conditionally shown)
                 wellPanel(
                   h4("Filter by Overlap Across Experiments"),
                   materialSwitch(inputId = "overlapping_genes_only",
                                 label = "Enable filter:",
                                 value = TRUE, status = "primary"),
                   numericInput("overlapping_genes_count", "Minimum number of experiments a gene must appear in top N:", value = 4),
                 )
    ),
    mainPanel(
      uiOutput("dynamicTabs")
    )
  )
)
