# ui.R
fluidPage(
  titlePanel("MUSIC app - Version 1.0"),
  sidebarLayout(
    sidebarPanel(
                 radioButtons("data_input", "Select MUSIC data set:",
                   choices = list("Genome Wide" = 1,"Focused Candidate" = 2)),
                 hr(),
                 sliderInput("plotHeight", "Plot height (# pixels): ",
                             value = 600, min = 100, max = 5000, step = 50
                 ),
                 sliderInput("plotWidth", "Plot width (# pixels):", 
                             value = 800, min = 100, max = 5000, step = 50
                 ),
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
                 ##XY plot
                 conditionalPanel(
                   condition = "input.tabs == 'XY'",
                   pickerInput(
                     inputId = "XYX", 
                     label = "Select outcome type (X-axis):",
                     choices = GENOME_WIDE_TYPES,
                     selected = "PQ_DELETION",
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
                                value = 200),
                   checkboxGroupInput(inputId = "marginalPlot",
                                      label = "Add boxplot:",
                                      choices = c("horizontal","vertical"),
                                      inline = T
                   ),
                   numericInput("highlightTop", "Highlight top 'x' genes (based on x-axis only!): ", value = 100),
                   checkboxGroupInput(inputId = "highlightShow", label = "Show  top up/down:", inline = T, choices = c("down","up"), selected =  c("down","up")),
                   checkboxInput(inputId = "overlapping_genes_only",
                                 label = "Only show genes present in top 'x' in multiple experiments:",
                                 value = TRUE),
                   numericInput("overlapping_genes_count", "Gene in top 'x' in how many experiments:", value = 4),
                 ),
    ),
    mainPanel(
      uiOutput("dynamicTabs")
    )
  )
)
