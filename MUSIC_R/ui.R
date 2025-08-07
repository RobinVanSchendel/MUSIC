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
    ),
    mainPanel(
      uiOutput("dynamicTabs")
    )
  )
)
