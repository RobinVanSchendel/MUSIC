# ui.R
fluidPage(
  titlePanel("MUSIC app - Version 1.0"),
  sidebarLayout(
    sidebarPanel(
                 radioButtons("data_input", "Select MUSIC data set:",
                   choices = 
                     list("Genome Wide" = 1,
                          "Focused Candidate" = 2
                     )
                   ,),
    ),
    mainPanel(
      uiOutput("dynamicTabs")
    )
  )
)
