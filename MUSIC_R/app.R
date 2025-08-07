source("global.R")
ui <- source("ui.R")$value
server <- source("server.R")$value

shinyApp(ui = ui, server = server)
