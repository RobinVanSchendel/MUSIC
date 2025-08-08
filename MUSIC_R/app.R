source("global.R")
ui <- source("ui.R")$value
server <- source("server.R")$value

##ensure to disable this in the production environment
devtools::load_all(reload = TRUE)

shinyApp(ui = ui, server = server)