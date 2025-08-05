## Clean-up the environment before a new run
gc()
rm(list = ls())

# Initialize global, server, and user interface
source("app/global.R", local = TRUE)
source("app/ui.R", local = TRUE)
source("app/server.R", local = TRUE)

# Run app
shinyApp(ui = ui, server = server)
