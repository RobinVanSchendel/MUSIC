prepareWaterfallData <- function(con, input) {
  
  
  tbl(con, "geneAlt") %>%
    ##TODO: alter this entire calculation!!
    filter(Alias == "MB01") %>%
    filter(Type %in% input$waterfall_type) %>%
    select(Gene, Alias, Type, fraction) %>%
    collect()
}