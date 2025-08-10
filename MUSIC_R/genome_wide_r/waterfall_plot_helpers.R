prepareWaterfallData <- function(con) {
  tbl(con, "geneAlt") %>%
    ##TODO: alter this entire calculation!!
    filter(Alias == "MB01") %>%
    filter(Type %in% GENOME_WIDE_TYPES) %>%
    select(Gene, Alias, Type, fraction) %>%
    collect()
}