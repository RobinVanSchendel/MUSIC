read_in_umap_outcome_data <- function(){
  df = read.table(file = UMAP_OUTCOME_FILE, header = T, sep = "\t", fill = T) %>%
    ##actually needed because something goes wrong with the labels otherwise
    mutate(Type = ifelse(Type == "B-DEL_10p","B-DEL_>10",Type))
  df
}