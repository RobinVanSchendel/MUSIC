read_in_umap_outcome_data <- function(){
  df = read.table(file = UMAP_OUTCOME_FILE, header = T, sep = "\t", fill = T) %>%
    ##actually needed because something goes wrong with the labels otherwise
    mutate(Type = ifelse(Type == "B-DEL_10p","B-DEL_>10",Type))
  df
}

read_in_umap_gene_specific_data <- function(con, genes){
  tbl(con, "outcomesGene") %>%
    filter(Gene %in% genes) %>%
    collect() %>%
    mutate(log2fraction = log2(fraction/mean),
           OutcomeAlias = paste(outcomeTop, Alias, sep = "|"))
}

##this is to get the correct data for the selected genes
build_umap_gene_specific_data <- function(
    umap_df,
    hm_data,
    genes
) {
  expand_grid(
    Gene = genes,
    OutcomeAlias = unique(umap_df$OutcomeAlias)
  ) %>%
    left_join(
      umap_df %>% select(X1, X2, OutcomeAlias),
      by = "OutcomeAlias"
    ) %>%
    left_join(
      hm_data %>% select(Gene, OutcomeAlias, Pvalue, log2fraction),
      by = c("Gene", "OutcomeAlias")
    ) 
}
