read_in_umap_gene_data <- function(){
  df = read.table(file = UMAP_GENE_FILE, header = T, sep = "\t", fill = T)
  df
}

##add labels if that is requested
addLabelsUMAP <- function(base_plot, df, input){
  if(input$UMAP_gene_add_labels){
    base_plot = base_plot + geom_text_repel(data=df %>% filter(Gene != "NonTargeting"), 
                                            aes(label = Gene), max.overlaps = 100)
  }
  base_plot
}