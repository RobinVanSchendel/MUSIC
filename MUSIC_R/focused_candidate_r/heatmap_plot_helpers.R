read_in_heatmap_data <- function(){
  df = read.table(file = HEATMAP_FILE, header = T, sep = "\t", fill = T) %>%
    mutate(Alias = gsub("MBsub","",Alias))
  df
}

convert_data_frame_to_matrix <- function(df, limits = c(-2,2), na.rm = T){
  ##probably the gene list needs to be included here
  wide_df <- df %>%
    ##remove NonTargeting
    filter(!grepl("NonTargeting", Gene)) %>%
    mutate(Outcome_Alias = paste(Outcome, Alias)) %>%
    select(Outcome_Alias, Barcode, log2fraction) %>%
    pivot_wider(
      names_from  = Barcode,
      values_from = log2fraction
    ) 
  
  # Convert to matrix
  mat <- wide_df %>%
    column_to_rownames("Outcome_Alias") %>%  # set row names
    as.matrix()
  
  #limit values
  mat <- pmax(pmin(mat, limits[2]), limits[1])
  
  # Replace NA with 0
  if(na.rm == T){
    mat[is.na(mat)] <- 0
  }
  
  return(mat)
}

##function that takes the a matrix as an input and returns a
##modified data.frame that is longer and uses heatmap.2 for its sorting
retrieve_long_format_sorted_by_heatmap <- function(matrix){
  
  #create a heatmap.2 from the matrix
  hm = heatmap.2(matrix)
  ##make the input ready for ggplot
  long_df <- as.data.frame(matrix) |>
    tibble::rownames_to_column("Outcome_Alias") |>
    tidyr::pivot_longer(
      -Outcome_Alias,
      names_to  = "Barcode",
      values_to = "log2fraction"
    ) %>%
    mutate(
      Barcode   = factor(Barcode, levels = rownames(hm$carpet)),
      Outcome_Alias = factor(Outcome_Alias, levels = colnames(hm$carpet))
    ) %>%
    mutate(
      x_num = as.numeric(Barcode),
      y_num = as.numeric(Outcome_Alias)
    )
  
  long_df
}