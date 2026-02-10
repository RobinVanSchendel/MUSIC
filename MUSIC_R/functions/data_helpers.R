getFoldRatesOutcome <- function(df, sd_mult = 3, return_summary = TRUE) {
  required_cols <- c("Alias", "outcomeTop", "Gene", "log2fraction", "counts")
  
  # Check for required columns
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check that at least one NonTargeting exists
  if (!any(grepl("NonTargeting", df$Gene))) {
    stop("No 'NonTargeting' entries found in the 'Gene' column.")
  }
  
  # Compute summary table
  df_summary <- df %>%
    ungroup() %>%
    group_by(Alias, outcomeTop) %>%
    summarise(
      foldRate1 = -sd_mult * sd(log2fraction[grepl("NonTargeting", Gene) & counts > 0], na.rm = TRUE),
      foldRate2 =  sd_mult * sd(log2fraction[grepl("NonTargeting", Gene) & counts > 0], na.rm = TRUE),
      .groups = "drop"
    )
  
  if (return_summary) {
    return(df_summary)
  } else {
    # Attach foldRates to every row
    df_out <- df %>%
      left_join(df_summary, by = c("Alias", "outcomeTop"))
    return(df_out)
  }
}



##get the P Score that we want to plot for
getPScoreFromPValue <- function(df, PValue = 0.001, two.sided = TRUE, column = "Gene"){
  n_genes <- n_distinct(df[[column]])
  factor <- ifelse(two.sided, 2, 1)
  p_score <- -log10(PValue / n_genes / factor)
  return(p_score)
}


###important function to determine the hits
getTopGenesVolcanoFocused <- function(df,
                                      highlightTop = 5,
                                      overlapping_barcodes_only = FALSE,
                                      overlapping_barcodes_count = 2,
                                      plotX = "log2fraction",
                                      pvalue_col = "PValue",
                                      outcome_col = "Outcome",
                                      direction_filter = c("ALL", "UP", "DOWN")) {
  
  # Match argument
  direction_filter <- match.arg(direction_filter)
  
  # Check required columns
  required_cols <- c(plotX, pvalue_col, "Alias", "Barcode", "Gene", outcome_col)
  missing_cols <- setdiff(required_cols, colnames(df))
  if(length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  
  # Pre-filtering: compute fold rates and threshold P-value
  calculatedPScore <- getPScoreFromPValue(df)
  df <- getFoldRatesOutcome(df, return_summary = F)
  
  # Filter events by P-value and foldRate thresholds
  df <- df %>%
    filter(.data[[pvalue_col]] > calculatedPScore)
  
  if(direction_filter == "ALL"){
    df <- df %>% filter(.data[[plotX]] < foldRate1 | .data[[plotX]] > foldRate2)
  } else if(direction_filter == "UP"){
    df <- df %>% filter( .data[[plotX]] > foldRate2)
  } else if(direction_filter == "DOWN"){
    df <- df %>% filter( .data[[plotX]] < foldRate1)
  }
  
  # Compute Manhattan distance for volcano plot ranking
  df <- df %>%
    mutate(
      manhattanDist = abs(.data[[plotX]]) + abs(.data[[pvalue_col]])
    )
  
  # Strongest UP hits by Manhattan distance
  dfTop <- df %>%
    filter(.data[[plotX]] > 0) %>%
    group_by(Alias) %>%
    slice_max(order_by = manhattanDist, n = highlightTop)
  
  # Strongest DOWN hits by Manhattan distance
  dfTopLow <- df %>%
    filter(.data[[plotX]] < 0) %>%
    group_by(Alias) %>%
    slice_max(order_by = manhattanDist, n = highlightTop)
  
  
  # Optionally filter by overlapping Barcodes
  if (overlapping_barcodes_only) {
    dfTopOverlap <- dfTop %>% group_by(Barcode) %>% count() %>% filter(n >= overlapping_barcodes_count)
    dfTop <- dfTop %>% filter(Barcode %in% dfTopOverlap$Barcode)
    
    dfTopLowOverlap <- dfTopLow %>% group_by(Barcode) %>% count() %>% filter(n >= overlapping_barcodes_count)
    dfTopLow <- dfTopLow %>% filter(Barcode %in% dfTopLowOverlap$Barcode)
  }
  
  # Combine up and down hits
  df_out <- bind_rows(dfTop, dfTopLow) %>% ungroup()
  
  # Assign direction and filter by overlapping genes + direction
  if(nrow(df_out) > 0) {
    df_keep <- df_out %>%
      mutate(direction = ifelse(.data[[plotX]] > 0, "UP", "DOWN")) %>%
      group_by(Barcode, .data[[outcome_col]], direction) %>%
      dplyr::count() %>%
      filter(n >= overlapping_barcodes_count)
    
    df_out <- df_out %>%
      filter(Barcode %in% df_keep$Barcode) %>%
      mutate(direction = ifelse(.data[[plotX]] > 0, "UP", "DOWN")) %>%
      ungroup()
  }
  
  # Apply optional direction filter
  if(direction_filter != "ALL") {
    df_out <- df_out %>% filter(direction == direction_filter)
  }
  
  return(df_out)
}

