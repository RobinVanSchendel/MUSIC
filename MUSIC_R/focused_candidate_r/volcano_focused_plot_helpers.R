prepareVolcanoFocusedData <- function(con, type) {
  tbl(con, "outcomes") %>%
    filter(outcomeTop == type) %>%
    collect() %>%
    mutate(log2fraction = log2(fraction/mean))
}

getTopGenesVolcano <- function(df, input, plotX) {
  ##at the moment on log2fraction
  dfTop <- df %>% group_by(Alias) %>% slice_max(order_by = log2fraction, n = input$highlightTop)
  dfTopLow <- df %>% group_by(Alias) %>% slice_min(order_by = log2fraction, n = input$highlightTop)
  
  if (input$overlapping_genes_only) {
    dfTopOverlap <- dfTop %>% group_by(Gene) %>% count() %>% filter(n >= input$overlapping_genes_count)
    dfTop <- dfTop %>% filter(Gene %in% dfTopOverlap$Gene)
    
    dfTopLowOverlap <- dfTopLow %>% group_by(Gene) %>% count() %>% filter(n >= input$overlapping_genes_count)
    dfTopLow <- dfTopLow %>% filter(Gene %in% dfTopLowOverlap$Gene)
  }
  
  list(up = dfTop, down = dfTopLow)
}

addHighlightsAndMarginalsVolcano <- function(plot, input, highlightPart, dfTop, dfTopLow) {
  if (!is.null(highlightPart) && nrow(highlightPart) > 0) {
    plot <- plot +
      geom_text_repel(data = highlightPart, aes(label = Gene), color = "blue", min.segment.length = unit(0, 'lines')) +
      geom_point(data = highlightPart, color = "blue", size = 4)
  }
  
  if (input$highlightTop > 0) {
    if ("up" %in% input$highlightShow) {
      plot <- plot +
        geom_point(data = dfTop, color = "red", size = 2) +
        geom_text_repel(data = dfTop, aes(label = Gene), color = "red", min.segment.length = unit(0, 'lines'), max.overlaps = Inf)
    }
    if ("down" %in% input$highlightShow) {
      plot <- plot +
        geom_point(data = dfTopLow, color = "red", size = 2) +
        geom_text_repel(data = dfTopLow, aes(label = Gene), color = "red", min.segment.length = unit(0, 'lines'), max.overlaps = Inf)
    }
  }
  plot
}