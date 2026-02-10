prepareVolcanoFocusedData <- function(con, type) {
  tbl(con, "outcomes") %>%
    filter(outcomeTop == type) %>%
    collect() %>%
    mutate(log2fraction = log2(fraction/mean))
}


getHighlightedFocusedGenes <- function(df, input) {
  df %>% filter(Gene %in% input$highlightGeneFocused)
}

addHighlightsAndMarginalsVolcanoFocused <- function(plot, input, highlightPart, dfTop) {
  if (!is.null(highlightPart) && nrow(highlightPart) > 0) {
    plot <- plot +
      geom_text_repel(data = highlightPart, aes(label = Gene), color = "blue", min.segment.length = unit(0, 'lines')) +
      geom_point(data = highlightPart, color = "blue", size = 4)
  }
  
  if (input$highlightTop > 0) {
    plot <- plot +
      geom_point(data = dfTop, color = "red", size = 2) +
      geom_text_repel(data = dfTop, aes(label = Gene), color = "red", min.segment.length = unit(0, 'lines'), max.overlaps = Inf)
  }
  plot
}

getDirectionFilter <- function(selected) {
  # Normalize to uppercase for consistency
  selected <- toupper(selected)
  
  if(all(c("UP","DOWN") %in% selected)) {
    return("ALL")
  } else if("UP" %in% selected) {
    return("UP")
  } else if("DOWN" %in% selected) {
    return("DOWN")
  } else {
    # Default fallback if nothing is selected
    return("ALL")
  }
}