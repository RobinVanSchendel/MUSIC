prepareVolcanoData <- function(con, type) {
  aliases = names(GENOME_WIDE_LABEL_MAP)
  tbl(con, "geneAlt") %>%
    filter(Type == type, Alias %in% aliases ) %>%
    mutate(log2fraction = log(fraction / mean) / log(2)) %>%
    collect() 
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

createHoverTooltip <- function(hover, df) {
  if (is.null(hover)) return(NULL)
  
  point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
  if (nrow(point) == 0) return(NULL)
  
  left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)
  
  left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)
  
  style <- paste0("position:absolute;
                  padding: 5px;
                  z-index:100; background-color: rgba(200, 200, 245, 0.65); ",
                  "left:", left_px + 20, "px; top:", top_px + 32, "px;")
  ##check if log2fraction is present
  if("log2fraction" %in% colnames(point)){
    wellPanel(
      style = style,
      p(HTML(paste0(
        "<b> Name: </b>", point$Gene, "<br/>",
        "<b> Fold change: </b>", round(point$log2fraction, 2), "<br/>",
        "<b> Fraction: </b>", round(point$fraction, 2), "<br/>",
        "<b> Significance: </b>", round(point$Pvalue, 2), "<br/>",
        "<b> Mutagenic events: </b>", point$counts, "/", point$trials, "<br/>"
      )))
    )
  } else{
    ##probably have to check a bit better
    wellPanel(
      style = style,
      p(HTML(paste0(
        "<b> Name: </b>", point$Gene, "<br/>"
      )))
    )
  }
  
}
