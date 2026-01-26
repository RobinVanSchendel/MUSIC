prepareXYData <- function(con, plotX, plotY) {
  aliases = names(GENOME_WIDE_LABEL_MAP)
  tbl(con, "geneAlt") %>%
    filter(Type %in% c(plotY, plotX), Alias %in% aliases) %>%
    select(Gene, Alias, Type, fraction) %>%
    collect() %>%
    tidyr::spread(key = "Type", value = "fraction")
}

getHighlightedGenes <- function(df, input) {
  df %>% filter(Gene %in% input$highlightGene)
}

getTopGenes <- function(df, input, plotX) {
  dfTop <- df %>% group_by(Alias) %>% slice_max(order_by = !!as.name(plotX), n = input$highlightTop)
  dfTopLow <- df %>% group_by(Alias) %>% slice_min(order_by = !!as.name(plotX), n = input$highlightTop)
  
  if (input$overlapping_genes_only) {
    dfTopOverlap <- dfTop %>% group_by(Gene) %>% count() %>% filter(n >= input$overlapping_genes_count)
    dfTop <- dfTop %>% filter(Gene %in% dfTopOverlap$Gene)
    
    dfTopLowOverlap <- dfTopLow %>% group_by(Gene) %>% count() %>% filter(n >= input$overlapping_genes_count)
    dfTopLow <- dfTopLow %>% filter(Gene %in% dfTopLowOverlap$Gene)
  }
  
  list(up = dfTop, down = dfTopLow)
}

buildXYPlot <- function(df, input, plotX, plotY) {
  plot <- ggplot(df, aes(x = !!as.name(plotX), y = !!as.name(plotY)))
  
  if (input$geom_which == "density_log10") {
    plot <- plot + geom_bin2d(bins = input$densityDots) + scale_fill_viridis_c(trans = "log10")
  } else if (input$geom_which == "density") {
    plot <- plot + geom_bin2d(bins = input$densityDots) + scale_fill_viridis_c()
  } else {
    plot <- plot + geom_scattermore(alpha = 0.12, pointsize = 2.0)
  }
  
  plot
}

addHighlightsAndMarginals <- function(plot, df, input, plotX, plotY, highlightPart, dfTop, dfTopLow) {
  if (!is.null(highlightPart) && nrow(highlightPart) > 0) {
    plot <- plot +
      geom_text_repel(data = highlightPart, aes(label = Gene), color = "blue", min.segment.length = unit(0, 'lines')) +
      geom_point(data = highlightPart, aes(x = !!as.name(plotX), y = !!as.name(plotY)), color = "blue", size = 4)
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
  
  if ("horizontal" %in% input$marginalPlot) {
    ymin <- min(df[[plotY]], na.rm = TRUE)
    ymax <- max(df[[plotY]], na.rm = TRUE)
    locY <- ymin - (ymax - ymin) / 8
    plot <- plot + geom_boxploth(aes(y = locY), width = (ymax - ymin) / 8)
  }
  
  if ("vertical" %in% input$marginalPlot) {
    xmin <- min(df[[plotX]], na.rm = TRUE)
    xmax <- max(df[[plotX]], na.rm = TRUE)
    locX <- xmin - (xmax - xmin) / 8
    plot <- plot + geom_boxplot(aes(x = locX), width = (xmax - xmin) / 8)
  }
  
  plot
}
