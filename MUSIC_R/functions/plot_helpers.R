
##function to create a tooltip based on the a hover
createHoverTooltip <- function(
    hover,
    df,
    display_cols = c("Gene", "log2fraction", "fraction", "Pvalue", "counts", "trials"),
    labels = c(
      Gene = "Name",
      log2fraction = "Fold change",
      fraction = "Fraction",
      Pvalue = "Significance",
      counts = "Mutagenic events",
      trials = "Trials"
    ),
    digits = 2
) {
  if (is.null(hover)) return(NULL)
  
  point <- nearPoints(df, hover, threshold = 10, maxpoints = 1, addDist = FALSE)
  if (nrow(point) == 0) return(NULL)
  
  ## Positioning
  left_pct <- (hover$x - hover$domain$left) /
    (hover$domain$right - hover$domain$left)
  top_pct <- (hover$domain$top - hover$y) /
    (hover$domain$top - hover$domain$bottom)
  
  left_px <- hover$range$left +
    left_pct * (hover$range$right - hover$range$left)
  top_px <- hover$range$top +
    top_pct * (hover$range$bottom - hover$range$top)
  
  style <- paste0(
    "position:absolute; padding:5px; z-index:100;",
    "background-color: rgba(200,200,245,0.65); ",
    "left:", left_px + 20, "px; top:", top_px + 32, "px;"
  )
  
  ## Build tooltip content
  lines <- character(0)
  for (col in display_cols) {
    if (!col %in% colnames(point)) next
    
    ##get the value
    value <- point[[col]]
    if (is.numeric(value)) value <- round(value, digits)
    
    ##find the label    
    label <- if (!is.null(labels) && col %in% names(labels)) {
      labels[[col]]
    } else {
      col
    }
    
    ##add the text
    lines <- c(
      lines,
      paste0("<b>", label, ":</b> ", value, "<br/>")
    )
  }
  
  wellPanel(
    style = style,
    p(HTML(paste(lines, collapse = "")))
  )
}