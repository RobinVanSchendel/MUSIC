updateZoomRanges <- function(brush_input, ranges) {
  if (!is.null(brush_input)) {
    ranges$x <- c(brush_input$xmin, brush_input$xmax)
    ranges$y <- c(brush_input$ymin, brush_input$ymax)
  } else {
    ranges$x <- NULL
    ranges$y <- NULL
  }
}
