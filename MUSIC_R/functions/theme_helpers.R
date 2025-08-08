generateTheme <- function(input) {
  themeObj <- theme(text = element_text(size = input$theme_font))
  
  if (input$theme_background) {
    themeObj <- themeObj + theme(panel.background = element_blank())
  }
  if (!input$theme_grid) {
    themeObj <- themeObj + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  }
  
  themeObj <- themeObj + theme(
    axis.line = element_line(colour = 'black', linewidth = input$theme_line_thickness),
    axis.ticks = element_line(colour = 'black', linewidth = input$theme_line_thickness)
  )
  
  if (input$theme_text_vertical) {
    themeObj <- themeObj + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }
  
  themeObj
}
