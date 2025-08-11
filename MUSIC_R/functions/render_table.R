library(DT)

makeInteractiveTable <- function(data, pageLength = -1) {
  datatable(
    data,
    options = list(
      pageLength = pageLength,
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
    ),
    extensions = 'Buttons',
    filter = 'top',
    selection = 'multiple',
    rownames = FALSE
  )
}
