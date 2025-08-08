renderGeneBarcodeTable <- function(data_input) {
  
  if(data_input == "1"){
    data = GENOME_WIDE_GENE_INFO
  } else if(data_input == "2"){
    data = FOCUSED_GENE_INFO
  }
  
  DT::datatable(data, caption = "Genes",
                rownames = FALSE, extensions = 'Buttons',
                options = list(
                  lengthMenu = list(c(10, 20, 100 - 1), c('10', '20', '100', 'All')),
                  pageLength = 20,
                  dom = 'ftlpiB',
                  buttons = c('copy', 'excelHtml5')
                )
  )
}
