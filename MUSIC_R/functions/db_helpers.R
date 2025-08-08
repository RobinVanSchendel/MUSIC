getDbConnection <- function(data_input) {
  if (data_input == "1") {
    genome_wide_conn()
  } else if (data_input == "2") {
    candidate_conn()
  } else {
    NULL
  }
}
