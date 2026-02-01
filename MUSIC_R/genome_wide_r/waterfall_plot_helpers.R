prepareWaterfallData <- function() {
  df = read.table(file = WATERFALL_FILE, header = T, sep = "\t") %>%
    mutate(
      Type = recode(
        Type,
        "Insertion (1bp)"  = "Insertion (1bp)",
        "Deletion (?15bp)" = "Deletion (≥15bp)",
        "POLQ-deletion"    = "POLQ-deletion",
        "HDR" = "Homology-directed repair",
        "TINS" = "Templated insertions"
      )
    ) %>%
    ##set as factor
    mutate(Type = factor(Type, levels = GENOME_WIDE_TYPES)) %>%
    ##mutate EJ here
    mutate(Pathway1 = ifelse(Pathway1 == "EJ", "c-NHEJ", Pathway1))
  df
}

##annotate the waterfall plot
annotate_waterfall_data <- function(df, sdLimit = WATERFALL_SD_LIMIT) {
  
  df %>%
    group_by(Type) %>%
    mutate(
      mean_lfc = mean(meanLFC, na.rm = TRUE),
      sd_lfc   = sd(meanLFC, na.rm = TRUE),
      RankDown = max(RankUp, na.rm = TRUE) - RankUp,
      hit = meanLFC >  mean_lfc + sdLimit * sd_lfc |
        meanLFC <  mean_lfc - sdLimit * sd_lfc |
        RankUp <= 20 |
        RankDown <= 20,
      assigned = Pathway1 != "unassigned"
    ) %>%
    ungroup()
}

waterfall_sd_lines <- function(df) {
  df %>%
    distinct(Type, mean_lfc, sd_lfc)
}