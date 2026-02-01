library(gridExtra)
library(gplots)
library(DT)
library(scattermore)
library(ggstance)
library(shinycssloaders)
library(ggExtra)
library(ggrepel)
library(dplyr)
library(tidyr)
library(dbplyr)
library(shinyWidgets)
library(RColorBrewer)
library(RSQLite)
library(stringr)
library(ggnewscale)
library(shinyalert)
library(tibble)
library(ggh4x)
library(ggplot2)
library(shiny)

source("functions/tabs.R")
source("functions/theme_helpers.R")
source("functions/zoom.R")
source("functions/render_table.R")

source("genome_wide_r/xy_plot_helpers.R")
source("genome_wide_r/gene_barcode.R")
source("genome_wide_r/volcano_plot_helpers.R")
source("genome_wide_r/waterfall_plot_helpers.R")

source("focused_candidate_r/volcano_focused_plot_helpers.R")
source("focused_candidate_r/umap_gene_plot_helpers.R")
source("focused_candidate_r/heatmap_plot_helpers.R")
source("focused_candidate_r/plot_dna_event.R")
source("focused_candidate_r/umap_outcome_plot_helpers.R")


GENOME_WIDE_DB = "Z:/Marco/Temp/Backup D/Temp/MBdata/Processed_Again/MBCrisprMBAgain_1_1.db" ##still change
CANDIDATE_DB = "Z:/Marco/CustomClonedLibrary/oPools seq files/MBCrisprMBSubscreeen_Full_1MHfSNVtop150.db"   # new, test if same output as originalis a

UMAP_GENE_FILE = "data/dfTest.txt"
UMAP_OUTCOME_FILE = "data/outcomeUmapcoordinates.txt"
HEATMAP_FILE = "data/20240930-145847_ForUMAP_testNT.txt"
WATERFALL_FILE = "data/TableS1_Genome-wide_MUSIC.txt"

##only used in waterfall plots for now
GENOME_WIDE_TYPES = c( "Insertion (1bp)","Deletion (≥15bp)",
                       "Delins","Homology-directed repair",
                       "POLQ-deletion" ,"Templated insertions")

##ensure the waterfall types are coloured
WATERFALL_COLORS = setNames(c(
  "#FF8C00", "#124E8B","#05878C","#EC71A7", "#B03060", "#E6332A"),
  GENOME_WIDE_TYPES
)

##Pathway colors
WATERFALL_PATHWAY_COLORS = c("unassigned" = "#CC0000", "HIT" = "#CC0000",
                             "BER/SSBR" = "khaki3",  "FORK QC" = "forestgreen",
                             "HR" = "purple",  "c-NHEJ" = "darkorange",       "NER" = "pink"  ,   "FA/ICL" ="hotpink1" ,
                             "53BP1" = "gold2" ,  "MMR"  = 'violet',   
                             "CONTROL" = "blue" , "TMEJ"  = "maroon" ,   "SSA"  = "brown"  ,   "RER" =  "magenta")

##for waterfall plot
WATERFALL_SD_LIMIT = 2.5

GENOME_WIDE_LABEL_MAP <- c(
  MB01 = "Target 1 - Replicate 1",
  MB02 = "Target 1 - Replicate 2",
  MB03 = "Target 2 - Replicate 1",
  MB04 = "Target 2 - Replicate 2",
  MB05 = "Target 3 - Replicate 1",
  MB06 = "Target 3 - Replicate 2"
)

FOCUSED_MAP <- c(
  MBsub1 = "Target 1 - Replicate 1",
  MBsub2 = "Target 1 - Replicate 2",
  MBsub3 = "Target 1 - Replicate 3"
)

##Define colours UMAP - gene
UMAP_GENE_COLORS = c("other" = "black", "HIT" = "red",
          "BER/SSBR" = "khaki",  "FORK QC" = "forestgreen",
          "HR" = "purple",  "EJ" = "orange",       "NER" = "pink"  ,   "FA/ICL" ="hotpink1" ,
          "53BP1" = "gold2" ,  "MMR"  = 'violet',   
          "CONTROL" = "blue" , "TMEJ"  = "maroon" ,   "SSA"  = "brown"  ,   "RER" =  "magenta")

HEATMAP_REPLICATE_COLORS = c("1" = "#BF9136", "2" = "#828330", "3" = "#226162")

UMAP_OUTCOME_COLORS = c(
  "B-DEL_10p" = "palegreen4",
  "B-DEL_>10" = "palegreen4",
  "B-DEL_4-10" = "palegreen3",
  "B-DEL_1-3" = "palegreen2",
  "D-DEL_10p" = "dodgerblue4",
  "D-DEL_4-10" = "dodgerblue3",   
  "D-DEL_1-3" = "dodgerblue",   
  "DELINS" = "darkslategray3",       
  "HDR" = "hotpink1",          
  "INSERTION_1bp" = "orange",
  "INSERTION" = "gold2",   
  "PQ_DELETION" = "maroon", 
  "TINS" = "red",
  "SNV" = "darkslategray4"
)


genome_wide_conn <- function(){
  if(!file.exists(GENOME_WIDE_DB)){
    warning(paste("File does not exist:",GENOME_WIDE_DB))
    return(NULL)
  }
  con <- dbConnect(RSQLite::SQLite(), dbname = GENOME_WIDE_DB)
}

candidate_conn <- function(){
  if(!file.exists(CANDIDATE_DB)){
    warning(paste("File does not exist:",CANDIDATE_DB))
    return(NULL)
  }
  con <- dbConnect(RSQLite::SQLite(), dbname = CANDIDATE_DB)
}

##collect the genes/barcodes to keep the data static
genome_wide_connection = genome_wide_conn()
if(!is.null(genome_wide_connection)){
  GENOME_WIDE_GENE_INFO <- tbl(genome_wide_connection, "countTable") %>% select(Gene, Barcode) %>%
    distinct() %>%
    group_by(Gene) %>%
    summarise(nrBarcodes = n()) %>% 
    collect()
} else{
  GENOME_WIDE_GENE_INFO = NULL
}

##collect the genes/barcodes to keep the data static
candidate_connection = candidate_conn()
if(!is.null(candidate_connection)){
  FOCUSED_GENE_INFO <- tbl(candidate_connection, "countTable") %>% select(Gene, Barcode) %>%
    distinct() %>%
    group_by(Gene) %>%
    summarise(nrBarcodes = n()) %>% 
    collect()
  
  FOCUSED_OUTCOMES_MEAN <- tbl(candidate_connection, "outcomes") %>% 
    select(outcomeTop, mean, Alias) %>%
    distinct() %>%
    group_by(outcomeTop) %>%
    summarise(mean = mean(mean)) %>%
    arrange(desc(mean)) %>%
    collect() 
  
  FOCUSED_OUTCOMES = setNames(FOCUSED_OUTCOMES_MEAN$outcomeTop, 
                              paste(FOCUSED_OUTCOMES_MEAN$outcomeTop, round(FOCUSED_OUTCOMES_MEAN$mean, 5)))
  
} else{
  FOCUSED_GENE_INFO = NULL
  FOCUSED_TYPES = NULL
}

getPlotWidth <- function(width_input) {
  if (is.null(width_input) || width_input == 0) {
    return("100%")
  } else {
    return(width_input)
  }
}




