library(gridExtra)
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
library(RSQLite)
library(ggplot2)
library(shiny)

source("functions/tabs.R")
source("functions/db_helpers.R")
source("functions/theme_helpers.R")
source("genome_wide_r/xy_plot_helpers.R")
source("genome_wide_r/gene_barcode.R")
source("genome_wide_r/volcano_plot_helpers.R")


GENOME_WIDE_DB = "Z:/Marco/Temp/Backup D/Temp/MBdata/Processed_Again/MBCrisprMBAgain_1_1.db" ##still change
CANDIDATE_DB = "Z:/Marco/CustomClonedLibrary/oPools seq files/MBCrisprMBSubscreeen_Full_1MHfSNVtop150.db"   # new, test if same output as originalis a

GENOME_WIDE_TYPES = c( "DELETION", "DELINS","HDR",
                       "INSERTION", "INSERTION_1bp","LARGE_DELETION",
                       "PQ_DELETION" ,"TINS","WT")

GENOME_WIDE_LABEL_MAP <- c(
  MB01 = "Target 1 - Replicate 1",
  MB02 = "Target 1 - Replicate 2",
  MB03 = "Target 2 - Replicate 1",
  MB04 = "Target 2 - Replicate 2",
  MB05 = "Target 3 - Replicate 1",
  MB06 = "Target 3 - Replicate 2",
  MB07 = "Target 1 - InhibitorX?",
  MB08 = "Target 1 - InhibitorY?"
)

genome_wide_conn <- function(){
  if(!file.exists(GENOME_WIDE_DB)){
    stop(paste("File does not exist:",GENOME_WIDE_DB))
  }
  con <- dbConnect(RSQLite::SQLite(), dbname = GENOME_WIDE_DB)
}

candidate_conn <- function(){
  if(!file.exists(CANDIDATE_DB)){
    stop(paste("File does not exist:",CANDIDATE_DB))
  }
  con <- dbConnect(RSQLite::SQLite(), dbname = CANDIDATE_DB)
}

##collect the genes/barcodes to keep the data static
GENOME_WIDE_GENE_INFO <- tbl(genome_wide_conn(), "countTable") %>% select(Gene, Barcode) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarise(nrBarcodes = n()) %>% 
  collect()

##collect the genes/barcodes to keep the data static
FOCUSED_GENE_INFO <- tbl(candidate_conn(), "countTable") %>% select(Gene, Barcode) %>%
  distinct() %>%
  group_by(Gene) %>%
  summarise(nrBarcodes = n()) %>% 
  collect()


