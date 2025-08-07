library(dplyr)
library(dbplyr)
library(RSQLite)
library(ggplot2)
library(shiny)

GENOME_WIDE_DB = "Z:/Marco/Temp/Backup D/Temp/MBdata/Processed_Again/MBCrisprMBAgain_1_1.db" ##still change
CANDIDATE_DB = "Z:/Marco/CustomClonedLibrary/oPools seq files/MBCrisprMBSubscreeen_Full_1MHfSNVtop150.db"   # new, test if same output as originalis a

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

