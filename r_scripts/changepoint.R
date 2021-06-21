#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]
results_table_name = args[2]
# 
# # CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"
results_table_name = "_intensity.csv.gz"

library(dplyr)
library(stringr)
library(parallel)
library(data.table)
library(tidyr)
library(ff)

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = read.csv(directories_list)
processing_path = directories_list$path[directories_list$contains == "processing"]
colocalized_path = file.path(processing_path, "06_Colocalization")
changepoint_path = file.path(processing_path, "07_Changepoint")
input_path = directories_list$path[directories_list$contains == "input"]
summary_path = file.path(input_path, "summary.csv")
# Get image list
file_list = read.csv(summary_path)
# Prep up
file_list <-
  file_list %>% 
  mutate(
    path = paste0(protein_relative_path, results_table_name),
    path = file.path(colocalized_path, path)
  ) %>% 
  filter(
    file.exists(path),
    cohort != "Calibrations"
  )
file_list <- file_list$path

# SETUP
# Load libraries
if("pacman" %in% rownames(installed.packages()) == FALSE)
{install.packages("pacman")}

# Library list
pacman::p_load(dplyr, tidyr, data.table, changepoint, parallel, R.utils) 

# Data frame must have x, y and id columns
ChangepointFx <- function(df){
  tryCatch({
    # Expand time list
    x_sequence <- min(df$x):max(df$x)
    # List of times to remove after running changepoint
    remove_x <- setdiff(x_sequence, df$x)
    remove_x <- !x_sequence %in%remove_x
    
    # Adjust df table
    df <- 
      df %>% 
      complete(
        # Expand times
        x = seq(min(x), max(x))
      ) %>% 
      fill(
        # Add missing y with previous value
        y,
        .direction = "down"
      ) %>% 
      mutate(
        # Log transform y
        y = log(y + 1)
      ) %>% 
      distinct()
    
    # Run changepoint
    m.pelt <- cpt.mean(df$y, method = "PELT", penalty = "Manual", pen.value = 0.3, minseglen = 3)
    # Get y
    y <- param.est(m.pelt)$mean
    # Transform y back to linear
    y <- exp(y)-1
    # Expand y to match time
    len <- seg.len(m.pelt)
    y <- rep(y, times = len)
    # Get ID
    id <- df$id[1]
    # Get times
    x <- df$x
    # Make new table
    new_df <- suppressMessages(bind_cols(x, y, id))
    names(new_df) <- c("x", "y", "id")
    # Keep only original times
    new_df <- new_df[remove_x, ]
    return(new_df)
  }, error = function(e) {print("Error with ChangepointFx")})}

# Run tables
EvaluateTable <- function(TableX){
  tryCatch({
  # Get paths
  TablePath <- file_list[TableX]
  ExportPath <- dirname(TablePath)
  ExportName <- basename(TablePath)
  ExportName <- gsub(results_table_name, "_changepoint.csv.gz", ExportName)
  ExportPath <- file.path(ExportPath, ExportName)
  
  # Get intensity table
  IntensityTable <- fread(TablePath)
  
  # Prep table up
  FilteredTable <-
    IntensityTable %>% 
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    mutate(
      LIFETIME = max(FRAME) - min(FRAME)
    ) %>% 
    filter(
      LIFETIME >= 6
    ) %>%
    group_by(
      UNIVERSAL_TRACK_ID
    ) %>% 
    arrange(
      LIFETIME
    ) %>% 
    mutate(
      x = FRAME,
      y = TOTAL_INTENSITY,
      id = UNIVERSAL_TRACK_ID
    ) %>% 
    ungroup() %>% 
    select(
      x,
      y,
      id
    ) %>% 
    group_split(
      id
    )
  if(NROW(FilteredTable) > 0){
    # Run changepoint
    CPResults <- mclapply(FilteredTable, ChangepointFx, mc.cores = detectCores())
    CPResults <- CPResults[(which(sapply(CPResults,is.list), arr.ind=TRUE))]
    CPResults <- rbindlist(CPResults)
    # Get changepoint parameters
    CPResults$id <- paste(CPResults$id, CPResults$x, sep = "...")
    CPResults$x <- NULL
    names(CPResults) <- c("CHANGEPOINT_TOTAL_INTENSITY", "UNIVERSAL_SPOT_ID")
    
    fwrite(CPResults, ExportPath, row.names = F, na = "")
  }
  return(ExportPath)
  }, error = function(e) {print(paste("Error with EvaluateTable. TableX =", TableX))})}
lapply(1:NROW(file_list), EvaluateTable)

# Get image list
file_list = read.csv(summary_path)
# Prep up
file_list <-
  file_list %>% 
  mutate(
    path = dirname(protein_relative_path),
    path = dirname(path),
    old_path = file.path(colocalized_path, path),
    new_path = file.path(changepoint_path, path),
    cohort = basename(dirname(path)),
    cohort = file.path(changepoint_path, cohort)
  ) %>% 
  distinct() %>% 
  select(
    cohort,
    old_path,
    new_path
  )

# Create re-colocalized path if it doesn't exist
if(!file.exists(changepoint_path)){
  dir.create(changepoint_path)
}
# Create cohort folders
lapply(unique(file_list$cohort), dir.create)

# Move images to changepoint folder
for(ImageX in 1:NROW(file_list)){
  old_path = file_list$old_path[ImageX]
  new_path = file_list$new_path[ImageX]
  file.move(old_path, new_path)
}
