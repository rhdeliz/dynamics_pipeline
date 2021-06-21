#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]
new_image_ending = args[2]
results_table_name = args[3]

# # CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
# 

# parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"
# new_image_ending = "_intensity_ref.tif"
# results_table_name = "_intensity.csv.gz"

library(dplyr)
library(parallel)
library(tidyr)
library(data.table)
library(ff)
library(igraph)

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = read.csv(directories_list)
processing_path = directories_list$path[directories_list$contains == "processing"]
extraction_path = file.path(processing_path, "05_IntensityExtraction")
colocalized_path = file.path(processing_path, "06_Colocalization")
input_path = directories_list$path[directories_list$contains == "input"]
summary_path = file.path(input_path, "summary.csv")
# Get image list
file_list = read.csv(summary_path)
# Get images table
image_list = NULL
image_list$table = paste0(file_list$protein_relative_path, results_table_name)
image_list$table = file.path(extraction_path, image_list$table)
# Get list of files that need colocalization
colocalization_list <- as_tibble(image_list)
# Get colocalization list
colocalization_list <-
  colocalization_list %>% 
  mutate(
    cell = dirname(table),
    image = dirname(cell),
    cohort = dirname(image),
    cohort = basename(cohort),
    cohort = file.path(colocalized_path, cohort)
  ) %>% 
  filter(
    file.exists(table)
  ) %>% 
  group_by(
    cell
  ) %>% 
  mutate(
    n = n()
  )

# Create re-colocalized path if it doesn't exist
if(!file.exists(colocalized_path)){
  dir.create(colocalized_path)
}
# Get cohort list
Cohorts <- unique(colocalization_list$cohort)
# Create cohort path if it doesn't exist
for(Cohort in Cohorts){
  if(!file.exists(Cohort)){
    dir.create(Cohort)
  }
}

# Get list of no colocalization needed
NoColocalizationNeeded <-
  colocalization_list %>% 
  filter(
    n == 1
  )

# Move images that don't need colocalization
ImagePath <- NoColocalizationNeeded$image
for(Image in ImagePath){
  old_path = Image
  Cohort = basename(dirname(Image))
  Image = basename(Image)
  new_path = file.path(colocalized_path, Cohort, Image)
  file.move(old_path, new_path)
}

# Get list of files that need colocalization
ColocalizationNeeded <-
  colocalization_list %>% 
  filter(
    n > 1
  )
# Screen which frames have multiple proteins
ScreeningFx <- function(CellX){
  tryCatch({
    # Get tables
    TableList <-
      ColocalizationNeeded %>% 
      filter(
        cell == Cells[CellX]
      )
    Tables <- lapply(TableList$table, fread)
    Tables <- rbindlist(Tables, fill = T)
    # Get only frames with multiple puncta
    Tables <-
      Tables %>% 
      arrange(
        PROTEIN,
        TRACK_ID,
        FRAME
      ) %>% 
      ungroup() %>% 
      group_by(
        FRAME
      ) %>% 
      mutate(
        N = NROW(unique(PROTEIN)),
      ) %>% 
      filter(
        N > 1
      )
    return(Tables)
  }, error = function(e){print(paste("ERROR with MoveNoColocalizationNeededFx CellX =", CellX))})
}
Cells <- unique(ColocalizationNeeded$cell)
ScreenedTables <- mclapply(1:NROW(Cells), ScreeningFx)
ScreenedTables <- ScreenedTables[(which(sapply(ScreenedTables,is.list), arr.ind=TRUE))]
ScreenedTables <- rbindlist(ScreenedTables, fill=TRUE)

SplitCellTables <-
  ScreenedTables %>% 
  mutate(
    CELL_PATH = dirname(RELATIVE_PATH),
    CELL_PATH = file.path(extraction_path, CELL_PATH),
    IMAGE_PATH = paste0(RELATIVE_PATH, new_image_ending),
    IMAGE_PATH = file.path(extraction_path, IMAGE_PATH)
  ) %>% 
  group_split(
    CELL_PATH
  )

GetImageTablePairs <- function(CellTable){
  tryCatch({
    # Get combinations
    Proteins = unique(CellTable$PROTEIN)
    Proteins = expand.grid(Proteins, Proteins)
    names(Proteins) = c("IMAGE", "TABLE")
    Proteins <-
      Proteins %>%
      filter(
        IMAGE != TABLE
      ) %>% 
      mutate(
        COMPLEMENTARY_PROTEIN = IMAGE,
        SOURCE_TABLE = TABLE,
        
        IMAGE = paste0(IMAGE, new_image_ending),
        IMAGE = file.path(CellTable$CELL_PATH[1], IMAGE),
        
        TABLE = paste0(TABLE, results_table_name),
        TABLE = file.path(CellTable$CELL_PATH[1], TABLE)
      ) %>% 
      filter(
        file.exists(IMAGE),
        file.exists(TABLE)
      ) %>% 
      group_by(
        SOURCE_TABLE
      ) %>% 
      mutate(
        COLUMN_NAME = 1:n(),
      )
    return(Proteins)
  }, error = function(e){print(paste("ERROR with GetImageTablePairs"))})
}
PairsList <- mclapply(SplitCellTables, GetImageTablePairs)
PairsList <- PairsList[(which(sapply(PairsList,is.list), arr.ind=TRUE))]
PairsList <- rbindlist(PairsList, fill = TRUE)

GetIntensities <- function(PairX){
  tryCatch({
    # Get parameters
    image_path <- PairsList$IMAGE[PairX]
    table_path <- PairsList$TABLE[PairX]
    CellTable <- fread(table_path)
    puncta_radius = CellTable$PUNCTA_DIAMETER[1]/2
    
    #Import images
    img = ijtiff::read_tif(image_path, frames = unique(CellTable$FRAME), msg = FALSE)
    GetFrameTable <- function(FrameX){
      tryCatch({
        img_frame = which(unique(CellTable$FRAME)==FrameX)
        Z = img[,,,img_frame]
        X = NROW(Z)
        X = rep(1:X, NCOL(Z))
        
        Y = NCOL(Z)
        Y = rep(1:Y, each = NROW(Z))
        
        Z = as.vector(Z)
        frame_img = cbind(X, Y, Z)
        frame_img = as.data.frame(frame_img)
        names(frame_img) = c("x", "y", "z")
        frame_img$t = FrameX
        
        FrameCellTable <- CellTable[CellTable$FRAME == FrameX,]
        
        Table <- NULL
        Table[[1]] = frame_img
        Table[[2]] = as.data.frame(FrameCellTable)
        
        return(Table)
      }, error = function(e){print(paste("     ERROR with GetFrameTable FrameX =", FrameX))})
    }
    N_FRAMES = unique(CellTable$FRAME)
    N_FRAMES = sort(N_FRAMES)
    Tables <- mclapply(N_FRAMES, GetFrameTable)
    # Get frames
    FrameFx <- function(TableX){
      tryCatch({
        frame_img = TableX[[1]]
        FrameShortFullTable = TableX[[2]]
        # Get pixel intensity
        SubPixelLocalization <- function(SpotX){
          # Get coordinates
          x.coordinate = FrameShortFullTable$POSITION_X[SpotX]
          x.coordinate = as.numeric(x.coordinate)
          y.coordinate = FrameShortFullTable$POSITION_Y[SpotX]
          y.coordinate = as.numeric(y.coordinate)
          id = FrameShortFullTable$UNIVERSAL_SPOT_ID[SpotX]
          # Filter out spot
          TotalIntensity <- frame_img[frame_img$x >= x.coordinate - puncta_radius - .5 &
                                        frame_img$x <= x.coordinate + puncta_radius + .5 &
                                        frame_img$y >= y.coordinate - puncta_radius - .5 &
                                        frame_img$y <= y.coordinate + puncta_radius + .5,]
          
          TotalIntensity$x.low <- TotalIntensity$x - .5
          TotalIntensity$x.high <- TotalIntensity$x + .5
          TotalIntensity$x.coord.low <- x.coordinate - puncta_radius
          TotalIntensity$x.coord.high <- x.coordinate + puncta_radius
          TotalIntensity$x.high.capture <-
            data.table::fifelse(
              TotalIntensity$x.coord.high > TotalIntensity$x.high,
              TotalIntensity$x.high,
              TotalIntensity$x.coord.high
            )
          TotalIntensity$x.low.capture <-
            data.table::fifelse(
              x.coordinate - puncta_radius < TotalIntensity$x.low,
              TotalIntensity$x.low,
              TotalIntensity$x.coord.low
            )
          TotalIntensity$x.capture <- TotalIntensity$x.high.capture - TotalIntensity$x.low.capture
          
          TotalIntensity$y.low <- TotalIntensity$y - .5
          TotalIntensity$y.high <- TotalIntensity$y + .5
          TotalIntensity$y.coord.low <- y.coordinate - puncta_radius
          TotalIntensity$y.coord.high <- y.coordinate + puncta_radius
          TotalIntensity$y.high.capture <-
            data.table::fifelse(
              TotalIntensity$y.coord.high > TotalIntensity$y.high,
              TotalIntensity$y.high,
              TotalIntensity$y.coord.high
            )
          TotalIntensity$y.low.capture <-
            data.table::fifelse(
              TotalIntensity$y.coord.low < TotalIntensity$y.low,
              TotalIntensity$y.low,
              TotalIntensity$y.coord.low
            )
          
          TotalIntensity$y.capture <- TotalIntensity$y.high.capture - TotalIntensity$y.low.capture
          
          TotalIntensity$SPOT_AREA <- TotalIntensity$x.capture * TotalIntensity$y.capture
          TotalIntensity$TOTAL_INTENSITY <- TotalIntensity$SPOT_AREA*TotalIntensity$z
          
          TotalIntensity = TotalIntensity[TotalIntensity$x.capture > 0,]
          TotalIntensity = TotalIntensity[TotalIntensity$y.capture > 0,]
          
          Result <- NULL
          Result$UNIVERSAL_SPOT_ID = id
          Result$TOTAL_INTENSITY = sum(TotalIntensity$TOTAL_INTENSITY)
          Result$STANDARD_DEVIATION = sd(TotalIntensity$TOTAL_INTENSITY)
          
          return(Result)
        }
        Intensities <- lapply(1:NROW(FrameShortFullTable), SubPixelLocalization)
        Intensities <- data.table::rbindlist(Intensities)
        return(Intensities)
      }, error = function(e){print(paste("     ERROR with FrameFx"))})
    }
    Intensities <- mclapply(Tables, FrameFx)
    Intensities <- Intensities[(which(sapply(Intensities,is.list), arr.ind=TRUE))]
    Intensities <- data.table::rbindlist(Intensities)
    # Add complementary protein name
    Intensities$COMPLEMENTARY_PROTEIN = PairsList$COMPLEMENTARY_PROTEIN[PairX]
    COLUMN_NUMBER <- PairsList$COLUMN_NAME[PairX]
    names(Intensities) <-c(
      "UNIVERSAL_SPOT_ID",
      paste0("COMPLEMENTARY_TOTAL_INTENSITY_", COLUMN_NUMBER),
      paste0("COMPLEMENTARY_STANDARD_DEVIATION_", COLUMN_NUMBER),
      paste0("COMPLEMENTARY_PROTEIN_", COLUMN_NUMBER)
      )
    
    TableName <- paste(PairsList$SOURCE_TABLE[PairX],  PairsList$COMPLEMENTARY_PROTEIN[PairX], "colocalization_intensity.csv.gz", sep = "_")
    TableName <- file.path(dirname(PairsList$IMAGE[PairX]), TableName)
    fwrite(Intensities, TableName, row.names = F, na = "")
    
    return(TableName)
  }, error = function(e){print(paste("ERROR with GetIntensities. PairX =", PairX))})
}
lapply(1:NROW(PairsList), GetIntensities)
