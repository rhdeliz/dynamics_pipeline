#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]
new_image_ending = args[2]
results_table_name = args[3]
# 
# # CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)


# parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"
# new_image_ending = "_intensity_ref.tif"
# results_table_name = "_intensity.csv.gz"

library(dplyr)
library(parallel)
library(tidyr)
library(data.table)
library(ff)

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = read.csv(directories_list)
processing_path = directories_list$path[directories_list$contains == "processing"]
tracking_path = file.path(processing_path, "04_Track")
extraction_path = file.path(processing_path, "05_IntensityExtraction")
constants_path = file.path(parameters_path, "constants.csv") 
input_path = directories_list$path[directories_list$contains == "input"]
summary_path = file.path(input_path, "summary.csv")

# Get puncta diameter
puncta_diameter = read.csv(constants_path)
puncta_diameter = puncta_diameter$value[puncta_diameter$parameter == "puncta_diameter"]
puncta_radius = puncta_diameter * 0.5
# Get image parameters table
file_list = read.csv(summary_path)
# Get images table
image_list = NULL
image_list$image = paste0(file_list$protein_relative_path, new_image_ending)
image_list$image = file.path(tracking_path, image_list$image)
image_list$xml = paste0(file_list$protein_relative_path, ".xml")
image_list$xml = file.path(tracking_path, image_list$xml)
image_list$output = paste0(file_list$protein_relative_path, results_table_name)
image_list$output = file.path(tracking_path, image_list$output)
image_list <- as_tibble(image_list)
# Keep only files that exist
image_list <- image_list[file.exists(image_list$image) & file.exists(image_list$xml),]

# Functions
XMLtoTableFx <- function(FileX){
  tryCatch({
    # Get cell parameters
    cell_parameters = file_list[FileX,]
    cell_parameters <- as_tibble(cell_parameters)
    # Rename columns
    names(cell_parameters) <- toupper(names(cell_parameters))
    cell_parameters <-
      cell_parameters%>% 
      mutate(
        RELATIVE_PATH = PROTEIN_RELATIVE_PATH,
        PROTEIN = PROTEIN_NAME,
        CELL = dirname(PROTEIN_RELATIVE_PATH),
        IMAGE = dirname(CELL),
        CELL = basename(CELL),
        CELL = gsub("Cell_", "", CELL),
        IMAGE = basename(IMAGE),
        CELL_POSITION_X = POSITION_X,
        CELL_POSITION_Y = POSITION_Y,
        CELL_AREA = AREA
      ) %>% 
      select(-c(
        PROTEIN_RELATIVE_PATH,
        PROTEIN_NAME,
        POSITION_X,
        POSITION_Y,
        AREA
      ))
    # Get cell data
    image_path = image_list$image[FileX]
    xml_path = image_list$xml[FileX]
    save_path = image_list$output[FileX]
    
    # Import XML
    xml_data <- XML::xmlParse(xml_path)
    xml_data <- XML::xmlToList(xml_data)
    
    # Get spot data
    Spots <- xml_data$Model$AllSpots
    SpotsFx <- function(SpotX) {
      tryCatch({
        Table <- Spots[SpotX]$SpotsInFrame
        Table <- Table[1:(NROW(Table)-1)]
        Table <- as.data.frame(Table)
        Table <- t(Table)
        Table <- as.data.frame(Table)
        
        if(NCOL(Table) == 1) {
          Table = NULL
        }
        return(Table)
      }, error = function(e){print(paste("ERROR with SpotsFx. SpotX =", SpotX))})
    }
    SpotsTable <- mclapply(1:NROW(Spots), SpotsFx, mc.cores = detectCores(logical = F))
    SpotsTable <- SpotsTable[(which(sapply(SpotsTable,is.list), arr.ind=TRUE))]
    SpotsTable <- data.table::rbindlist(SpotsTable)
    
    # Get tracks data
    Tracks <- xml_data$Model$AllTracks
    TracksFx <- function(TrackX) {
      tryCatch({
        Table <- Tracks[TrackX]$Track
        Table <- Table[1:(NROW(Table)-1)]
        Table <- as.data.frame(Table)
        Table <- t(Table)
        Table <- as.data.frame(Table)
        
        Table$TRACK_ID = TrackX - 1
        
        return(Table)
      }, error = function(e){print(paste("ERROR with TracksFx TrackX =", TrackX))})
    }
    TracksTable <- mclapply(1:NROW(Tracks), TracksFx, mc.cores = detectCores(logical = F))
    TracksTable <- TracksTable[(which(sapply(TracksTable,is.list), arr.ind=TRUE))]
    TracksTable <- data.table::rbindlist(TracksTable)
    
    SPOTS1 <-
      TracksTable %>%
      mutate(
        ID = SPOT_SOURCE_ID
      ) %>%
      select(
        ID,
        TRACK_ID
      )
    
    SPOTS2 <-
      TracksTable %>%
      mutate(
        ID = SPOT_TARGET_ID
      ) %>%
      select(
        ID,
        TRACK_ID
      )
    
    SPOTS <- bind_rows(SPOTS1, SPOTS2) %>% distinct()
    SpotsTable <- merge(SpotsTable, SPOTS, by = "ID")
    remove(SPOTS, SPOTS1, SPOTS2)
    # Add one to account for running average
    SpotsTable$FRAME <- as.numeric(SpotsTable$FRAME) + 1
    SpotsTable$POSITION_X <- as.numeric(SpotsTable$POSITION_X)
    SpotsTable$POSITION_Y <- as.numeric(SpotsTable$POSITION_Y)
    # Generate table to find missing spots
    MissingSpotsTables <-
      SpotsTable %>% 
      group_by(
        TRACK_ID
      ) %>% 
      mutate(
        RANGE = max(FRAME) - min(FRAME),
        N = n() - 1
      ) %>% 
      filter(
        RANGE != N
      ) %>% 
      ungroup() %>% 
      select(
        TRACK_ID,
        FRAME,
        POSITION_X,
        POSITION_Y
      ) %>% 
      arrange(
        TRACK_ID,
        FRAME,
      ) %>% 
      group_split(
        TRACK_ID
      )
    
    # Add missing puncta only if it's not a calibration
    # Dark-phases are normal in fluorophores
    if(cell_parameters$COHORT != "Calibrations"){
      MissingSpotsFx <- function(TableX){
        tryCatch({
          Table <- TableX
          Table$MISSING = FALSE
          # Get actual frames
          ActualFrames <- unique(Table$FRAME)
          # Get predicted frames
          PredictedFrames <- min(Table$FRAME):max(Table$FRAME)
          # Find missing frames
          MissingFrames <- NULL
          MissingFrames$FRAME <- PredictedFrames[!PredictedFrames %in% ActualFrames]
          MissingFrames$TRACK_ID <- rep(Table$TRACK_ID[1], NROW(MissingFrames$FRAME))
          MissingFrames$MISSING <- rep(TRUE, NROW(MissingFrames$FRAME))
          # Create new table
          MissingFrames <- bind_rows(Table, MissingFrames)
          
          # Fill coordinates
          MissingFrames <-
            MissingFrames %>% 
            arrange(
              FRAME
            ) %>% 
            mutate(
              PREVIOUS_POSITION_X = POSITION_X,
              NEXT_POSITION_X = POSITION_X,
              PREVIOUS_POSITION_Y = POSITION_Y,
              NEXT_POSITION_Y = POSITION_Y
            ) %>% 
            fill(
              PREVIOUS_POSITION_X,
              PREVIOUS_POSITION_Y,
              .direction = "down"
            ) %>% 
            fill(
              NEXT_POSITION_X,
              NEXT_POSITION_Y,
              .direction = "up"
            ) %>% 
            filter(
              MISSING == TRUE
            ) %>% 
            group_by(
              PREVIOUS_POSITION_X,
              PREVIOUS_POSITION_Y,
              NEXT_POSITION_X,
              NEXT_POSITION_Y
            ) %>% 
            mutate(
              POSITION = 1:n() / (n() + 1),
              POSITION_X = (NEXT_POSITION_X - PREVIOUS_POSITION_X)*POSITION + PREVIOUS_POSITION_X,
              POSITION_Y = (NEXT_POSITION_Y - PREVIOUS_POSITION_Y)*POSITION + PREVIOUS_POSITION_Y
            ) %>% 
            ungroup() %>% 
            select(
              TRACK_ID,
              FRAME,
              POSITION_X,
              POSITION_Y
            )
          
        }, error = function(e){print(paste("     ERROR with MissingSpotsFx"))})
      }
      MissingSpotsTables <- mclapply(MissingSpotsTables, MissingSpotsFx)
      MissingSpotsTables <- rbindlist(MissingSpotsTables)
      # Add missing coordinates
      SpotsTable <- bind_rows(SpotsTable, MissingSpotsTables) 
    }
    
    # Combine cell parameters with spot data
    select_cell_parameters <-
      cell_parameters %>% 
      select(
        IMAGE, CELL, PROTEIN
      )
    # Add parameters to spots table
    ShortFullTable <- bind_cols(select_cell_parameters, SpotsTable)
    ShortFullTable <- 
      ShortFullTable %>% 
      mutate(
        UNIVERSAL_TRACK_ID = paste(IMAGE, CELL, PROTEIN, TRACK_ID, sep = "..."),
        UNIVERSAL_SPOT_ID = paste(UNIVERSAL_TRACK_ID, FRAME, sep = "..."),
        TRACKING_TOTAL_INTENSITY = TOTAL_INTENSITY,
        TRACKING_STANDARD_DEVIATION = STANDARD_DEVIATION
      ) %>% 
      select(-c(
        # Conflicting with merger
        IMAGE, CELL, PROTEIN, TOTAL_INTENSITY, STANDARD_DEVIATION,
        # Conflicting with analysis
        MAX_INTENSITY, MEDIAN_INTENSITY, VISIBILITY, MEAN_INTENSITY, RADIUS, MANUAL_COLOR, MIN_INTENSITY, POSITION_Z, POSITION_T
      ))
    #Import images
    img = ijtiff::read_tif(image_path, frames = unique(ShortFullTable$FRAME), msg = FALSE)
    GetFrameTable <- function(FrameX){
      tryCatch({
        img_frame = which(unique(ShortFullTable$FRAME)==FrameX)
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
        
        FrameShortFullTable <- ShortFullTable[ShortFullTable$FRAME == FrameX,]
        
        Table <- NULL
        Table[[1]] = frame_img
        Table[[2]] = as.data.frame(FrameShortFullTable)
        
        return(Table)
      }, error = function(e){print(paste("     ERROR with GetFrameTable FrameX =", FrameX))})
    }
    N_FRAMES = unique(ShortFullTable$FRAME)
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
          Result$SPOT_AREA = sum(TotalIntensity$SPOT_AREA)
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
    
    # Merge xml and image data
    CompleteIntensities <- merge(ShortFullTable, Intensities, by = "UNIVERSAL_SPOT_ID")
    CompleteIntensities <- bind_cols(CompleteIntensities, cell_parameters)
    CompleteIntensities <-
      CompleteIntensities %>%
      mutate(
        POSITION_X = as.numeric(POSITION_X),
        ABSOLUTE_POSITION_X = CELL_POSITION_X + POSITION_X,
        ABSOLUTE_POSITION_X = ABSOLUTE_POSITION_X*CALIBRATION_UM,
        
        POSITION_Y = as.numeric(POSITION_Y),
        ABSOLUTE_POSITION_Y = CELL_POSITION_Y + POSITION_Y,
        ABSOLUTE_POSITION_Y = ABSOLUTE_POSITION_Y*CALIBRATION_UM,
      ) %>% 
      arrange(
        TRACK_ID,
        FRAME
      )
      
    # Write table
    CompleteIntensities$ID = NULL
    CompleteIntensities$name = NULL
 
    data.table::fwrite(CompleteIntensities, save_path, row.names = F, na = "")

    if(FileX %% (detectCores(logical = FALSE)/2) == 0){
      Progress = FileX/NROW(image_list)
      Progress = Progress*100
      Progress = round(Progress)
      Progress = paste0("     ", Progress, "% complete")
      print(Progress)
    }
    return(save_path)
  }, error = function(e){print(paste("ERROR with XMLtoTableFx FileX =", FileX))})
}
ResultsPath <- lapply(1:NROW(image_list), XMLtoTableFx)
# Make vector of paths
ResultsPath <- as.vector(ResultsPath)
ResultsPath <- unlist(ResultsPath)
# Remove error rows
ResultsPath <- ResultsPath[!grepl("ERROR with XMLtoTableFx FileX =", ResultsPath)]
# Get images path
ImagePath <- dirname(dirname(ResultsPath))
ImagePath <- unique(ImagePath)
# Get cohort names
Cohorts <- basename(dirname(ImagePath))
Cohorts <- file.path(extraction_path, Cohorts)
Cohorts <- unique(Cohorts)
# Create re-extraction path if it doesn't exist
if(!file.exists(extraction_path)){
  dir.create(extraction_path)
}
# Create cohort path if it doesn't exist
for(Cohort in Cohorts){
  if(!file.exists(Cohort)){
    dir.create(Cohort)
  }
}

# Move Image
for(Image in ImagePath){
  old_path = Image
  Cohort = basename(dirname(Image))
  Image = basename(Image)
  new_path = file.path(extraction_path, Cohort, Image)
  file.move(old_path, new_path)
}

# Delete cohort if empty
for(Cohort in Cohorts){
  if(NROW(list.files(Cohorts)) == 0){
    unlink(Cohort)
  }
}
