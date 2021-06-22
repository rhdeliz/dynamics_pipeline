#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]

# CLEAN ENVIRONMENT----
# remove(list = ls())
# gc(reset = TRUE)
# pacman::p_unload(pacman::p_loaded(), character.only = TRUE)
# 
# parameters_path = "/Users/u_deliz/Desktop/NewPipeline/Input/parameter_tables/"

# Ending of files
intensity_ending = "_intensity.csv.gz"
colocalization_intensity_ending = "_colocalization_intensity.csv.gz"
colocalization_coordinates_filename = "colocalization_coordinates.csv.gz"
changepoint_ending = "_changepoint.csv.gz"

library(dplyr)
library(stringr)
library(parallel)
library(data.table)
library(tidyr)

# Get directories
directories_list = file.path(parameters_path, "directories.csv")
directories_list = read.csv(directories_list)
processing_path = directories_list$path[directories_list$contains == "processing"]
changepoint_path = file.path(processing_path, "07_Changepoint")
input_path = directories_list$path[directories_list$contains == "input"]
output_path = directories_list$path[directories_list$contains == "output"]
summary_path = file.path(input_path, "summary.csv")
# Get image list
file_list = read.csv(summary_path)
# Prep up
file_list <-
  file_list %>% 
  mutate(
    image = dirname(dirname(protein_relative_path)),
    image = basename(image),
    date = substr(image, 0, 8),
    date =  as.Date(date, format = '%Y%m%d'),
    exposure = ifelse(word(exposure, -1) == "ms", as.numeric(word(exposure, 1))/1000, ifelse(word(exposure, -1) == "s", as.numeric(word(exposure, 1)), exposure)),
    direction = direction*pi/180,
    angle = angle*pi/180
  )

# Separate calibrations and the rest
calibration_list = file_list[file_list$cohort=="Calibrations",]
image_list = file_list[file_list$cohort!="Calibrations",]

# Pairing
GetCalibrationImages <- function(ImageX){
  tryCatch({
    # Get calibration date
    IMAGE_DATE = image_list$date[ImageX]
    CALIBRATION_DATES = calibration_list$date
    DATE_DIFFERENCE = abs(CALIBRATION_DATES - IMAGE_DATE)
    CALIBRATION_DATE = CALIBRATION_DATES[which.min(DATE_DIFFERENCE)]
    
    # Calibration File
    CalibrationFile <-
      calibration_list %>% 
      filter(
        channel == image_list$channel[ImageX],
        date == CALIBRATION_DATE
      )
    
    if(NROW(CalibrationFile)>0){
      # Get calibration
      # Nearest laser
      PowerFiltered = CalibrationFile[which.min(abs(CalibrationFile$power - image_list$power[ImageX])),]
      # Nearest exposure
      ExposureFiltered = PowerFiltered[which.min(abs(PowerFiltered$exposure - image_list$exposure[ImageX])),]
      # Nearest direction
      DirectionFiltered = PowerFiltered[which.min(abs(tan(ExposureFiltered$direction-image_list$direction[ImageX]))),]
      AngleFiltered = DirectionFiltered[which.min(abs(tan(DirectionFiltered$angle-image_list$angle[ImageX]))),]
      CalibrationImage = AngleFiltered$image[NROW(AngleFiltered)]
      CalibrationProtein =  AngleFiltered$protein_name[1]
      
    } else{
      # Get calibration for other dates
      CalibrationFile <-
        calibration_list %>% 
        filter(
          channel == image_list$channel[ImageX]
        )
      
      if(NROW(CalibrationFile)>0){
        # Nearest laser
        PowerFiltered = CalibrationFile[which.min(abs(CalibrationFile$power - image_list$power[ImageX])),]
        # Nearest exposure
        ExposureFiltered = PowerFiltered[which.min(abs(PowerFiltered$exposure - image_list$exposure[ImageX])),]
        # Nearest direction
        DirectionFiltered = PowerFiltered[which.min(abs(tan(ExposureFiltered$direction-image_list$direction[ImageX]))),]
        AngleFiltered = DirectionFiltered[which.min(abs(tan(DirectionFiltered$angle-image_list$angle[ImageX]))),]
        # Nearest date
        DateFiltered = AngleFiltered[which.min(abs(AngleFiltered$date - image_list$date[ImageX])),]
        # Last picture taken
        CalibrationImage = DateFiltered$image[NROW(DateFiltered)]
        CalibrationProtein =  DateFiltered$protein_name[1]
      } else{
        CalibrationImage = NA
        CalibrationProtein = NA
      }
    }
    # Create export table
    ExportTable <- image_list[ImageX,]
    ExportTable$CALIBRATION_IMAGE = CalibrationImage
    ExportTable$FLUOROPHORE = CalibrationProtein
    
    return(ExportTable)
  }, error = function(e){print(paste("ERROR with MoveNoColocalizationNeededFx ImageX =", ImageX))})
}
PairedList <- mclapply(1:NROW(image_list), GetCalibrationImages)
PairedList <- rbindlist(PairedList, fill = TRUE)
PairedList <- PairedList %>% distinct()

CalibrationImages <-
  PairedList %>% 
  select(
    CALIBRATION_IMAGE,
    FLUOROPHORE
  ) %>%
  distinct() %>% 
  drop_na()

GetCalibrationIntensity <- function(ImageX){
  tryCatch({
    CALIBRATION_IMAGE = CalibrationImages$CALIBRATION_IMAGE[ImageX]
    FLUOROPHORE = CalibrationImages$FLUOROPHORE[ImageX]
    CALIBRATION_IMAGE = file.path(changepoint_path, "Calibrations", CALIBRATION_IMAGE, "Cell_1", paste0(FLUOROPHORE, "_intensity.csv.gz"))
    
    IntensitiesTable <- fread(CALIBRATION_IMAGE)
    
    IntensitiesTable <-
      IntensitiesTable %>% 
      mutate(
        CALIBRATION_IMAGE = IMAGE
      ) %>% 
      filter(
        SPOT_AREA == PUNCTA_DIAMETER^2
      ) %>% 
      group_by(
        CALIBRATION_IMAGE
      ) %>% 
      summarize(
        CALIBRATION_STANDARD_DEVIATION = sd(TOTAL_INTENSITY),
        CALIBRATION_TOTAL_INTENSITY = mean(TOTAL_INTENSITY)
      )
    return(IntensitiesTable)
  }, error = function(e){print(paste("ERROR with GetCalibrationIntensity ImageX =", ImageX))})
}
CalibrationIntensities <- mclapply(1:NROW(CalibrationImages), GetCalibrationIntensity)
CalibrationIntensities <- CalibrationIntensities[(which(sapply(CalibrationIntensities,is.list), arr.ind=TRUE))]
CalibrationIntensities <- rbindlist(CalibrationIntensities, fill = TRUE)
# Pair image and calibration
PairedCalibrations <- merge(PairedList, CalibrationIntensities, by = "CALIBRATION_IMAGE", all = TRUE)
names(PairedCalibrations) <- toupper(names(PairedCalibrations))
PairedCalibrations <-
  PairedCalibrations %>%
  mutate(
    CELL_PATH = dirname(PROTEIN_RELATIVE_PATH),
    CALIBRATION_TOTAL_INTENSITY = ifelse(is.na(CALIBRATION_TOTAL_INTENSITY), 1, CALIBRATION_TOTAL_INTENSITY),
    FLUOROPHORE = ifelse(is.na(FLUOROPHORE), CHANNEL, FLUOROPHORE)
  ) %>%
  distinct()

CombineCellTables <- function(CellX){
  tryCatch({
    # Get tables
    CellPairedCalibrations <-
      PairedCalibrations %>% filter(
        CELL_PATH ==  Cells[CellX]
      ) %>% 
      mutate(
        PROTEIN = PROTEIN_NAME 
      ) %>% 
      select(
        PROTEIN,
        CALIBRATION_IMAGE,
        CALIBRATION_TOTAL_INTENSITY,
        CALIBRATION_STANDARD_DEVIATION
      )
    
    # Get cell path
    CellPath <- file.path(changepoint_path, Cells[CellX])
    
    # Get protein combinations
    Proteins = unique(CellPairedCalibrations$PROTEIN)
    Proteins = expand.grid(Proteins, Proteins)
    names(Proteins) = c("IMAGE", "TABLE")
    Proteins <-
      Proteins %>%
      filter(
        IMAGE != TABLE
      ) %>% 
      mutate(
        TABLE = paste0(IMAGE, "_", TABLE, colocalization_intensity_ending)
      )
    
    # Get tables
    IntensityTables <- paste0(CellPairedCalibrations$PROTEIN, intensity_ending)
    IntensityTables <- file.path(CellPath, IntensityTables)
    IntensityTables <- IntensityTables[file.exists(IntensityTables)]
    if(NROW(IntensityTables)>0){
      
      IntensityTables <- lapply(IntensityTables, fread)
      IntensityTables <- rbindlist(IntensityTables, fill = TRUE)
      # Pair calibration
      IntensityTables <- merge(IntensityTables, CellPairedCalibrations, by = "PROTEIN", all = TRUE)
      
      # Add colocalization intensities (image and coordinates flipped)
      ColocalizationIntensityTables <- file.path(CellPath, Proteins$TABLE)
      ColocalizationIntensityTables <- ColocalizationIntensityTables[file.exists(ColocalizationIntensityTables)]
      ColocalizationIntensityTables <- if(NROW(ColocalizationIntensityTables)>0){
        ColocalizationIntensityTables <- lapply(ColocalizationIntensityTables, fread)
        ColocalizationIntensityTables <- rbindlist(ColocalizationIntensityTables, fill = TRUE)
        
        # Normalize colocalization intensity
        ComplementaryProteins <- names(ColocalizationIntensityTables)
        ComplementaryProteinsIndex <- substr(ComplementaryProteins, 0, nchar("COMPLEMENTARY_PROTEIN_"))
        ComplementaryProteinsIndex <- which(ComplementaryProteinsIndex == "COMPLEMENTARY_PROTEIN_")
        ComplementaryProteins <- ComplementaryProteins[ComplementaryProteinsIndex]
        ComplementaryProteins <- grep("COMPLEMENTARY_PROTEIN_", ComplementaryProteins)
        
        ColocalizationIntensityNormalization <- function(ProteinX){
          tryCatch({
            # Get protein name
            ProteinColumnName <- paste0("COMPLEMENTARY_PROTEIN_", ProteinX)
            ProteinName = which(names(ColocalizationIntensityTables) == ProteinColumnName)
            ProteinName <- ColocalizationIntensityTables[, ..ProteinName]
            # Get calibration intensity
            ProteinCalibration <- merge(ProteinName, CellPairedCalibrations, by.x = ProteinColumnName, by.y = "PROTEIN")
            
            # Get intensity
            IntensityColumnName <- paste0("COMPLEMENTARY_TOTAL_INTENSITY_", ProteinX)
            ProteinIntensity = which(names(ColocalizationIntensityTables) == IntensityColumnName)
            ProteinIntensity <- ColocalizationIntensityTables[, ..ProteinIntensity]
            ProteinIntensity <- ProteinIntensity/ProteinCalibration$CALIBRATION_TOTAL_INTENSITY
            names(ProteinIntensity) <- paste0("COMPLEMENTARY_NORMALIZED_INTENSITY_", ProteinX)
            return(ProteinIntensity)
          }, error = function(e){print(paste("ERROR with ColocalizationIntensityNormalization ProteinX =", ProteinX))})
        }
        ColocalizationNormalizedIntensities <- lapply(ComplementaryProteins, ColocalizationIntensityNormalization)
        ColocalizationNormalizedIntensities <- ColocalizationNormalizedIntensities[(which(sapply(ColocalizationNormalizedIntensities,is.list), arr.ind=TRUE))]
        ColocalizationNormalizedIntensities <- rbindlist(ColocalizationNormalizedIntensities, fill = TRUE)
        ColocalizationIntensityTables <- cbind(ColocalizationIntensityTables, ColocalizationNormalizedIntensities)
          
        IntensityTables <- merge(IntensityTables, ColocalizationIntensityTables, by = "UNIVERSAL_SPOT_ID", all = TRUE)
      }
      
      # Add colocalization data (coordinates paired)
      ColocalizationCoordinates <- file.path(CellPath, colocalization_coordinates_filename)
      ColocalizationCoordinates <- ColocalizationCoordinates[file.exists(ColocalizationCoordinates)]
      ColocalizationCoordinates <- if(NROW(ColocalizationCoordinates)>0){
        ColocalizationCoordinates <- fread(ColocalizationCoordinates)
        ColocalizationCoordinates$RELATIVE_PATH <- NULL
        IntensityTables <- merge(IntensityTables, ColocalizationCoordinates, by = "UNIVERSAL_SPOT_ID", all = TRUE)
      }
      
      # Add changepoint analysis
      ChangepointTables <- paste0(CellPairedCalibrations$PROTEIN, changepoint_ending)
      ChangepointTables <- file.path(CellPath, ChangepointTables)
      ChangepointTables <- ChangepointTables[file.exists(ChangepointTables)]
      ChangepointTables <- if(NROW(ChangepointTables)>0){
        ChangepointTables <- lapply(ChangepointTables, fread)
        ChangepointTables <- ChangepointTables[(which(sapply(ChangepointTables,is.list), arr.ind=TRUE))]
        ChangepointTables <- rbindlist(ChangepointTables, fill = TRUE)
        IntensityTables <- merge(IntensityTables, ChangepointTables, by = "UNIVERSAL_SPOT_ID", all = TRUE)
      }
      
      # Determine nearest neighbor to spot of the same protein
      NearestNeighborTable <-
        IntensityTables %>% 
        group_by(
          PROTEIN,
          FRAME
        ) %>% 
        mutate(
          N = n()
        ) %>% 
        filter(
          N > 1
        ) %>% 
        select(
          PROTEIN,
          FRAME,
          PUNCTA_DIAMETER,
          UNIVERSAL_SPOT_ID,
          POSITION_X,
          POSITION_Y
        ) %>% 
        ungroup() %>% 
        group_split(
          PROTEIN,
          FRAME
        )
      # Nearest neighbor
      RADIUS <- max(IntensityTables$PUNCTA_DIAMETER, na.rm = T) * 3
      NearestNeighborSearch <- function(ProteinFrameX){
        # Get nearest spot
        CoordiantesTable <- ProteinFrameX %>% select(POSITION_X, POSITION_Y)
        DistanceTable <- RANN::nn2(CoordiantesTable, k = 2)
        DistanceTable <- DistanceTable$nn.dists[,2]
        # Get number of spots nearby
        ClusterTable <- RANN::nn2(CoordiantesTable, searchtype = c("radius"), radius = RADIUS, k = NROW(ProteinFrameX))
        ClusterTable <- ClusterTable$nn.dists
        ClusterTable <- ClusterTable > 0 & ClusterTable <  1e+153 
        ClusterTable <- rowSums(ClusterTable)
        # Put results together
        DistanceResults <- ProteinFrameX %>% select(UNIVERSAL_SPOT_ID)
        DistanceResults$NEAREST_SPOT <- DistanceTable
        DistanceResults$SPOTS_WITHIN_RADIUS = ClusterTable
        DistanceResults$SPOT_RADIUS_LIMIT = RADIUS
        return(DistanceResults)
      }
      NeighborResults <- lapply(NearestNeighborTable, NearestNeighborSearch)
      NeighborResults <- NeighborResults[(which(sapply(NeighborResults,is.list), arr.ind=TRUE))]
      NeighborResults <- rbindlist(NeighborResults, fill = TRUE)
      IntensityTables <- merge(IntensityTables, NeighborResults, by = "UNIVERSAL_SPOT_ID", all = TRUE)
      # Add missing values
      MissingIndex <- which(is.na(IntensityTables$NEAREST_SPOT))
      IntensityTables$NEAREST_SPOT[MissingIndex] <- Inf
      MissingIndex <- which(is.na(IntensityTables$SPOTS_WITHIN_RADIUS))
      IntensityTables$SPOTS_WITHIN_RADIUS[MissingIndex] <- 0
      MissingIndex <- which(is.na(IntensityTables$SPOT_RADIUS_LIMIT))
      IntensityTables$SPOT_RADIUS_LIMIT[MissingIndex] <- RADIUS
      
      # Get landing frame
      LANDING_FRAME = min(IntensityTables$FRAME)
      
      # Calculate basic parameters and save
      IntensityTables <-
        IntensityTables %>% 
        filter(
          # Filter small spots (probably at edges)
          SPOT_AREA >= PUNCTA_DIAMETER^2
        ) %>% 
        group_by(
          UNIVERSAL_TRACK_ID
        ) %>% 
        mutate(
          # Get min/max spot xy-coordinates
          MIN_X = min(ABSOLUTE_POSITION_X/CALIBRATION_UM),
          MAX_X = max(ABSOLUTE_POSITION_X/CALIBRATION_UM),
          MIN_Y = min(ABSOLUTE_POSITION_Y/CALIBRATION_UM),
          MAX_Y = max(ABSOLUTE_POSITION_Y/CALIBRATION_UM)
        ) %>% 
        filter(
          # Remove spots near edges
          MIN_X > PUNCTA_DIAMETER*2,
          MIN_Y > PUNCTA_DIAMETER*2,
          MAX_X < WIDTH - PUNCTA_DIAMETER*2,
          MAX_Y < HEIGHT - PUNCTA_DIAMETER*2
        ) %>% 
        select(-c(
          MIN_X,
          MAX_X,
          MIN_Y,
          MAX_Y
        )) %>% 
        group_by(
          PROTEIN
        ) %>% 
        mutate(
          # Normalize intensities
          NORMALIZED_INTENSITY = TOTAL_INTENSITY/CALIBRATION_TOTAL_INTENSITY,
          NORMALIZED_CHANGEPOINT_INTENSITY = CHANGEPOINT_TOTAL_INTENSITY/CALIBRATION_TOTAL_INTENSITY
        ) %>% 
        group_by(
          UNIVERSAL_TRACK_ID
        ) %>% 
        mutate(
          # Ligand density every 0.5 steps using log base 10
          LIGAND_DENSITY_CAT =
            # Round to nearest half log
            round(
              # Classifies ligands based on log using base 3.162 (10^.5)
              log(
                LIGAND_DENSITY,
                base = (10^.5)),
              digits = 0
            ) * 0.5,
          # Convert back to linear
          LIGAND_DENSITY_CAT = signif(10^LIGAND_DENSITY_CAT, 2),
          # Get frame number
          FRAMES_ADJUSTED = FRAME - min(FRAME),
          # Get puncta lifetime
          LIFETIME = max(FRAMES_ADJUSTED)/FRAME_RATE,
          # Get actual time
          TIME = FRAMES_ADJUSTED/FRAME_RATE,
          # Time since first spot
          TIME_SINCE_LANDING = (FRAME - LANDING_FRAME)/FRAME_RATE,
          # Max total intensity of track
          MAX_NORMALIZED_INTENSITY = max(NORMALIZED_INTENSITY, na.rm = TRUE),
          # To be used later for overall delta and for categorizing de-novo and disassembling
          STARTING_NORMALIZED_INTENSITY = NORMALIZED_INTENSITY[1],
          # Ending intensity
          ENDING_NORMALIZED_INTENSITY = NORMALIZED_INTENSITY[n()],
          # Overall change in intensity from start to max
          START_TO_MAX_INTENSITY = MAX_NORMALIZED_INTENSITY - STARTING_NORMALIZED_INTENSITY,
          # For pointing out which frame contains the max intensity
          MAX_INTENSITY_TIME = ifelse(MAX_NORMALIZED_INTENSITY == NORMALIZED_INTENSITY, TIME, NA),
          # Get first frame to reach track max intensity in case of duplicates
          MAX_INTENSITY_TIME = min(MAX_INTENSITY_TIME, na.rm = TRUE)
        )
      
      # Write file
      DestinationPath <- file.path(CellPath, "Analysis.csv.gz")
      fwrite(IntensityTables, DestinationPath, row.names = F, na = "")
      return(DestinationPath)
    }
  }, error = function(e){print(paste("ERROR with CombineCellTables CellX =", CellX))})
}
Cells <- unique(PairedCalibrations$CELL_PATH)
CellAnalysis <- mclapply(1:NROW(Cells), CombineCellTables)
CellAnalysis <- unlist(CellAnalysis)
CellAnalysis <- CellAnalysis[file.exists(CellAnalysis)]
# 
# # Combine all cell tables
# Images <- dirname(dirname(CellAnalysis))
# Images <- unique(Images)
# CombineImageTables <- function(ImageX){
#   tryCatch({
#     # Get cell table paths
#     CellsList <- CellAnalysis[dirname(dirname(CellAnalysis)) == Images[ImageX]]
#     # Get cell tables
#     CellsList <- lapply(CellsList, fread)
#     CellsList <- rbindlist(CellsList, fill = TRUE)
#     # Save combined tables
#     DestinationPath <- file.path(Images[ImageX], "Analysis.csv.gz")
#     fwrite(CellsList, DestinationPath, row.names = F, na = "")
#     
#     # Make image summary for future analysis
#     CellSummary <-
#       CellsList %>% 
#       group_by(
#         LIGAND_DENSITY_CAT,
#         COHORT,
#         IMAGE,
#         CELL,
#         PROTEIN
#       ) %>% 
#       mutate(
#         CELL_AREA = CELL_AREA*CALIBRATION_UM,
#         FRAMES = max(TIME_SINCE_LANDING*FRAME_RATE),
#       ) %>% 
#       group_by(
#         LIGAND_DENSITY_CAT,
#         COHORT,
#         IMAGE,
#         CELL,
#         PROTEIN,
#         UNIVERSAL_TRACK_ID
#       ) %>% 
#       mutate(
#         SPOTS = n()
#       ) %>% 
#       filter(
#         FRAMES_ADJUSTED == 0
#       ) %>% 
#       group_by(
#         LIGAND_DENSITY_CAT,
#         COHORT,
#         IMAGE,
#         CELL,
#         PROTEIN,
#         CELL_AREA
#       ) %>% 
#       summarize(
#         SPOTS = n(),
#         # SPOTS_PER_FRAME = SPOTS/FRAMES,
#         # # SPOTS_PER_AREA_PER_FRAME = SPOTS_PER_FRAME/CELL_AREA,
#         LIFETIME = mean(LIFETIME, na.rm = T),
#         STARTING_NORMALIZED_INTENSITY = mean(STARTING_NORMALIZED_INTENSITY, na.rm = T),
#         MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY, na.rm = T),
#         START_TO_MAX_INTENSITY = mean(START_TO_MAX_INTENSITY, na.rm = T)
#       )
#     
#     ProteinSummary <-
#       CellSummary %>% 
#       group_by(
#         LIGAND_DENSITY_CAT,
#         COHORT,
#         IMAGE,
#         PROTEIN
#       ) %>% 
#       summarize(
#         CELLS = n(),
#         SPOTS = sum(SPOTS),
#         SPOTS_PER_FRAME = mean(SPOTS_PER_FRAME),
#         SPOTS_PER_AREA_PER_FRAME = mean(SPOTS_PER_AREA_PER_FRAME),
#         LIFETIME = mean(LIFETIME, na.rm = T),
#         STARTING_NORMALIZED_INTENSITY = mean(STARTING_NORMALIZED_INTENSITY, na.rm = T),
#         MAX_NORMALIZED_INTENSITY = mean(MAX_NORMALIZED_INTENSITY, na.rm = T),
#         START_TO_MAX_INTENSITY = mean(START_TO_MAX_INTENSITY, na.rm = T)
#       )
#     
#     return(ImageSummary)
#   }, error = function(e){print(paste("ERROR with CombineImageTables ImageX =", ImageX))})
# }
# Images <- lapply(1:NROW(Images), CombineImageTables)
# 
# # Make table of analyzed images of a date

