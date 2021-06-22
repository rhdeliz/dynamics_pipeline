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
# Break tables up by frame``
ScreenedFrameTables <-
  ScreenedTables %>% 
  group_split(
    FRAME
  )

# Compute cost (i.e., find nearest neighbors for each puncta)
CostFx <- function(FrameTable){
  tryCatch({
    # Get distance
    Radius <- max(FrameTable$PUNCTA_DIAMETER, na.rm = T)/2
    Radius = Radius*3
    # Get minimum association time
    AssocTime <- max(FrameTable$TRACKMATE_FRAME_GAP, na.rm = T) + 1
    
    # Get number of proteins
    Proteins <- unique(FrameTable$PROTEIN)
    # Get combinations
    Combinations <- combn(Proteins, 2)
    Tests <- 1:NCOL(Combinations)
    
    FindNeighbors <- function(PairX){
      tryCatch({
        # Proteins to run
        ReferenceProtein <- Combinations[1,PairX]
        QueryProtein <- Combinations[2,PairX]
        # Get tables
        ReferenceProtein <- FrameTable %>% filter(PROTEIN == ReferenceProtein)
        QueryProtein <- FrameTable %>% filter(PROTEIN == QueryProtein)
        ReferenceProteinCoordinates <- ReferenceProtein %>% select(POSITION_X, POSITION_Y)
        QueryProteinCoordinates <- QueryProtein %>% select(POSITION_X, POSITION_Y)
        
        Distance <- RANN::nn2(QueryProteinCoordinates,  ReferenceProteinCoordinates, searchtype = c("radius"), radius = Radius,k = 1)
        
        GetSpotAndDistanceFx <- function(DistanceX){
          tryCatch({
            QuerySpotName <- Distance$nn.idx[DistanceX]
            if(QuerySpotName != 0){
              Result <- NULL
              Result$REFERENCE_SPOT <- ReferenceProtein$UNIVERSAL_SPOT_ID[DistanceX]
              Result$REFERENCE_TRACK <- ReferenceProtein$UNIVERSAL_TRACK_ID[DistanceX]
              Result$QUERY_SPOT <- QueryProtein$UNIVERSAL_SPOT_ID[QuerySpotName]
              Result$QUERY_TRACK <- QueryProtein$UNIVERSAL_TRACK_ID[QuerySpotName]
              
              Result$COLOCALIZATION_DISTANCE <- Distance$nn.dists[DistanceX]
              Result$ASSOCIATION_TIME_THRESHOLD <- AssocTime
              # Result$UNIVERSAL_SPOT_ID <- c(QueryProtein$UNIVERSAL_SPOT_ID[QuerySpotName], ReferenceProtein$UNIVERSAL_SPOT_ID[DistanceX])
              # Result$COLOCALIZATION_SPOT <- c(ReferenceProtein$UNIVERSAL_SPOT_ID[DistanceX], QueryProtein$UNIVERSAL_SPOT_ID[QuerySpotName])
              # Result$DISTANCE <- rep(Distance$nn.dists[DistanceX], 2)
              Result <- as.data.frame(Result)
              return(Result)
            }
          }, error = function(e){print(paste("ERROR with GetSpotAndDistanceFx. DistanceX =", DistanceX))})
        }
        Distance <- lapply(1:NROW(Distance$nn.idx), GetSpotAndDistanceFx)
        Distance <- Distance[(which(sapply(Distance,is.list), arr.ind=TRUE))]
        Distance <- rbindlist(Distance)
        return(Distance)
      }, error = function(e){print(paste("ERROR with FindNeighbors. PairX = ", PairX))})
    }
    TableDistances <- lapply(Tests, FindNeighbors)
    TableDistances <- rbindlist(TableDistances)
    
    # Add cell identifier
    if(NROW(TableDistances) > 0){
      TableDistances$IMAGE_CELL = paste(FrameTable$IMAGE[1], FrameTable$CELL[1], sep = "...")
      TableDistances$RELATIVE_PATH = FrameTable$RELATIVE_PATH[1]
    }
    
    return(TableDistances)
  }, error = function(e){print(paste("ERROR with PairingFx"))})
}
CostTable <- mclapply(ScreenedFrameTables, CostFx)
CostTable <- CostTable[(which(sapply(CostTable,is.list), arr.ind=TRUE))]
CostTable <- rbindlist(CostTable, fill = TRUE)
# Split by cell
CostTable <-
  CostTable %>% 
  group_split(
    IMAGE_CELL
  )

# Compute edges and groups
EdgeFx <- function(CellTable){
  tryCatch({
    # Get edges (i.e., the colocalized spot pair)
    Edges <-
      CellTable %>% 
      group_by(
        REFERENCE_TRACK,
        QUERY_TRACK
      ) %>%
      mutate(
        ASSOCIATION_TIME = n()
      ) %>% 
      filter(
        ASSOCIATION_TIME > ASSOCIATION_TIME_THRESHOLD
      ) %>% 
      select(
        REFERENCE_TRACK,
        QUERY_TRACK
      ) %>% 
      distinct()
    
    # Get group numbers 
    names(Edges) <- c("from", "to")
    SubGroups <- igraph::graph_from_data_frame(Edges, directed = F)
    SubGroups <- igraph::clusters(SubGroups)$membership
    SubGroups <- as.data.frame(SubGroups)
    
    # Create new table from result
    Groups <- NULL
    Groups$UNIVERSAL_TRACK_ID <- rownames(SubGroups)
    Groups$COLOCALIZATION_GROUP <- SubGroups$SubGroups
    Groups <- as_tibble(Groups)
    remove(SubGroups, Edges)
    
    # Combine
    Groups <- merge(CellTable, Groups, by.x = "REFERENCE_TRACK", by.y = "UNIVERSAL_TRACK_ID")
    
    # Split reference and query
    GroupsReference <-
      Groups %>%
      mutate(
        UNIVERSAL_SPOT_ID = REFERENCE_SPOT,
        COLOCALIZATION_SPOT = QUERY_SPOT
      ) %>% 
      select(-c(
        REFERENCE_TRACK,
        REFERENCE_SPOT,
        QUERY_TRACK,
        QUERY_SPOT
      ))
    GroupsQuery <-
      Groups %>%
      mutate(
        UNIVERSAL_SPOT_ID = QUERY_SPOT,
        COLOCALIZATION_SPOT = REFERENCE_SPOT
      ) %>% 
      select(-c(
        REFERENCE_TRACK,
        REFERENCE_SPOT,
        QUERY_TRACK,
        QUERY_SPOT
      ))
    # Combine
    Groups <- bind_rows(GroupsReference, GroupsQuery)
    
    return(Groups)
  }, error = function(e){print(paste("ERROR with EdgeFx"))})
}
Edges <- mclapply(CostTable, EdgeFx)
Edges <- Edges[(which(sapply(Edges,is.list), arr.ind=TRUE))]
Edges <- rbindlist(Edges, fill = TRUE)

# Save edges data
Edges <-
  Edges %>% 
  group_split(
    IMAGE_CELL
  )

# Save tables
SaveFx <- function(EdgesCellTable){
  tryCatch({
    EdgesCellTable$IMAGE_CELL = NULL
    # Create file path
    save_path = dirname(EdgesCellTable$RELATIVE_PATH[1])
    save_path = file.path(extraction_path, save_path, "colocalization_coordinates.csv.gz")
    # Save
    fwrite(EdgesCellTable, save_path, row.names = F, na = "")
    return(save_path)
  }, error = function(e){print(paste("ERROR with SaveFx"))})
}
lapply(Edges, SaveFx)

# Move images that needed colocalization
ImagePath <- ColocalizationNeeded$image
for(Image in ImagePath){
  old_path = Image
  Cohort = basename(dirname(Image))
  Image = basename(Image)
  new_path = file.path(colocalized_path, Cohort, Image)
  file.move(old_path, new_path)
}
