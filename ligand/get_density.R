#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

parameters_path = args[1]
ligand_images_path = args[2]
output_path = args[3]

library(ggplot2)
library(dplyr)
library(ggdark)
library(parallel)
library(XML)
library(fitdistrplus)
library(data.table)
library(scales)
library(ff)

# Output
output_path <- file.path(output_path, "Ligand")
if(!file.exists(output_path)){
  dir.create(output_path)
}

# Import table
LigandList <- file.path(parameters_path, "ligand.csv")
LigandList <- read.csv(LigandList)

GetDensityFx <- function(ImageX){
  tryCatch({
    # Get paths
    FolderName <- tools::file_path_sans_ext(LigandList$image[ImageX])
    DataPath <- file.path(ligand_images_path, FolderName)
    MetadataTable <- file.path(DataPath, "metadata.csv")
    XMLTable <- paste0(LigandList$protein_name[ImageX], ".xml")
    XMLTable <- file.path(DataPath, XMLTable)
    
    # Get variables
    DILUTION_FACTOR = LigandList$dilution[ImageX]
    
    MetadataTable <- read.csv(MetadataTable)
    
    PIXEL_SIZE = which(MetadataTable$parameter=="calibration_um")
    PIXEL_SIZE = MetadataTable$value[PIXEL_SIZE]
    PIXEL_SIZE = as.numeric(PIXEL_SIZE)
    
    X_SIZE = which(MetadataTable$parameter=="width")
    X_SIZE = MetadataTable$value[X_SIZE]
    X_SIZE = as.numeric(X_SIZE)
    X_SIZE = X_SIZE*PIXEL_SIZE
    
    Y_SIZE = which(MetadataTable$parameter=="height")
    Y_SIZE = MetadataTable$value[Y_SIZE]
    Y_SIZE = as.numeric(Y_SIZE)
    Y_SIZE = Y_SIZE*PIXEL_SIZE
    
    # Get XML
    xml_data <- xmlParse(XMLTable)
    xml_data <- xmlToList(xml_data)
    
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
    
    # Get number of spots per frame
    LigandTable <-
      SpotsTable %>% 
      group_by(
        FRAME
      ) %>% 
      summarize(
        N = n()
      ) %>% 
      mutate(
        FRAME = as.numeric(FRAME),
        FRAME = FRAME + 1
      )
    
    # Get upper boundary
    UpperLimitFx <- function(FrameX){
      
      LigandTable <-
        LigandTable %>%
        filter(
          FRAME <= FrameX
        )
      
      R2 = cor(LigandTable$N, LigandTable$FRAME)
      R2 = as_tibble(R2)
      return(R2)
    }
    MAX_FRAME = max(LigandTable$FRAME)
    Results <- mclapply(25:MAX_FRAME, UpperLimitFx)
    Results <- rbindlist(Results)
    # Get min R
    MIN_R2 = min(Results)
    ROW = which(Results$value == MIN_R2)
    UPPER_FRAME_LIMIT = ROW + 25
    
    LowerLimitFx <- function(FrameX){
      
      LigandTable <-
        LigandTable %>%
        filter(
          FRAME <= UPPER_FRAME_LIMIT,
          FRAME >= FrameX
        )
      
      R2 = cor(LigandTable$N, LigandTable$FRAME)
      R2 = as_tibble(R2)
      return(R2)
    }
    Results <- mclapply(1:(UPPER_FRAME_LIMIT-1), LowerLimitFx)
    Results <- rbindlist(Results)
    # Get min R
    R2 = min(Results)
    LOWER_FRAME_LIMIT = which(Results$value == R2)
    # Redefine ligand density
    LigandTable <-
      LigandTable %>% 
      filter(
        FRAME >= LOWER_FRAME_LIMIT,
        FRAME <= UPPER_FRAME_LIMIT
      )
    
    # Get density
    fit <- lm(log(N) ~ FRAME, data = LigandTable)
    INTERCEPT <- fit$coefficients[1]
    INTERCEPT <- exp(INTERCEPT)
    LIGAND_DENSITY =  as.numeric(INTERCEPT)*DILUTION_FACTOR/X_SIZE/Y_SIZE
    
    # Draw fit line
    FRAMES <- LigandTable$FRAME
    N <- exp(predict(fit,list(FRAME=FRAMES)))
    FittedLigandTable <- as.data.frame(N)
    FittedLigandTable$FRAME <- FRAMES
    
    # Plot
    ggplot() +
      geom_path(
        data = FittedLigandTable,
        aes(
          x = FRAME,
          y = N
        ),
        color ="red"
      ) +
      geom_point(
        data = LigandTable,
        aes(
          x = FRAME,
          y = N
        )
      ) +
      scale_y_continuous(
        # Log transform axis. The secondary axis reinforces that it's in log as it'll say 2^x
        trans = "log2",
        sec.axis = sec_axis(
          trans = ~.,
          breaks = trans_breaks("log2", function(x) 2^x),
          labels = trans_format("log2", math_format(2^.x))
        )
      ) +
      dark_theme_classic() +
      labs(
        x = "Frame",
        y = "Spots",
        title = paste(round(INTERCEPT),"molecules, R =", round(R2, 2), "\nLigand Density =", round(LIGAND_DENSITY, 2), "mol µm^-2")
      ) +
      ggsave(
        file.path(DataPath, "Fit.pdf"),
        height = 3,
        width = 4
      )
    
    # Revert ggdark
    ggdark::invert_geom_defaults()
    
    # Create output table
    OutputTable <- NULL
    OutputTable$DATE = substr(FolderName, 0, 8)
    OutputTable$IMAGE = FolderName
    OutputTable$LIGAND_DENSITY = LIGAND_DENSITY
    
    OutputTable$LIGAND = LigandList$ligand[ImageX]
    OutputTable$DILUTION = LigandList$dilution[ImageX]
    OutputTable$PROTEIN_NAME = LigandList$protein_name[ImageX]
    
    OutputTable$MOLECULES = as.numeric(INTERCEPT)
    OutputTable$R2 = R2
    OutputTable$UPPER_FRAME_LIMIT = UPPER_FRAME_LIMIT
    OutputTable$LOWER_FRAME_LIMIT = LOWER_FRAME_LIMIT
    
    
    # Save it to csv
    write.csv(OutputTable, file.path(DataPath, "ligand_density.csv"), row.names = F)

    return(FolderName)
    
  }, error = function(e){print(paste("ERROR with GetDensityFx ImageX =", ImageX))})
}
Results <- mclapply(1:NROW(LigandList), GetDensityFx)

# Move files
for(ImageX in 1:NROW(Results)){
  old_path = file.path(ligand_images_path, Results[ImageX])
  new_path = file.path(output_path, Results[ImageX])
  file.move(old_path, new_path)
}


