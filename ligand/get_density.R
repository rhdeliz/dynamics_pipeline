DILUTION_FACTOR = 1000
Y_SIZE = 88
X_SIZE = 88
INPUT_TABLE_PATH = "/Users/u_deliz/Desktop/IL-1.csv"

library(fitdistrplus)
library(ggplot2)
library(dplyr)
library(ggdark)

# Import table
LigandTable <- read.csv(INPUT_TABLE_PATH)
setwd(dirname(INPUT_TABLE_PATH))
names(LigandTable) <- "FRAME"

# Get number of spots per frame
LigandTable <-
  LigandTable %>%
  group_by(
    FRAME
  ) %>% 
  summarize(
    N = n()
  ) %>% 
  mutate(
    FRAME = FRAME + 1
  )

FitFx <- function(FrameX){
  
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
Results <- mclapply(25:MAX_FRAME, FitFx)
Results <- rbindlist(Results)
# Get min R
MIN_R2 = min(Results)
ROW = which(Results$value == MIN_R2)
FRAME_LIMIT = ROW + 25

LigandTable <-
  LigandTable %>% 
  filter(
    FRAME <= FRAME_LIMIT
  )

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
    title = paste(round(INTERCEPT),"molecules, R =", round(R2, 2), "\nLigand Density =", round(LIGAND_DENSITY, 2))
  ) +
  ggsave(
    "Fit.pdf",
    height = 3,
    width = 4
  )

# Revert ggdark
ggdark::invert_geom_defaults()

# Create output table
OutputTable <- NULL
OutputTable$LIGAND_DENSITY = LIGAND_DENSITY
OutputTable$MOLECULES = as.numeric(INTERCEPT)
OutputTable$R2 = R2
OutputTable$FRAME_LIMIT = FRAME_LIMIT
# Save it to csv
write.csv(OutputTable, "LigandDensity.csv", row.names = F)