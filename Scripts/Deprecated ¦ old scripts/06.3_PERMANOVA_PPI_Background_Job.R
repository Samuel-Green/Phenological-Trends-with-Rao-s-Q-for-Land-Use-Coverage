library(vegan)
library(dplyr)
library(here)

Knepp_DF_YoI.PPI <- readRDS("~/user_storage/phenological-trends-with-rao-s-q-for-land-use-coverage/Results 📈📉/Knepp Estate/Knepp_DF_YoI-PPI.rds")

Results <- here::here("Results 📈📉")
Knepp_Results <- file.path(Results, "Knepp Estate")

## Subsample the dataframe
# Take a sample of 10,000 rows to prevent distance-matrix memory crashes
# Many land cover classes are ecologically insignificant in this dataset and reduce model performance
# Therefore, I'm specifying the important land cover types and then subsetting to include those only

target_classes <- c(
  "Broadleaved woodland", 
  "Coniferous woodland",
  "Arable and horticulture", 
  "Improved grassland", 
  "Urban",
  "Suburban" 
) # Add or remove classes as required

# Dynamically calculate the sample size per class

total_target_pixels <- 10000
n_classes <- length(target_classes)

# Use floor() to round down so we don't accidentally ask for a fraction of a pixel

pixels_per_class <- floor(total_target_pixels / n_classes) 

message("Sampling up to ", pixels_per_class, " pixels across ", n_classes, " classes.")

## And now we loop over the dataframe to run the PERMANOVAe upon it

Knepp_PERMANOVA_YoI.PPI <- lapply(
  names(Knepp_DF_YoI.PPI),
  function(y) {
    
    df <- Knepp_DF_YoI.PPI[[y]]
    
    # 3. Filter out the fringe/empty classes
    
    df_filtered <- df[df$LandCover %in% target_classes, ]
    
    # 4. Stratified Subsampling using dynamic values
    
    set.seed(123)
    df_subset <- df_filtered %>%
      group_by(LandCover) %>%
      # slice_sample automatically handles groups that are smaller than pixels_per_class!
      slice_sample(n = pixels_per_class, replace = FALSE) %>% 
      ungroup() %>%
      as.data.frame()
    
    list(
      
      # These are structured in the same way as the PERMANOVAe from the other case studies
      
      Shannon = adonis2(
        df_subset$ShannonH ~ df_subset$LandCover,
        #data = df_subset, # IDK why, but R cannot find ShannonH unless defined by a $
        permutations = 999,
        parallel = parallel::detectCores() - 2
      ),
      
      Rao_Classic = adonis2(
        df_subset$RaosQ_Classic ~ df_subset$LandCover,
        #data = df_subset,
        permutations = 999,
        parallel = parallel::detectCores() - 2
      ),
      
      Rao_TWDTW = adonis2(
        df_subset$RaosQ_TWDTW ~ df_subset$LandCover,
        #data = df_subset,
        permutations = 999,
        parallel = parallel::detectCores() - 2
      )
    )
  }
)

names(Knepp_PERMANOVA_YoI.PPI) <- names(Knepp_DF_YoI.PPI) # Name the list elements

saveRDS(Knepp_PERMANOVA_YoI.PPI, file = file.path(Knepp_Results, "Knepp_PERMANOVA_YoI_PPI.rds"))

# Finally, putting the results into their own dataframe

Knepp_PERMANOVA_Summary.PPI <- do.call(rbind, lapply(
  names(Knepp_PERMANOVA_YoI.PPI),
  function(y) {
    
    res <- Knepp_PERMANOVA_YoI.PPI[[y]]
    
    data.frame(
      Year = y,
      Index = c("Shannon", "Classic Rao", "TWDTW Rao"),
      R2 = c(
        res$Shannon$R2[1],
        res$Rao_Classic$R2[1],
        res$Rao_TWDTW$R2[1]
      ),
      F = c(
        res$Shannon$F[1],
        res$Rao_Classic$F[1],
        res$Rao_TWDTW$F[1]
      ),
      p_value = c(
        res$Shannon$`Pr(>F)`[1],
        res$Rao_Classic$`Pr(>F)`[1],
        res$Rao_TWDTW$`Pr(>F)`[1]
      )
    )
  }
))

saveRDS(Knepp_PERMANOVA_Summary.PPI, file = file.path(Knepp_Results, "Knepp_PERMANOVA_Summary_PPI.rds"))

message("Knepp Estate PPI analysis complete.") # End of script ####