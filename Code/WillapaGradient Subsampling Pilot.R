######################
# Evaluating CAMA subsampling for biometrics 11/2025
# Emily Grason
######################


########### Packages
library(ggplot2)
library(dplyr)


########### Data
wb.25 <- read.csv("~/Documents/GitHub/CAMA-Subsampling/Data/Willapa Gradient Survey 2025.CrabSize.csv")
wb.25$Site_Name <- as.factor(wb.25$Site_Name)

########### Plots

######Full Sample
ggplot(wb.25, aes(x = "", y = CW_mm)) +
  geom_violin() + 
  geom_jitter(aes(width = 0.2, 
              alpha = 0.5)) + 
  theme_bw() +
  theme_bw(base_size = 16) +
  guides(alpha = "none") +
  ylab("Crab Size (mm)") 

mean(wb.25$CW_mm < 35)
mean(wb.25$CW_mm >70)



#########Compare Subsamples to mean
subsample_stats_multi <- function(data, subsample_sizes, n_reps) {
  data <- as.numeric(data)
  
  # Total ("full") statistics
  data <- as.numeric(data)
  total_mean <- mean(data, na.rm = TRUE)
  total_max <- max(data, na.rm = TRUE)
  total_min <- min(data, na.rm = TRUE)
  total_var  <- var(data, na.rm = TRUE)
  total_small <- mean(data < 35, na.rm = TRUE)
  total_large <- mean(data >70, na.rm = TRUE)
  
  
  # Loop over each subsample size
  results_list <- lapply(subsample_sizes, function(size) {
    
    # Replicate subsamples
    replicate_results <- replicate(n_reps, {
      subsample <- sample(data, size = size, replace = FALSE)
      
      # Compute subsample stats
      s_mean <- mean(subsample)
      s_max  <- max(subsample)
      s_min  <- min(subsample)
      s_var  <- var(subsample)
      s_small <- mean(subsample < 35)
      s_large <- mean(subsample > 70)
      
      # Return vector of differences
      c(
        mean_diff = s_mean - total_mean,
        max_diff  = s_max  - total_max,
        min_diff  = s_min  - total_min,
        var_diff  = s_var  - total_var,
        small_dif = s_small - total_small,
        large_dif = s_large - total_large
      )
    })
    
    # Convert each replicate to a row
    replicate_df <- as.data.frame(t(replicate_results))
    replicate_df$subsample_size <- size
    replicate_df$replicate <- seq_len(n_reps)
    
    replicate_df
  })
  
  # Combine all results
  results_df <- do.call(rbind, results_list)
  
  # Return tidy output
  return(results_df)
}

#Run sampling function on data
ntot <- nrow(wb.25)
psamp <- c(0.01, 0.05, 0.10, 0.20, 0.25, 0.50, 0.75, 1)
nsamp <- as.integer(ntot * psamp)
n.iter <- 1000

pooled.sample <- subsample_stats_multi(wb.25$CW_mm, nsamp, n.iter)


## Plots
ggplot(delta.mu, aes(x = subsample_size/ntot, y = mean_diff)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Mean (mm)")

ggplot(delta.mu, aes(x = subsample_size/ntot, y = max_diff)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Max (mm)")

ggplot(delta.mu, aes(x = subsample_size/ntot, y = min_diff)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True min (mm)")

ggplot(delta.mu, aes(x = subsample_size/ntot, y = var_diff)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Variance (mm)")

ggplot(delta.mu, aes(x = subsample_size/ntot, y = small_dif)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Proportion of Small (mm)")

ggplot(delta.mu, aes(x = subsample_size/ntot, y = large_dif)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Proportion of Large (mm)")


# Calculate percentages outside range
summary_df <- delta.mu %>%
  group_by(as.factor(subsample_size)) %>%
  summarize(
    total_count = n(),
    outside_range_count = sum(difference < -3 | difference > 3),
    percent_outside = outside_range_count / total_count
  )



######Split by Site
ggplot(wb.25, aes(x = Site_Name, y = CW_mm, group = Site_Name,)) +
  geom_violin() + 
  geom_jitter(aes(x = Site_Name, y = CW_mm,
                  width = 0.2, 
                  alpha = 0.5, 
                  group = Site_Name,
                  color = Site_Name)) + 
  theme_bw() +
  theme_bw(base_size = 16) +
  guides(alpha = "none") +
  ylab("Crab Size (mm)") +
  guides(color=guide_legend("Site"))