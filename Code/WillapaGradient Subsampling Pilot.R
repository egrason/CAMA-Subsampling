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
subsample_stats_multi <- function(data, subsample_sizes, n_reps, ntrap_vec) {
  data <- as.numeric(data)
  
  # Total ("full") statistics
  total_mean <- mean(data, na.rm = TRUE)
  total_max  <- max(data, na.rm = TRUE)
  total_min  <- min(data, na.rm = TRUE)
  total_var  <- var(data, na.rm = TRUE) * (((length(data)) - 1) / length(data))
  total_small <- mean(data < 35, na.rm = TRUE)
  total_large <- mean(data > 70, na.rm = TRUE)
  total_small_n <- sum(data < 35, na.rm = TRUE)
  total_large_n <- sum(data > 70, na.rm = TRUE)
  
  # Loop over each subsample size (using index!)
  results_list <- lapply(seq_along(subsample_sizes), function(i) {
    
    size <- subsample_sizes[i]
    ntrap_val <- ntrap_vec[i]   # assign corresponding ntrap
    
    replicate_results <- replicate(n_reps, {
      subsample <- sample(data, size = size, replace = FALSE)
      
      s_mean <- mean(subsample)
      s_max  <- max(subsample)
      s_min  <- min(subsample)
      s_var  <- var(subsample) * (((length(subsample)) - 1) / length(subsample))
      s_small <- mean(subsample < 35)
      s_large <- mean(subsample > 70)
      s_small_n <- sum(subsample < 35)
      s_large_n <- sum(subsample > 70)
      
      c(
        mean_diff  = s_mean - total_mean,
        max_diff   = s_max  - total_max,
        min_diff   = s_min  - total_min,
        var_diff   = s_var  - total_var,
        small_diff = s_small - total_small,
        large_diff = s_large - total_large,
        small_diff_n = s_small_n - total_small_n,
        large_diff_n = s_large_n - total_large_n
      )
    })
    
    df <- as.data.frame(t(replicate_results))
    df$subsample_size <- size
    df$ntrap          <- ntrap_val
    df$replicate      <- seq_len(n_reps)
    
    df
  })
  
  do.call(rbind, results_list)
}

#Run sampling function on data
ntot <- nrow(wb.25)
psamp <- c(0.01, 0.05, 0.10, 0.20, 0.25, 0.50, 0.75, 1)
nsamp <- as.integer(ntot * psamp)
ntrap <- c(2, 12, 24, 48, 60, 120, 180, 240)
n.iter <- 1000

pooled.sample <- subsample_stats_multi(wb.25$CW_mm, nsamp, n.iter, ntrap)


## Plots
ggplot(pooled.sample, aes(x = subsample_size/ntot, y = mean_diff)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Mean (mm)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = max_diff)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Max (mm)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = min_diff)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Min (mm)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = var_diff)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Variance (mm)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = small_diff)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Proportion of Small (mm)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = large_diff)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Proportion of Large (mm)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = 100*small_diff_n/240)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Small Crab CPUE (per 100)")

ggplot(pooled.sample, aes(x = subsample_size/ntot, y = 100*large_diff_n/240)) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Large Crab CPUE (per 100)")

#total_small_n <- sum(data < 35, na.rm = TRUE)
#total_large_n <- sum(data > 70, na.rm = TRUE)


ggplot(pooled.sample, aes(x = subsample_size/ntot, 
                          y = 100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
                            (100*(sum(wb.25$CW_mm < 35)/240)) 
                        )) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Percent Difference from True Small CPUE")

ggplot(pooled.sample, aes(x = subsample_size/ntot, 
                          y = 100 * ((100*((sum(wb.25$CW_mm >70)) + large_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm > 70)/240))) /
                            (100*(sum(wb.25$CW_mm > 70)/240)) 
                          )) +
  geom_jitter(alpha = 0.5) + 
  theme_bw(base_size = 16) +
  labs(x = "Proportion of Crabs Measured", y = "Percent Difference from True Large CPUE")



# Calculate percentages outside range
mean_summary <- pooled.sample %>%
  group_by(as.factor(subsample_size)) %>%
  summarize(
    total_count = n(),
    outside_range_count = sum(mean_diff < -3 | mean_diff > 3),
    percent_outside = outside_range_count / total_count
  )



CPUE_small_summary <- pooled.sample %>%
  group_by(as.factor(subsample_size)) %>%
  summarize(
    total_count = n(),
    outside_range_count = sum(
      100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
        (100*(sum(wb.25$CW_mm < 35)/240))  < -5 | 
      100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
        (100*(sum(wb.25$CW_mm < 35)/240))   > 5),
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