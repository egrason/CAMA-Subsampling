######################
# Evaluating CAMA subsampling for biometrics 11/2025
# Emily Grason
######################


########### Packages
library(ggplot2)


########### Data
wb.25 <- read.csv("~/Documents/GitHub/CAMA-Subsampling/Data/Willapa Gradient Survey 2025.CrabSize.csv")


########### Plots


######Full Sample
ggplot(wb.25, aes(x = "", y = CW_mm)) +
  geom_violin() + 
  geom_jitter(width = 0.2, alpha = 0.5) + 
  theme_bw()


ntot <- nrow(wb.25)
psamp <- c(0.01, 0.05, 0.10, 0.20, 0.25, 0.50, 0.75, 1)
nsamp <- as.integer(ntot * psamp)
n.iter <- 1000


# Function to subsample a dataset for multiple subsample sizes
subsample_differences_multi <- function(data, subsample_sizes, n_reps) {
  data <- as.numeric(data)
  total_mean <- mean(data, na.rm = TRUE)
  
  # For each subsample size, compute replicate differences
  results <- lapply(subsample_sizes, function(size) {
    differences <- replicate(n_reps, {
      subsample <- sample(data, size = size, replace = FALSE)
      subsample_mean <- mean(subsample)
      subsample_mean - total_mean
    })
    
    data.frame(
      subsample_size = size,
      difference = differences
    )
  })
  
  # Combine results into a single data frame
  results_df <- do.call(rbind, results)
  return(results_df)
}

#Run Sampling function on data
delta.mu <- subsample_differences_multi(wb.25$CW_mm, nsamp, n.iter)

ggplot(delta.mu, aes(x = subsample_size/ntot, y = difference)) +
  geom_point(alpha = 0.5) + 
  theme_bw() +
  labs(x = "Proportion of Crabs Measured", y = "Difference from True Mean (mm)")
