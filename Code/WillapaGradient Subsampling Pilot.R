######################
# Evaluating CAMA subsampling for biometrics 11/2025
# Emily Grason
######################


########### Packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)


########### Data
wb.25 <- read.csv("~/Documents/GitHub/CAMA-Subsampling/Data/Willapa Gradient Survey 2025.CrabSize.csv")
wb.25$Site_Name <- as.factor(wb.25$Site_Name)


########### Example Plots

n <- 701

# --- 1. Unimodal normal distribution ---
normal_data <- data.frame(
  CW_mm = rnorm(n, mean = 53, sd = 2)   # sd controls variance
)

# --- 2. Bimodal distribution ---
bimodal_data <- data.frame(
  CW_mm = c(
    rnorm(n/2, mean = 73, sd = 5),
    rnorm(n/2, mean = 33, sd = 5)
  )
)

# --- 3. Uniform distribution ---
uniform_data <- data.frame(
  CW_mm = runif(n, min = 53 - 20, max = 53 + 20)
)

make_plot <- function(df, title) {
  ggplot(df, aes(x = "", y = CW_mm)) +
    geom_violin() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    geom_hline(
      yintercept = mean(df$CW_mm),
      color = "black",
      linewidth = 1,
      linetype = "solid"
    ) +
    ylim(0, 100) +    
    theme_bw(base_size = 16) +
    ylab("Crab Size (mm)")
}

p1 <- make_plot(normal_data)
p2 <- make_plot(bimodal_data)
p3 <- make_plot(uniform_data)

p1
p2
p3


###########################
#########Pooled Sample (single bucket)
###########################

######Full Sample Plot - Pooled
ggplot(wb.25, aes(x = "", y = CW_mm)) +
  geom_violin() + 
  geom_jitter(width = 0.2, 
              alpha = 0.5) + 
  geom_hline(
    yintercept = mean(wb.25$CW_mm, na.rm = TRUE),
    color = "black",
    linewidth = 1,
    linetype = "solid"
  ) +
  theme_bw() +
  theme_bw(base_size = 16) +
  guides(alpha = "none") +
  ylab("Crab Size (mm)") 

n_small <- mean(wb.25$CW_mm < 35)
n_large <- mean(wb.25$CW_mm >70)

set.seed(123)   # optional for reproducibility

wb.25$highlight <- FALSE
wb.25$highlight[sample(nrow(wb.25), 100)] <- TRUE

ggplot(wb.25, aes(x = "", y = CW_mm)) +
  geom_violin() +
  geom_jitter(
    aes(color = highlight),
    width = 0.2,
    alpha = 0.5
  ) +
  scale_color_manual(
    values = c("TRUE" = "red", "FALSE" = "black"),
    guide = "none"
  ) +
  theme_bw(base_size = 16) +
  ylab("Crab Size (mm)")

###########################
#########Compare Subsamples of varying sizes to Total Sample
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
psamp <- c(0.01, 0.05, 0.10, 100/701, 0.20, 0.25, 0.50, 0.75, 1)
nsamp <- as.integer(ntot * psamp)
ntrap <- c(2, 12, 24, 34, 48, 60, 120, 180, 240)
n.iter <- 1000

pooled.sample <- subsample_stats_multi(wb.25$CW_mm, nsamp, n.iter, ntrap)

###########################
## Plots

## Absolute difference plots
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


###########################
# Calculate percentages outside range

##Single factor, Single Boundary - Simple
mean_summary <- pooled.sample %>%
  group_by(as.factor(subsample_size)) %>%
  summarize(
    total_count = n(),
    outside_range_count = sum(small_diff < -0.1 | small_diff > 0.1),
    percent_outside = outside_range_count / total_count
  )

##Single factor, Single Boundary - CPUE
CPUE_small_summary <- pooled.sample %>%
  group_by(as.factor(subsample_size)) %>%
  summarize(
    total_count = n(),
    outside_range_count = sum(
      100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
        (100*(sum(wb.25$CW_mm < 35)/240))  < -b | 
      100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
        (100*(sum(wb.25$CW_mm < 35)/240))   > b),
    percent_outside = outside_range_count / total_count
  )


##Single factor - Multiple boundaries - CPUE
boundaries <- c(5, 10, 25)
out_of_range <- map_dfr(boundaries, function(b) {
  pooled.sample %>%
    group_by(subsample_size) %>%
    summarise(
      boundary = b,
      n_outside = sum(
        100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
          (100*(sum(wb.25$CW_mm < 35)/240))  < -b | 
          100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
          (100*(sum(wb.25$CW_mm < 35)/240))   > b),
      n_total   = n(),
      prop_outside = n_outside / n_total
    )
})

##Single Factor - Variable range - Simple
#"mean_diff", "min_diff", "max_diff", "var_diff", "small_diff", 
#"large_diff" 

# Choose variable
variable_to_test <- "large_diff_n"

# Continuous boundaries
boundary_seq <- seq(1, 50, by = 1)

# Compute prop_outside for each boundary × subsample_size
#out_of_range_surface <- tibble(boundary = boundary_seq) %>%
#  mutate(
#    results = map(boundary, function(bnd) {
#      pooled.sample %>%
#        group_by(subsample_size) %>%
#        summarise(
#          prop_outside = mean(get(variable_to_test) < -bnd | get(variable_to_test) > bnd),
#          .groups = "drop"
#        )
#    })
#  ) %>%
#  unnest(results)




##Single Factor - Variable range - CPUE
# CPUE vars: "small_diff_n", "large_diff_n"

# Choose variable
variable_to_test <- "small_diff_n"

# Continuous boundaries
boundary_seq <- seq(1, 50, by = 1)

out_of_range_surface <- tibble(boundary = boundary_seq) %>%
  mutate(
    results = map(boundary, function(bnd) {
      pooled.sample %>%
        group_by(subsample_size) %>%
        summarise(
          prop_outside = mean(
            100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
              (100*(sum(wb.25$CW_mm < 35)/240))  < -bnd | 
              100 * ((100*((sum(wb.25$CW_mm < 35)) + small_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm < 35)/240))) /
              (100*(sum(wb.25$CW_mm < 35)/240))   > bnd),
          .groups = "drop"
        )
    })
  ) %>%
  unnest(results)


#Plot for either of above
total_crabs <- 701  # total number of crabs
subsample_sizes <- sort(unique(out_of_range_surface$subsample_size))  # all subsample sizes used

ggplot(out_of_range_surface, aes(x = subsample_size / total_crabs, 
                                 y = boundary, 
                                 fill = prop_outside)) +

  # Vertical guide lines at each subsample size
  geom_vline(
    xintercept = subsample_sizes / total_crabs,
    color = "grey80",
    linewidth = 0.3
  ) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", 
                       name = "Proportion\nOutside",
                       limits = c(0,1)) +
  
  # Primary x-axis: proportion of crabs measured
  scale_x_continuous(
    name = "Proportion of Crabs Measured",
    
    # Secondary axis: show only the actual subsample sizes used
    sec.axis = sec_axis(
      ~ . * total_crabs, 
      breaks = subsample_sizes
    )
  ) +
  
  labs(
    y = "Boundary (Percent Difference)",
    title = paste("Proportion Outside ±Boundary for", variable_to_test)
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x.bottom = element_text(color = "black"),
    axis.text.x.top = element_text(color = "black", angle = 45, hjust = 1, vjust = 1)
  )



###########################
######### By Trap (sequential)
###########################



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