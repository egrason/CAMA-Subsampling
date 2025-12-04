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
variable_to_test <- "var_diff"

# Continuous boundaries
boundary_seq <- seq(1, 50, by = 1)

# Compute prop_outside for each boundary × subsample_size
out_of_range_surface <- tibble(boundary = boundary_seq) %>%
  mutate(
    results = map(boundary, function(bnd) {
      pooled.sample %>%
        group_by(subsample_size) %>%
        summarise(
          prop_outside = mean(get(variable_to_test) < -bnd | get(variable_to_test) > bnd),
          .groups = "drop"
        )
    })
  ) %>%
  unnest(results)




##Single Factor - Variable range - CPUE
# CPUE vars: "small_diff_n", "large_diff_n"

# Choose variable
variable_to_test <- "large_diff_n"

# Continuous boundaries
boundary_seq <- seq(1, 50, by = 1)

out_of_range_surface <- tibble(boundary = boundary_seq) %>%
  mutate(
    results = map(boundary, function(bnd) {
      pooled.sample %>%
        group_by(subsample_size) %>%
        summarise(
          prop_outside = mean(
            100 * ((100*((sum(wb.25$CW_mm >70)) + large_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm >70)/240))) /
              (100*(sum(wb.25$CW_mm >70)/240))  < -bnd | 
              100 * ((100*((sum(wb.25$CW_mm >70)) + large_diff_n) / ntrap) - (100*(sum(wb.25$CW_mm >70)/240))) /
              (100*(sum(wb.25$CW_mm >70)/240))   > bnd),
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


#Plot just for 14.3%/100 crabs

target_n <- 100

filtered_surface <- out_of_range_surface %>%
  filter(subsample_size == target_n)

filtered_surface <- filtered_surface %>%
  mutate(prop_pct = prop_outside * 100)

# Find the boundary where percent = 5%
bnd_5pct <- filtered_surface$boundary[which.min(abs(filtered_surface$prop_pct - 5))]

ggplot(filtered_surface,
       aes(y = boundary, fill = prop_pct)) +
  
  # filled ribbon spanning fake x-range
  geom_ribbon(aes(xmin = -1, xmax = 1), color = NA) +
  
  # color ramp fixed 0–100%
  scale_fill_viridis_c(
    option = "plasma",
    name = "% Outside",
    limits = c(0, 100)
  ) +
  
  # 5% horizontal line
  geom_hline(yintercept = bnd_5pct, 
             linetype = "dashed", 
             color = "orange", 
             linewidth = 0.7) +
  
  scale_x_continuous(expand = c(0, 0)) +
  labs(
    title = "Proportion Outside Boundary at 14.3% Subsample",
    y = "Boundary (mm Difference)",
    x = NULL
  ) +
  
  theme_bw(14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

###########################
######### By Trap (sequential)
###########################

wb.25$Trap_Process_Order <- as.numeric(wb.25$Day+ wb.25$TrapOrder)

target_n <- 100   # target min number of crabs

wb_first100 <- wb.25 %>%
  group_by(Trap_Process_Order) %>%
  mutate(group_n = n()) %>%              # size of each group
  ungroup() %>%
  distinct(Trap_Process_Order, group_n) %>%  # one row per group
  mutate(cum_n = cumsum(group_n)) %>%        # cumulative count
  filter(cum_n <= target_n |                  # keep all full groups before threshold
           lag(cum_n, default = 0) < target_n) %>%  # keep first group that pushes past target
  inner_join(wb.25, by = "Trap_Process_Order")   # get all rows from selected groups

mean(wb_first100$CW_mm)
max(wb_first100$CW_mm)
min(wb_first100$CW_mm)
var(wb_first100$CW_mm)

n_first_small <- sum(wb_first100$CW_mm < 35)
n_first_large <- sum(wb_first100$CW_mm >70)

n_first_small *100 / (240 * (100/701))


###########################
######### By Trap (random)
###########################

sample_whole_traps <- function(df, group_col, target_n, seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  group_col <- rlang::ensym(group_col)   # convert to symbol for tidy eval
  
  # count rows per group, randomize group order, accumulate, and select groups
  selected_groups <- df %>%
    group_by(!!group_col) %>%
    mutate(group_n = n()) %>%
    ungroup() %>%
    distinct(!!group_col, group_n) %>%
    sample_frac(1) %>%                       # randomize order
    mutate(cum_n = cumsum(group_n)) %>%
    filter(cum_n <= target_n |               # keep all groups before target
             lag(cum_n, default = 0) < target_n) %>%  # include first group > target
    pull(!!group_col)
  
  # return full rows from selected groups
  df %>% filter((!!group_col) %in% selected_groups)
}

repeat_sampling_bind <- function(df, group_col, target_n, reps = 1000, base_seed = NULL) {
  
  if (!is.null(base_seed)) set.seed(base_seed)
  
  dplyr::bind_rows(
    lapply(seq_len(reps), function(i) {
      iter_seed <- if (!is.null(base_seed)) base_seed + i else NULL
      
      sample_whole_traps(
        df = df,
        group_col = {{group_col}},
        target_n = target_n,
        seed = iter_seed
      ) %>% mutate(iteration = i)
    })
  )
}

randomtrap_samples_df <- repeat_sampling_bind(
  df = wb.25,
  group_col = Trap_Process_Order,
  target_n = 100,
  reps = 1000,
  base_seed = 123
)

iteration_summary <- randomtrap_samples_df %>%
  group_by(iteration) %>%
  summarise(
    mean_CW = mean(CW_mm, na.rm = TRUE),
    var_CW  = var(CW_mm,  na.rm = TRUE),
    max_CW  = max(CW_mm,  na.rm = TRUE),
    min_CW  = min(CW_mm,  na.rm = TRUE)
  )

full_summary <- wb.25 %>%
  summarise(
    full_mean_CW = mean(CW_mm, na.rm = TRUE),
    full_var_CW  = var(CW_mm,  na.rm = TRUE),
    full_max_CW  = max(CW_mm,  na.rm = TRUE),
    full_min_CW  = min(CW_mm,  na.rm = TRUE)
  )

iteration_with_diff <- iteration_summary %>%
  bind_cols(full_summary) %>%   # adds full-sample values
  mutate(
    diff_mean_CW = mean_CW - full_mean_CW,
    diff_var_CW  = var_CW  - full_var_CW,
    diff_max_CW  = max_CW  - full_max_CW,
    diff_min_CW  = min_CW  - full_min_CW
  )

####Plots 1

#Compute percent differences from full_summary
iteration_with_pctdiff <- iteration_summary %>%
  bind_cols(full_summary) %>%
  mutate(
    pctdiff_mean_CW = 100 * (mean_CW - full_mean_CW) / full_mean_CW,
    pctdiff_var_CW  = 100 * (var_CW  - full_var_CW)  / full_var_CW,
    pctdiff_max_CW  = 100 * (max_CW  - full_max_CW)  / full_max_CW,
    pctdiff_min_CW  = 100 * (min_CW  - full_min_CW)  / full_min_CW
  )

#Convert to long format for easier processing
iter_pct_long <- iteration_with_pctdiff %>%
  select(
    iteration,
    pctdiff_mean_CW,
    pctdiff_var_CW,
    pctdiff_max_CW,
    pctdiff_min_CW
  ) %>%
  tidyr::pivot_longer(
    cols = -iteration,
    names_to = "variable",
    values_to = "pct_diff"
  )

#Define continuous boundary ranges for percent difference
range_defs <- tibble::tribble(
  ~variable,           ~min_bound, ~max_bound,
  "pctdiff_mean_CW",        0,         15,
  "pctdiff_var_CW",         0,         50,
  "pctdiff_max_CW",         0,         20,
  "pctdiff_min_CW",         0,         200
)

#Create a continuous grid of boundary values
bound_grid <- range_defs %>%
  rowwise() %>%
  mutate(bound_seq = list(seq(min_bound, max_bound, length.out = 100))) %>%
  tidyr::unnest(bound_seq) %>%
  rename(bound = bound_seq)

#Combine percent-differences with boundaries
surface_df <- map_dfr(
  split(iter_pct_long, iter_pct_long$variable),
  function(df_var) {
    
    this_var <- unique(df_var$variable)
    
    # Filter correct bound row for this variable, then DROP the 'variable' column
    bounds_var <- bound_grid %>%
      filter(variable == this_var) %>%
      select(-variable)   # <<< IMPORTANT FIX
    
    # Cartesian product for this variable only
    tidyr::crossing(df_var, bounds_var) %>%
      group_by(variable, bound) %>%
      summarise(
        percent_outside = mean(abs(pct_diff) > bound) * 100,
        .groups = "drop"
      )
  }
)

make_surface_strip_plot <- function(df, varname) {
  
  # Ensure the data is sorted by boundary
  df <- df %>% arrange(bound)
  
  # Interpolate the boundary where percent_outside = 5%
  # approx() returns the y-value (boundary) corresponding to xout = 5
  if (any(df$percent_outside >= 5)) {
    line_bound <- approx(
      x = df$percent_outside,
      y = df$bound,
      xout = 5
    )$y
  } else {
    line_bound <- NA
  }
  
  ggplot(df, aes(
    y = bound,
    x = 1,
    fill = percent_outside
  )) +
    geom_tile(width = 0.5) +
    scale_fill_viridis_c(
      option = "plasma",
      limits = c(0, 100),
      oob = scales::squish
    ) +
    # Draw horizontal line at the boundary corresponding to 5% outside
    {if (!is.na(line_bound)) geom_hline(yintercept = line_bound, color = "blue", linetype = "dashed", linewidth = 1)} +
    theme_bw(base_size = 14) +
    labs(
      title = paste("Percent Outside vs Boundary:", varname),
      y = "Boundary (Percent Difference)",
      x = NULL,
      fill = "% Outside"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
}

variables <- unique(surface_df$variable)

plot_list <- lapply(
  variables,
  function(v) {
    df_var <- surface_df %>% filter(variable == v)
    make_surface_strip_plot(df_var, v)
  }
)

names(plot_list) <- variables

plot_list$pctdiff_mean_CW
plot_list$pctdiff_max_CW
plot_list$pctdiff_var_CW
plot_list$pctdiff_min_CW



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