# Libraries -----------------------------------
library(data.table)
library(pwr)
library(ggplot2)
library(tidyr)
library(ReplicationSuccess)
library(osfr)
library(effectsize)
library(dplyr)
library(ggrepel)
library(pROC)
library(PRROC)
library(readxl)

# Import empyrical effect sizes --------------------------------------------

data_dir <- "empyrical_effect_size_datasets"

#### Empyrical ES dataset: Bonapersona 2021 #### 
# (https://www.nature.com/articles/s41593-020-00792-3)
# field: neuroscience and metabolism
# login to OSF via token
# osf_auth("6IKuU0C6eBogBCNmilT4s972LnDZJrQlfggilZEbG7wxmal6Ooifewwb18A5SkgIwETEpd")
# 
# file <- osf_retrieve_file("8cqwa")
# 
# osf_download(file,
#             path = "C:/Users/paspe/Desktop/scripts/ResearchTrajectory-master/input",
#             conflicts = "overwrite")

dataset_name <- "Bonapersona_2021"

meta <- read.csv(
  file.path(data_dir, "meta_effectsize_bonapersona.csv")
)

save_dir <- file.path("simulation_results", dataset_name)
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

save_path <- file.path(save_dir, paste0("effect_sizes_", dataset_name, ".pdf"))

# #### Empyrical ES dataset: Carneiro 2018 (Pessimistic) ####
# # # (doi: 10.1371/journal.pone.0196258.)
# # # Note: this df has smaller and negative ES (S error). Outcome: rodent fear conditioning.
# 
# dataset_name <- "Carneiro_2018"
# 
# load(file.path(data_dir, "es_data_carneiro.RData"))
# meta <- ES_data_Carneiro
# data.table::setnames(meta, "ES_d", "yi")
# 
# save_dir <- file.path("simulation_results", dataset_name)
# if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
# 
# save_path <- file.path(save_dir, paste0("effect_sizes_", dataset_name, ".pdf"))
# 
# #### Empyrical ES dataset: Rosso 2022 ####
# # Anxiety dataset, https://doi.org/10.1016/j.neubiorev.2022.104928
# dataset_name <- "Rosso_2022"
# 
# meta <- readxl::read_excel(
#   file.path(data_dir, "SR1_0_forR_rosso.xlsx"),
#   sheet = "default"
# )
# 
# save_dir <- file.path("simulation_results", dataset_name)
# if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
# 
# save_path <- file.path(save_dir, paste0("effect_sizes_", dataset_name, ".pdf"))

#### ES distribution plot ####

effect_sizes <- 
  ggplot(meta, aes(x = yi)) +
  geom_histogram(
    bins = 60,
    color = "white",
    linewidth = 0.3,
    fill = "grey40"
  ) +
  labs(
    x = expression(paste("Hedge's ", italic("g"), " (log scale)")),
    y = "Number of effect sizes",
    title = NULL
  ) +
  scale_x_continuous(
    trans = "pseudo_log",
    breaks = c(-20, -7, -3, -1.5, -0.8, 0, 0.8, 1.5, 3, 7, 20),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_minimal(base_size = 18, base_family = "Helvetica") +
  theme(
    plot.title = element_text(
      size = 22,
      face = "bold",
      hjust = 0.5
    ),
    axis.title.x = element_text(
      size = 20,
      margin = margin(t = 18)
    ),
    axis.title.y = element_text(
      size = 20,
      margin = margin(r = 12)
    ),
    axis.text = element_text(size = 10, color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
    axis.ticks.x = element_line(color = "black", linewidth = 0.4),
    axis.ticks.length.x = unit(4, "pt"),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
  )

effect_sizes

plot_name <- paste0("Effect_size_", dataset_name, ".png")
ggsave(file.path(save_path, plot_name), 
       effect_sizes, width = 12, height = 8, dpi = 150)



# Functions for study simulation ------------------------------------------


# Function to simulate a single study and return all statistics
simulate_study_complete <- function(true_effect, n_per_group, study_type = "exploratory") {
  # Generate data
  control_group <- rnorm(n_per_group, mean = 0, sd = 1)
  treatment_group <- rnorm(n_per_group, mean = true_effect, sd = 1)
  
  # Perform t-test
  t_result <- t.test(treatment_group, control_group, var.equal = TRUE)
  
  # Calculate effect sizes
  pooled_sd <- sqrt(((n_per_group - 1) * var(control_group) + 
                       (n_per_group - 1) * var(treatment_group)) / 
                      (2 * n_per_group - 2))
  
  cohens_d <- (mean(treatment_group) - mean(control_group)) / pooled_sd
  
  # Hedges' g correction
  j <- 1 - (3 / (4 * (2 * n_per_group) - 9))
  hedges_g <- cohens_d * j
  
  # Standard error for Hedges' g
  se_g <- sqrt((2 * n_per_group) / (n_per_group^2) + (hedges_g^2) / (4 * n_per_group))
  
  return(list(
    observed_g = hedges_g,
    se_g = se_g,
    p_value = t_result$p.value,
    t_statistic = t_result$statistic,
    ci_lower = hedges_g - 1.96 * se_g,
    ci_upper = hedges_g + 1.96 * se_g,
    sample_size = n_per_group
  ))
}

# Full research trajectory
run_research_trajectory <- function(true_effect, exploratory_n, shrinkage_factor) {
  
  # Run exploratory study and filter statistically significant ones
  exploratory_results <- simulate_study_complete(true_effect, exploratory_n, "exploratory")
  
  if (exploratory_results$p_value < alpha) {
    sig_exploratory_results <- exploratory_results
  } else {
    # Return early for non-significant exploratory results
    return(list(
      exploratory = exploratory_results,
      confirmatory = list(observed_g = NA, se_g = NA, p_value = NA, 
                          ci_lower = NA, ci_upper = NA, sample_size = NA),
      confirmatory_ss = NA,
      true_confirmatory_effect = NA
    ))
  }
  
  # Calculate confirmatory sample size based on OBSERVED effect
  observed_effect <- sig_exploratory_results$observed_g
  
  # Edge case: observed effect is extremely small, resulting in very high sample size
  if (abs(observed_effect) < 0.05) {
    return(list(
      exploratory = sig_exploratory_results,
      confirmatory = list(observed_g = NA, se_g = NA, p_value = NA, 
                          ci_lower = NA, ci_upper = NA, sample_size = NA),
      confirmatory_ss = NA,
      true_confirmatory_effect = NA
    ))
  }
  
  # Calculate required sample size for 80% power
  confirmatory_n <- tryCatch({
    ceiling(pwr.t.test(d = abs(observed_effect), power = 0.8, 
                       sig.level = alpha, type = "two.sample")$n)
  }, error = function(e) NA)
  
  # Cap sample size at realistic maximum (50)
  if (!is.na(confirmatory_n) && confirmatory_n > 50) {
    confirmatory_n <- 50
  }
  
  if (is.na(confirmatory_n) || confirmatory_n < 5) {
    return(list(
      exploratory = sig_exploratory_results,
      confirmatory = list(observed_g = NA, se_g = NA, p_value = NA, 
                          ci_lower = NA, ci_upper = NA, sample_size = NA),
      confirmatory_ss = confirmatory_n,
      true_confirmatory_effect = NA
    ))
  }
  
  # Calculate true effect for confirmatory study by applying the shrinakge factor:
  # Note: shrinkage factors won't be modeled here as thy can have many sources like winner's curse, publication bias, flexible experimental methods/analyses.
  true_confirmatory_effect <- true_effect * (1 - shrinkage_factor)
  
  # Run confirmatory study
  confirmatory_results <- simulate_study_complete(true_confirmatory_effect, 
                                                  confirmatory_n, "confirmatory")
  
  return(list(
    exploratory = sig_exploratory_results,
    confirmatory = confirmatory_results,
    confirmatory_ss = confirmatory_n,
    true_confirmatory_effect = true_confirmatory_effect
  ))
}

# Function to extract results into a clean data.table with clear naming
extract_results <- function(results_list) {
  dt_list <- lapply(results_list, function(x) {
    data.table(
      sim_id = x$sim_id,
      exploratory_n = x$exploratory_n,
      shrinkage = x$shrinkage,
      true_effect = x$true_effect,
      
      # Exploratory results (clear naming structure)
      exploratory_p = x$exploratory$p_value,
      exploratory_ss = x$exploratory$sample_size,
      exploratory_ci_lower = x$exploratory$ci_lower,
      exploratory_ci_upper = x$exploratory$ci_upper,
      exploratory_true_g = x$true_effect,  # True effect used for exploratory
      exploratory_observed_g = x$exploratory$observed_g,
      exploratory_se_g = x$exploratory$se_g,
      
      # Confirmatory results (clear naming structure)
      confirmatory_p = x$confirmatory$p_value,
      confirmatory_ss = x$confirmatory$sample_size,
      confirmatory_ci_lower = x$confirmatory$ci_lower,
      confirmatory_ci_upper = x$confirmatory$ci_upper,
      confirmatory_true_g = x$true_confirmatory_effect,  # True effect used for confirmatory
      confirmatory_observed_g = x$confirmatory$observed_g,
      confirmatory_se_g = x$confirmatory$se_g
    )
  })
  
  rbindlist(dt_list)
}

# Simulation parameters ---------------------------------------------------

set.seed(42)
n_simulations <- 10000
exploratory_sample_sizes <- c(5, 10, 15, 20)
alpha <- 0.05
exp_ss_heatmap <- 10  # select ss to plot in success rate heatmap
shrinkage_levels <- c(0.0, 0.2, 0.5, 0.8, 0.99)
SESOI_g <- 0.4

# Use realistic effect size distribution
true_effects <- sample(meta$yi[abs(meta$yi) >= 0.1 & abs(meta$yi) <= 2], 
                       n_simulations, replace = TRUE)

# Create simulation grid
simulation_grid <- expand.grid(
  sim_id = 1:n_simulations,
  exploratory_n = exploratory_sample_sizes,
  shrinkage = shrinkage_levels,
  stringsAsFactors = FALSE
)

# Add true effects
simulation_grid$true_effect <- rep(true_effects, length.out = nrow(simulation_grid))


# Simulation START --------------------------------------------------------

# Run simulations
results_list <- vector("list", nrow(simulation_grid))

for (i in seq_len(nrow(simulation_grid))) {
  if (i %% 100 == 0) cat("Processing simulation", i, "of", nrow(simulation_grid), "\n")
  
  # For each row, run the research trajectory
  row <- simulation_grid[i, ]
  results_list[[i]] <- run_research_trajectory(
    true_effect = row$true_effect,
    exploratory_n = row$exploratory_n,
    shrinkage_factor = row$shrinkage
  )
  
  # Add metadata
  results_list[[i]]$sim_id <- row$sim_id
  results_list[[i]]$exploratory_n <- row$exploratory_n
  results_list[[i]]$shrinkage <- row$shrinkage
  results_list[[i]]$true_effect <- row$true_effect
}

# Extract results
simulation_results <- extract_results(results_list)

# Filter only plausible effect sizes (< 2 and not NA)
simulation_results <- simulation_results[
  abs(exploratory_observed_g) <= 2 & abs(confirmatory_observed_g) <= 2
]

# Effect size classification
simulation_results <- simulation_results[,
                                         effect_class := fifelse(abs(true_effect) < 0.2, "Negligible",
                                                                 fifelse(abs(true_effect) < 0.5, "Small", 
                                                                         fifelse(abs(true_effect) < 0.8, "Medium", "Large")))
]

# Replication methods -----------------------------------------------------


# Start with basic criteria of replication success

# METHOD 1: Direction consistency
simulation_results[, direction_consistent := sign(exploratory_observed_g) == sign(confirmatory_observed_g)]

# Add more methods that build on direction consistency
simulation_results[, `:=`(
  # METHOD 2: CI overlap  and same sign
  ci_overlap = !(exploratory_ci_upper < confirmatory_ci_lower | 
                   confirmatory_ci_upper < exploratory_ci_lower) & direction_consistent,
  
  # METHOD 3: Exploratory effect in confirmatory CI  and same sign
  exp_in_conf_ci = exploratory_observed_g >= confirmatory_ci_lower & 
    exploratory_observed_g <= confirmatory_ci_upper & direction_consistent,
  
  # METHOD 4: two trial rule, both exploratory and confirmatory significant and same sign
  both_significant = exploratory_p < alpha & confirmatory_p < alpha & direction_consistent
)]

# METHOD 5: Sceptical p-value.
# doi: 10.1111/rssa.12493
# From the original study, it establish a sceptical prior centered at 0 (no effect), the posterior has a CI that has to touch the lower 
# bound with 0 so through a reverse Bayes we back calculate from the original study and this posterior (blue) a sufficiently sceptical prior 
# (meaning that the CI has to be the result of combining the original study and the posterior fixed at 0). We use this sufficiently sceptical
# prior to assess the replication study. The larger the effect of the original study the narrower will be the CI of the sceptical prior that 
# has to counterbalance the claim of the original study. 

simulation_results[, sceptical_p := ifelse(                                   # !is.finite() captures NA, NaN, and Inf.
  !is.finite(exploratory_observed_g) | !is.finite(exploratory_se_g) | 
    !is.finite(confirmatory_observed_g) | !is.finite(confirmatory_se_g) |
    exploratory_se_g == 0 | confirmatory_se_g == 0,
  NA,
  {
    z_o <- exploratory_observed_g / exploratory_se_g        # z values (g / se) of the exploratory and confirmatory studies
    z_r <- confirmatory_observed_g / confirmatory_se_g
    var_ratio <- (exploratory_se_g / confirmatory_se_g)^2   # ratio of the standard errors, larger c means the confirmation has higher precision
    
    # pSceptical fx. Type is the calibration of the p-value. It can be "golden", "nominal" or "controlled" depending on how conservative we want the threshold for replication success to be.
    # Golden: success only when replication ES is at least as large as the original
    # Controlled: controls the overall type I error
    # Nominal: no recalibration, the p-value is interpreted as is
    p_s <- pSceptical(z_o, z_r, c = var_ratio, alternative = "two.sided", type = "golden") 
    p_s
  }
), by = sim_id]
simulation_results[, sceptical_significant := sceptical_p < alpha]


## METHOD 6: Small Telescopes. 
## The logic is: Small telescope = small power. If the confirmatory has d below d33 of the original study than the original study wouldn't have 
# enough power (small telescope) to detect such small effect.
# The exploratory studies are usually underpowered, the effect size that corresponds to 33% power d33 will be relatively large so it will be hard 
# for the confirmatory effect size to not be smaller than d33

# compute the d33
simulation_results[, d33 := mapply(function(n) {
  if (is.na(n)) return(NA_real_)
  pwr.t.test(n = n, sig.level = alpha, power = 0.33,
             type = "two.sample", alternative = "two.sided")$d       # extract d33 from the exploratory sample sizes
}, exploratory_ss)]

# get the z score and compute the lower ci
z_crit <- qnorm(1 - alpha)
simulation_results[, lower_ci := confirmatory_observed_g - z_crit * confirmatory_se_g]

# test superiority of the confirmatory d to the exploratory d33
simulation_results[, small_telescope := direction_consistent & 
                     (
                       (exploratory_observed_g > 0 & confirmatory_ci_lower >  d33) |
                         (exploratory_observed_g < 0 & confirmatory_ci_upper < -d33)
                     )
]

## METHOD 7: Smallest Detectable Effect (SDE)
## Test if the confirmatory g is superior than the smallest detectable effect size of the exploratory at 80% power

# compute the SDE first
simulation_results[, sde := mapply(function(n) {
  if (is.na(n)) return(NA_real_)
  pwr.t.test(n = n, sig.level = alpha, power = 0.8,
             type = "two.sample", alternative = "two.sided")$d
}, confirmatory_ss)]


# Apply Hedges' correction to SDE
simulation_results[, sde_g := sde * (1 - (3 / (4 * (2 * confirmatory_ss) - 9)))]

# Standard error for SDE effect size (equal group sizes)
simulation_results[, se_sde := sqrt((2 * confirmatory_ss)/(confirmatory_ss^2) + (sde_g^2 / (4 * confirmatory_ss)))]

# Check if confirmatory ES exceeds SDE threshold
simulation_results[, sde_confirmed := direction_consistent & 
                     (
                       (exploratory_observed_g > 0 & confirmatory_ci_lower >  sde_g) |
                         (exploratory_observed_g < 0 & confirmatory_ci_upper < -sde_g)
                     )
]


## METHOD 8: Smallest Effect Size of Interest (SESOI)
# Test whether the effect size of the confirmatory is above a pre-specified (i.e. clinically meaningful) threshold

# criteria 1: the SESOI has to fall in the CI of the confirmatory g
# simulation_results[, sesoi_in_ci := (SESOI_g >= confirmatory_ci_lower & SESOI_g <= confirmatory_ci_upper) |
#                      (-SESOI_g >= confirmatory_ci_lower & -SESOI_g <= confirmatory_ci_upper)]

# criteria 2: the confirmatory CI has to be above the the SESOI
simulation_results[, ci_above_sesoi := confirmatory_ci_lower >= SESOI_g | confirmatory_ci_upper <= -SESOI_g]



# Simulation results ------------------------------------------------------

simulation_results[, expl_effect_size_class := factor(effect_class,
                                                      levels = c("Negligible", "Small", "Medium", "Large")
)]


replication_criteria <- c("direction_consistent", 
                          "ci_overlap", 
                          "exp_in_conf_ci",
                          "both_significant", 
                          "sceptical_significant", 
                          "small_telescope", 
                          "sde_confirmed", 
                          "ci_above_sesoi")

# Clean criterion names
criterion_names <- c(
  "direction_consistent" = "Direction Consistent",
  "ci_overlap" = "CI Overlap", 
  "exp_in_conf_ci" = "Exp in Conf CI",
  "both_significant" = "Both Significant",
  "sceptical_significant" = "Sceptical p-value", 
  "small_telescope" = "Small Telescope",
  "sde_confirmed" = "SDE Confirmed",
  "ci_above_sesoi" = "Conf CI Above SESOI"
)

# Summarise by shrinkage and effect class
summary_data_all <- simulation_results[,
                                   lapply(.SD, function(x) mean(as.logical(x), na.rm = TRUE)),
                                   by = .(shrinkage, expl_effect_size_class, exploratory_n),
                                   .SDcols = replication_criteria
]

# Filter for the exploratory ss of interest
summary_data <- summary_data_all[exploratory_n == exp_ss_heatmap]


# Prepare data for plotting
plot_data <- melt(summary_data,
                  id.vars = c("shrinkage", "expl_effect_size_class", "exploratory_n"),
                  variable.name = "Criterion",
                  value.name = "SuccessRate"
)
plot_data <- setDT(plot_data)

plot_data[, g_shrinkage := paste0(round(shrinkage * 100), "%")]

# Remove negligible as it's likely not replicated even if significant
plot_data <- plot_data[expl_effect_size_class != "Negligible"]  

plot_data[, Criterion_clean := factor(
  criterion_names[Criterion],
  levels = criterion_names[replication_criteria]
)]

#### Heatmap of the success rate ####
combined_heatmap <- ggplot(plot_data, aes(x = g_shrinkage, y = Criterion_clean, fill = SuccessRate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", SuccessRate)), size = 3.0) +
  scale_fill_gradient(low = "white", high = "#0072B2", limits = c(0, 1), name = "Success Rate") +
  facet_wrap(~ expl_effect_size_class, ncol = 3) +
  labs(
    title = "Replication Success Rate by Effect Size Class",
    x = "Replication Criterion",
    y = "Effect size shrinkage"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 15, face = "bold")
  )

print(combined_heatmap)

# Save the heatmap for the pre-sepcified exploratory sample size
heatmap_name <- paste0("combined_heatmap_success_rate_all_classes_expn", exp_ss_heatmap, ".png")
ggsave(file.path(save_path, heatmap_name), 
       combined_heatmap, width = 12, height = 8, dpi = 150)

#### Type I error ####
# Calculate False Positive Rate (FPR) for 99% shrinkage cases
# FPR = rate of claiming replication success when true effect is around 0
fpr_data <- simulation_results[shrinkage == 0.99,
                               lapply(.SD, function(x) mean(as.logical(x), na.rm = TRUE)),
                               by = .(exploratory_n),
                               .SDcols = replication_criteria]

# Prepare FPR data for plotting
fpr_plot_data <- melt(fpr_data,
                      id.vars = "exploratory_n",
                      variable.name = "Criterion",
                      value.name = "FPR"
)
fpr_plot_data <- setDT(fpr_plot_data)


# FPR bar plot per method (pooling exploratory sample sizes)
pr_pooled_data <- fpr_plot_data %>%
  group_by(Criterion) %>%
  summarise(mean_FPR = mean(FPR, na.rm = TRUE), .groups = 'drop')

fpr_plot_pooled_n <- ggplot(pr_pooled_data, aes(x = mean_FPR, y = Criterion)) +
  geom_col(alpha = 0.8, width = 0.7) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  labs(
    title = "False Positive Rate by Replication Criterion",
    subtitle = "Averaged across exploratory sample sizes (99% shrinkage scenarios)",
    x = "Replication Criterion",
    y = "False Positive Rate",
    fill = NULL
  ) +
  theme_minimal(base_size = 24) +
  theme(
    legend.position = "none",  # Remove legend since x-axis labels show criteria
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 15, color = "gray60")
  )

print(fpr_plot_pooled_n)


# Additional analysis: FPR by effect size class
fpr_by_class <- simulation_results[shrinkage == 0.99,
                                   lapply(.SD, function(x) mean(as.logical(x), na.rm = TRUE)),
                                   by = .(exploratory_n, expl_effect_size_class),
                                   .SDcols = replication_criteria]
print(fpr_by_class, digits = 3)

# Create comprehensive FPR heatmap
fpr_heatmap_data <- melt(fpr_by_class,
                         id.vars = c("exploratory_n", "expl_effect_size_class"),
                         variable.name = "Criterion",
                         value.name = "FPR")

fpr_heatmap <- ggplot(fpr_heatmap_data, aes(x = Criterion, y = factor(exploratory_n), fill = FPR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", FPR)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "#D55E00", limits = c(0, 1), 
                      name = "False Positive Rate") +
  facet_wrap(~ expl_effect_size_class, ncol = 2) +
  labs(
    title = "False Positive Rate Heatmap by Effect Size Class",
    subtitle = "Based on 99% shrinkage scenarios (true confirmatory effect ≈ 0)",
    x = "Replication Criterion",
    y = "Exploratory Sample Size (n per group)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray60"),
    strip.text = element_text(size = 11, face = "bold")
  )

print(fpr_heatmap)

# Save FPR heatmap
ggsave(file.path(save_path, "fpr_heatmap_by_class.png"), fpr_heatmap, 
       width = 12, height = 8, dpi = 300)


#### Shrinkage sensitivity ####
shrinkage_sensitivity <- simulation_results[, {
  sensitivity_results <- list()
  
  for(crit in replication_criteria) {
    crit_values <- as.numeric(get(crit))
    
    # Remove NA values for correlation
    valid_indices <- !is.na(crit_values) & !is.na(shrinkage)
    
    if(sum(valid_indices) > 10) {  # Need sufficient data points
      correlation <- cor(shrinkage[valid_indices], crit_values[valid_indices])
      
      # Linear model to get slope
      lm_model <- lm(crit_values[valid_indices] ~ shrinkage[valid_indices])
      slope <- coef(lm_model)[2]
      
      sensitivity_results[[crit]] <- data.table(
        criterion = crit,
        correlation = correlation,
        slope = slope
      )
    }
  }
  
  rbindlist(sensitivity_results)
}, by = .(exploratory_ss, expl_effect_size_class)]

# Remove negligible
shrinkage_sensitivity <- shrinkage_sensitivity[expl_effect_size_class != "Negligible"]

# Average across conditions
overall_sensitivity <- shrinkage_sensitivity[, .(
  mean_correlation = mean(correlation, na.rm = TRUE),
  mean_slope = mean(slope, na.rm = TRUE),
  se_correlation = sd(correlation, na.rm = TRUE) / sqrt(.N)
), by = criterion]

overall_sensitivity[, criterion_clean := criterion_names[criterion]]
overall_sensitivity <- overall_sensitivity[order(mean_correlation)]  # More negative = more sensitive

print(overall_sensitivity[, .(criterion_clean, mean_correlation, se_correlation)], digits = 3)

# Add clean criterion names
shrinkage_sensitivity[, criterion_clean := criterion_names[criterion]]

# Plot correlation between shrinkage and success rates
sensitivity_plot <- ggplot(shrinkage_sensitivity, 
                           aes(
                             x = correlation, 
                             y = reorder(criterion_clean, correlation),
                             fill = expl_effect_size_class)
                           ) +
  geom_boxplot(width = 0.6) +   
  facet_wrap(~expl_effect_size_class, ncol = 3) +
  scale_fill_brewer(palette = "Set4", type = "qual", name = "Exploratory Effect Size") +
  labs(
    title = "Shrinkage Sensitivity by Effect Size Class",
    x = "Correlation (Shrinkage ~ Replication Success)",
    y = "Replication Criterion"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold")
  )

print(sensitivity_plot)

ggsave(file.path(save_path, "shrinkage_sensitivity.png"), sensitivity_plot, 
       width = 12, height = 8, dpi = 300)




#### Precision-Recall ####

# Define ground truth: 
# True Positive = 0% or 20% shrinkage (genuine replications)
# True Negative = 99% shrinkage (failed replications)
pr_data <- simulation_results[shrinkage  %in% c(0, 0.2, 0.99)]

# Create binary ground truth labels
pr_data[, ground_truth := ifelse(shrinkage %in% c(0, 0.2), 1, 0)]

# Calculate precision-recall for each criterion and condition
pr_results <- list()

for(effect_class in unique(pr_data$expl_effect_size_class)) {
  for(sample_size in unique(pr_data$exploratory_n)) {
    
    subset_data <- pr_data[expl_effect_size_class == effect_class & 
                             exploratory_n == sample_size]
    
    if(nrow(subset_data) == 0) next
    
    for(crit in replication_criteria) {
      
      # Get predictions and ground truth
      predictions <- as.numeric(subset_data[[crit]])
      ground_truth <- subset_data$ground_truth
      
      # Remove NA values
      valid_idx <- !is.na(predictions) & !is.na(ground_truth)
      if(sum(valid_idx) < 10) next  # Need sufficient data
      
      predictions <- predictions[valid_idx]
      ground_truth <- ground_truth[valid_idx]
      
      # Calculate precision-recall curve
      if(length(unique(ground_truth)) > 1 && length(unique(predictions)) > 1) {
        
        # For binary predictions, calculate precision and recall directly
        if(length(unique(predictions)) == 2) {
          
          # Calculate confusion matrix components
          tp <- sum(predictions == 1 & ground_truth == 1)
          fp <- sum(predictions == 1 & ground_truth == 0)
          fn <- sum(predictions == 0 & ground_truth == 1)
          tn <- sum(predictions == 0 & ground_truth == 0)
          
          precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
          recall <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
          f1_score <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
          
          # Also calculate AUC-PR using PRROC if possible
          if(var(predictions) > 0) {
            pr_curve <- PRROC::pr.curve(scores.class0 = predictions[ground_truth == 1],
                                        scores.class1 = predictions[ground_truth == 0],
                                        curve = TRUE)
            auc_pr <- pr_curve$auc.integral
          } else {
            auc_pr <- NA
          }
          
        } else {
          # For continuous predictions, use PRROC
          pr_curve <- PRROC::pr.curve(scores.class0 = predictions[ground_truth == 1],
                                      scores.class1 = predictions[ground_truth == 0],
                                      curve = TRUE)
          
          precision <- mean(pr_curve$curve[,2], na.rm = TRUE)
          recall <- mean(pr_curve$curve[,1], na.rm = TRUE)
          f1_score <- ifelse(precision + recall > 0, 2 * precision * recall / (precision + recall), 0)
          auc_pr <- pr_curve$auc.integral
        }
        
        pr_results[[length(pr_results) + 1]] <- data.table(
          criterion = crit,
          effect_class = effect_class,
          sample_size = sample_size,
          precision = precision,
          recall = recall,
          f1_score = f1_score,
          auc_pr = auc_pr,
          n_obs = length(predictions),
          n_positive = sum(ground_truth == 1),
          n_negative = sum(ground_truth == 0)
        )
      }
    }
  }
}

# Combine results
pr_summary <- rbindlist(pr_results)
pr_summary[, criterion_clean := criterion_names[criterion]]

# Calculate overall performance across conditions
overall_pr <- pr_summary[, .(
  mean_precision = mean(precision, na.rm = TRUE),
  mean_recall = mean(recall, na.rm = TRUE),
  mean_f1 = mean(f1_score, na.rm = TRUE),
  mean_auc_pr = mean(auc_pr, na.rm = TRUE),
  se_precision = sd(precision, na.rm = TRUE) / sqrt(.N),
  se_recall = sd(recall, na.rm = TRUE) / sqrt(.N),
  se_f1 = sd(f1_score, na.rm = TRUE) / sqrt(.N)
), by = .(criterion, criterion_clean)]

# Sort by F1 score
overall_pr <- overall_pr[order(-mean_f1)]

print("(True Positives: 0% or 20% shrinkage; True Negatives: 99% shrinkage)")
print(overall_pr[, .(criterion_clean, mean_precision, mean_recall, mean_f1, mean_auc_pr)], digits = 3)


# Plot F1 Score comparison
f1_plot <- ggplot(overall_pr, aes(x = reorder(criterion_clean, mean_f1), y = mean_f1)) +
  geom_col(fill = "#0066B5", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean_f1 - se_f1, ymax = mean_f1 + se_f1), 
                width = 0.2, alpha = 0.7) +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = "F1 Score by Replication Criterion",
    subtitle = "Higher F1 indicates better balance of precision and recall",
    x = "Replication Criterion",
    y = "F1 Score",
    caption = "Error bars show standard error across effect size classes and sample sizes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "gray60")
  )

print(f1_plot)

ggsave(file.path(save_path, "f1_score_comparison.png"), f1_plot, 
       width = 10, height = 8, dpi = 300)

# Precision vs Recall scatter plot
pr_scatter <- ggplot(overall_pr, aes(x = mean_recall, y = mean_precision)) +
  geom_point(size = 4, alpha = 0.8, color = "#0066B5") +
  geom_text_repel(aes(label = criterion_clean), size = 4, max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = "Precision vs Recall",
    x = "Recall (Sensitivity)",
    y = "Precision (Positive Predictive Value)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 15)),
    plot.subtitle = element_text(size = 15, color = "gray60")
  )

print(pr_scatter)

ggsave(file.path(save_path, "precision_recall_scatter.png"), pr_scatter, 
       width = 10, height = 8, dpi = 300)



#### Summary table ####

summary_table <- overall_pr[, .(
  Criterion = criterion_clean,
  `Precision (%)` = round(mean_precision * 100, 1),
  `Recall (%)` = round(mean_recall * 100, 1),
  `F1 Score (%)` = round(mean_f1 * 100, 1),
  `AUC-PR (%)` = round(mean_auc_pr * 100, 1)
)][order(-`AUC-PR (%)`)]

print(summary_table)

# Best performing criterion
best_criterion <- summary_table[1, Criterion]
cat(sprintf("\nBest performing criterion: %s\n", best_criterion))
cat(sprintf("F1 Score: %.1f%%\n", summary_table[1, `F1 Score (%)`]))
cat(sprintf("Precision: %.1f%%, Recall: %.1f%%\n", 
            summary_table[1, `Precision (%)`], 
            summary_table[1, `Recall (%)`]))

