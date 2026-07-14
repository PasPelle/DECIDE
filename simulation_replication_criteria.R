# Libraries and Fx -----------------------------------
source(here("packages.R"))
source(here("all_functions.R"))

# Import empyrical effect sizes --------------------------------------------
# data_dir <- "empyrical_effect_size_datasets"
panels_dir <- file.path(save_dir_decide, "panels")
if (!dir.exists(panels_dir)) dir.create(panels_dir, recursive = TRUE)


#### Empyrical ES dataset: Bonapersona 2021 #### 
# (https://www.nature.com/articles/s41593-020-00792-3)
# field: neuroscience and metabolism
# login to OSF via token
# osf_auth("6IKuU0C6eBogBCNmilT4s972LnDZJrQlfggilZEbG7wxmal6Ooifewwb18A5SkgIwETEpd")
# #
# file <- osf_retrieve_file("8cqwa")
# 
# osf_download(file,
#             path = here("data"),
#             conflicts = "overwrite")

dataset_name <- "bonapersona_2021"

data_dir <- here("data", "empyrical_effect_size_datasets")

meta <- read.csv(
  file.path(data_dir, "meta_effectsize_bonapersona.csv")
)
# 
# save_dir <- file.path("simulation_results", dataset_name)
# if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

#### Empyrical ES dataset: Carneiro 2018 (Pessimistic) ####
# # (doi: 10.1371/journal.pone.0196258.)
# # Note: this df has smaller and negative ES (S error). Outcome: rodent fear conditioning.
# dataset_name <- "carneiro_2018"
# 
# load(here("data", "empyrical_effect_size_datasets", "es_data_carneiro.RData"))
# 
# setnames(ES_data_Carneiro, "ES_d", "yi")
# 
# meta <- ES_data_Carneiro
# save_dir <- here("simulation_results", dataset_name)
# if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# #### Empyrical ES dataset: Rosso 2022 ####
# Anxiety dataset, https://doi.org/10.1016/j.neubiorev.2022.104928
# dataset_name <- "rosso_2022"
# 
# meta <- read_excel(
#   file.path(data_dir, "SR1_0_forR_rosso.xlsx"),
#   sheet = "default"
# )
# # 
# save_dir <- file.path("simulation_results", dataset_name)
# if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)


# Global parameters
alpha <- 0.05
exp_ss_heatmap <- 10  # select ss for specific heatmap plots
k_labs <- 3


# SIMULATION --------------------------------------------------------------
RUN_SIMULATION <- FALSE  # set to TRUE to re-run

if (RUN_SIMULATION) {
  # simulation loop
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
    theme_minimal(base_size = 18, base_family = "Arial") +
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
  ggsave(file.path(save_dir, plot_name), 
         effect_sizes, width = 12, height = 8, dpi = 150)
  saveRDS(effect_sizes,
          here("results", "panels", paste0("effect_sizes_", dataset_name, ".rds")))
  

  # Simulation parameters ---------------------------------------------------
  
  set.seed(42)
  n_simulations <- 50000 # set to 50k for final Figure
  exploratory_sample_sizes <- c(5, 10, 15, 20)
  shrinkage_levels <- c(0.0, 0.2, 0.5, 0.8, 0.99)
  # SESOI_g <- 0.4
  
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
  
  results_list <- vector("list", nrow(simulation_grid))
  
  for (i in seq_len(nrow(simulation_grid))) {
    if (i %% 100 == 0) cat("Processing simulation", i, "of", nrow(simulation_grid), "\n") # print every 100 iteration
    
    # For each row, run the research trajectory
    row <- simulation_grid[i, ]
    results_list[[i]] <- run_research_trajectory_multilab(
      true_effect = row$true_effect,
      exploratory_n = row$exploratory_n,
      shrinkage_factor = row$shrinkage,
      k_labs = k_labs,                                  # comment out for single lab simulation
      # multilab_scale_factor = 2.8                     
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
  simulation_results[, .N, by = .(effect_class, exploratory_n, shrinkage)][order(N)]
  
  # Save raw simulation results
  saveRDS(simulation_results,
          here("results", paste0("simulation_results_raw_multilab_", dataset_name, ".rds")))
  
} else {
  simulation_results <- readRDS(
    file.path(save_dir, paste0("simulation_results_raw_multilab", dataset_name, ".rds")))
}

#### Count surviving simulations after exploratory p < 0.05 filter ####
# (these are already filtered in simulation_results since non-sig ones have NA confirmatory)
filter_counts <- simulation_results[effect_class %in% c("Small", "Medium", "Large"),
                                    .N,
                                    by = .(effect_class, exploratory_n)]

filter_counts[, effect_class := factor(effect_class, 
                                       levels = c("Small", "Medium", "Large"))]

n_conf_experiments_per_effect_size <- ggplot(filter_counts, aes(x = factor(exploratory_n), y = N, fill = effect_class)) +
  geom_col(position = "dodge") +
  # annotate("text", x = 0.6, y = 520, label = "N = 500 threshold", 
  #          color = "red", size = 3.5, hjust = 0) +
  scale_fill_brewer(palette = "Paired", name = "Effect size class") +
  labs(
    x = "Exploratory sample size (n per group)",
    y = "Simulations surviving exploratory filter",
    title = NULL
  ) +
  theme_prism(base_size = 14) +
  theme(legend.position = "right")
n_conf_experiments_per_effect_size

# ggsave(file.path(save_dir, "n_conf_experiments_per_effect_size.png"), 
#        n_conf_experiments_per_effect_size, 
#        width = 12, height = 8, dpi = 300)


# Run replication methods fx ----------------------------------------------

# Create the columns expected by compute_replication_flags()

simulation_results[, `:=`(
  id = sim_id,
  
  # exploratory
  g_exp     = exploratory_observed_g,
  se_exp    = exploratory_se_g,
  ci_lo_exp = exploratory_ci_lower,
  ci_hi_exp = exploratory_ci_upper,
  p_exp     = exploratory_p,
  n1_exp    = exploratory_ss,
  n2_exp    = exploratory_ss,
  
  # confirmatory
  g_conf     = confirmatory_observed_g,
  se_conf    = confirmatory_se_g,
  ci_lo_conf = confirmatory_ci_lower,
  ci_hi_conf = confirmatory_ci_upper,
  p_conf     = confirmatory_p
)]

# Provide SDE_pooled (required by de SDE rule) single-lab simulation
# simulation_results[, SDE_pooled := compute_sde_g(
#   n_per_group = confirmatory_ss,
#   alpha = alpha,
#   power = 0.8
# )]

# Provide SDE_pooled (required by de SDE rule) multi-lab simulation
simulation_results[, SDE_pooled := compute_sde_g(
  n_per_group = confirmatory_ss / k_labs,
  alpha = alpha,
  power = 0.8
)]

# Optional: SESOI
# simulation_results[, SESOI_g := 0.4]

# Run unified replication criteria computation
simulation_results <- compute_replication_flags(
  dt = simulation_results,
  cfg = list(
    alpha = alpha,
    z_pi = 1.96,
    power_d33 = 0.33,
    power_sde = 0.8,
    sceptical_type = "golden",
    sceptical_alternative = "two.sided"
  ),
  cols = list(
    id = "id",
    g_exp = "g_exp", se_exp = "se_exp", ci_lo_exp = "ci_lo_exp", ci_hi_exp = "ci_hi_exp", p_exp = "p_exp",
    n1_exp = "n1_exp", n2_exp = "n2_exp",
    g_conf = "g_conf", se_conf = "se_conf", ci_lo_conf = "ci_lo_conf", ci_hi_conf = "ci_hi_conf", p_conf = "p_conf",
    SDE_pooled = "SDE_pooled",
    SESOI_g = "SESOI_g"
    # tau2 intentionally not provided in the simulation -> PI becomes NA (monolab)
  ),
  compute_thresholds = TRUE,
  d33_col = "d33",
  g33_col = "g33"
)

flag_cols <- c(
  "direction_agreement",
  "ci_agreement",
  "ttest_sig",
  "sceptical_p",
  "sceptical_sig",
  "d33",
  "g33",
  "small_telescope_confirmed",
  "sde_confirmed",
  "ci_above_sesoi"
  # "exploratory_within_confirmatory_PI"
)

# Simulation results ------------------------------------------------------

simulation_results[, expl_effect_size_class := factor(effect_class,
                                                      levels = c("Negligible", "Small", "Medium", "Large")
)]


replication_criteria <- c(
  "sde_confirmed",
  "small_telescope_confirmed",
  "ttest_sig",
  "sceptical_sig",
  "direction_agreement",
  "ci_agreement"
)

# Clean criterion names
labels <- c(
  "sde_confirmed" = "Confirmatory ES > mDES",
  "small_telescope_confirmed" = "Small Telescopes",
  "ttest_sig" = "Significant & Same Direction",
  "sceptical_sig" = "Sceptical p-value",
  "direction_agreement" = "Same Direction",
  "ci_agreement" = "Exploratory ES within c-CI"
)

# Summarise by shrinkage and effect class
summary_data_all <- simulation_results[,
                                   lapply(.SD, function(x) mean(as.logical(x), na.rm = TRUE)),
                                   by = .(shrinkage, expl_effect_size_class, exploratory_n),
                                   .SDcols = replication_criteria
]



#### Heatmap of the success rate ####

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
  labels[Criterion],
  levels = labels[replication_criteria]
)]

criterion_order <- c(
  "Same Direction",
  "Small Telescopes",
  "Significant & Same Direction",
  "Sceptical p-value",
  "Exploratory ES within c-CI",
  "Confirmatory ES > mDES"
)
plot_data <- plot_data[labels[Criterion] %in% criterion_order]
plot_data[, Criterion_clean := factor(
  labels[Criterion],
  levels = rev(criterion_order)  # rev() makes first item appear at the top on ggplot y-axis
)]

# Plot for a pre-specified exploratory ss of interest
exp_ss_heatmap_plot <- ggplot(plot_data, aes(x = g_shrinkage, y = Criterion_clean, fill = SuccessRate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", SuccessRate)), size = 5) +
  scale_fill_gradient(low = "white", high = "#434343", limits = c(0, 1), name = "Success Rate") +
  facet_wrap(~ expl_effect_size_class, ncol = 3) +
  coord_equal() +
  labs(
    title = NULL,
    x = "Effect Size Shrinkage",
    y = NULL
  ) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(size = 20, face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing = unit(0.35, "lines"),
    plot.margin = margin(t = 0, r = 5, b = 0, l = 5)
  ) +
  theme_prism()

print(exp_ss_heatmap_plot)

# Save the heatmap for the pre-sepcified exploratory sample size
# heatmap_name <- paste0("combined_heatmap_success_rate_all_classes_expn", exp_ss_heatmap, ".png")
# ggsave(file.path(save_dir, heatmap_name), 
       # exp_ss_heatmap_plot, width = 15, height = 5, dpi = 300)

## Combined heatmap
# Average success rates across all exploratory sample sizes
summary_data_pooled <- summary_data_all[,
                                        lapply(.SD, mean, na.rm = TRUE),
                                        by = .(shrinkage, expl_effect_size_class),
                                        .SDcols = replication_criteria
]

# Prepare for plotting (same pipeline as before)
plot_data_pooled <- melt(summary_data_pooled,
                         id.vars = c("shrinkage", "expl_effect_size_class"),
                         variable.name = "Criterion",
                         value.name = "SuccessRate"
)
setDT(plot_data_pooled)
plot_data_pooled[, g_shrinkage := paste0(round(shrinkage * 100), "%")]
plot_data_pooled <- plot_data_pooled[expl_effect_size_class != "Negligible"]
plot_data_pooled[, Criterion_clean := factor(
  labels[Criterion],
  levels = rev(criterion_order)
)]

combined_heatmap_pooled <- ggplot(plot_data_pooled, 
                                  aes(x = g_shrinkage, y = Criterion_clean, fill = SuccessRate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", SuccessRate)), size = 5) +
  scale_fill_gradient(low = "white", high = "#434343", limits = c(0, 1), name = "Success Rate") +
  facet_wrap(~ expl_effect_size_class, ncol = 3) +
  coord_equal() +
  labs(title = NULL, x = "Effect Size Shrinkage", y = NULL) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text  = element_text(size = 15, face = "bold"),
    panel.spacing = unit(0.35, "lines")
  ) +
  theme_prism()

print(combined_heatmap_pooled)

ggsave(file.path(save_dir, paste0("combined_heatmap_success_rate_all_classes_pooled_n_", dataset_name, ".png")),
       combined_heatmap_pooled, width = 15, height = 5, dpi = 300)
saveRDS(combined_heatmap_pooled,
        here("results", "panels", paste0("combined_heatmap_pooled_", dataset_name, ".rds")))


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

setDT(pr_pooled_data)

# Create clean labels
pr_pooled_data[
  ,
  Criterion_clean := labels[Criterion]
]

# Compute ordered levels (smallest → largest FPR)
ordered_levels <- pr_pooled_data[
  order(mean_FPR),
  unique(Criterion_clean)
]

# Reassign factor with correct order
pr_pooled_data[
  ,
  Criterion_clean := factor(
    Criterion_clean,
    levels = ordered_levels
  )
]

fpr_plot_pooled_n <- ggplot(pr_pooled_data, aes(x = mean_FPR, y = Criterion_clean)) +
  geom_col(alpha = 0.8, width = 0.7) +
  labs(
    title = NULL, #"False Positive Rate by Replication Criterion"
    subtitle = NULL, # "Averaged across exploratory sample sizes (99% shrinkage scenarios)"
    x = "False Positive Rate",
    y = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
    # plot.title = element_text(size = 20, face = "bold"),
    # plot.subtitle = element_text(size = 15, color = "gray60")
  ) +
  theme_prism()

print(fpr_plot_pooled_n)

fpr_plot_pooled_n_name <- paste0("fpr_plot_pooled_n_", dataset_name, ".png")
ggsave(file.path(save_dir, fpr_plot_pooled_n_name), 
       fpr_plot_pooled_n, width = 8, height = 6, dpi = 300)

saveRDS(fpr_plot_pooled_n,
        here("results", "panels", paste0("fpr_plot_pooled_n_", dataset_name, ".rds")))

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
setDT(fpr_heatmap_data)

fpr_heatmap_data[, Criterion_clean := labels[Criterion]]

fpr_heatmap <- ggplot(fpr_heatmap_data, aes(x = Criterion_clean, y = factor(exploratory_n), fill = FPR)) +
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
ggsave(file.path(save_dir, paste0("fpr_heatmap_by_class_", dataset_name, ".png")), fpr_heatmap, 
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
}, by = .(exploratory_n, expl_effect_size_class)]

# Remove negligible
shrinkage_sensitivity <- shrinkage_sensitivity[expl_effect_size_class != "Negligible"]

# Average across conditions
overall_sensitivity <- shrinkage_sensitivity[, .(
  mean_correlation = mean(correlation, na.rm = TRUE),
  mean_slope = mean(slope, na.rm = TRUE),
  se_correlation = sd(correlation, na.rm = TRUE) / sqrt(.N)
), by = criterion]

overall_sensitivity[, criterion_clean := labels[criterion]]
overall_sensitivity <- overall_sensitivity[order(mean_correlation)]  # More negative = more sensitive

print(overall_sensitivity[, .(criterion_clean, mean_correlation, se_correlation)], digits = 3)

# Add clean criterion names
shrinkage_sensitivity[, criterion_clean := labels[criterion]]

# Plot correlation between shrinkage and success rates
sensitivity_plot <- ggplot(shrinkage_sensitivity, 
                           aes(
                             x = correlation, 
                             y = reorder(criterion_clean, correlation),
                             fill = expl_effect_size_class)
                           ) +
  geom_boxplot(width = 0.6) +   
  facet_wrap(~expl_effect_size_class, ncol = 3) +
  scale_fill_brewer(palette = "BuPu", type = "qual", name = "Exploratory Effect Size") +
  labs(
    title = NULL,
    x = "Correlation (Shrinkage ~ Replication Success)",
    y = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme_prism() +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold"),
    legend.position = "none"
  )


print(sensitivity_plot)

ggsave(file.path(save_dir, paste0("shrinkage_sensitivity_", dataset_name, ".png")), sensitivity_plot, 
       width = 12, height = 8, dpi = 300)

saveRDS(sensitivity_plot,
        here("results", "panels", paste0("shrinkage_sensitivity_", dataset_name, ".rds")))

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
pr_summary[, criterion_clean := labels[criterion]]

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

ggsave(file.path(save_dir, paste0("f1_score_comparison_", dataset_name, ".png")), f1_plot, 
       width = 10, height = 8, dpi = 300)

# Precision vs Recall scatter plot
pr_scatter <- ggplot(overall_pr, aes(x = mean_recall, y = mean_precision)) +
  geom_point(size = 4, alpha = 1, color = "#CCCCCC") +
  geom_text_repel(aes(label = criterion_clean), size = 4, max.overlaps = Inf) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  labs(
    title = NULL,
    x = "Sensitivity (Recall)",
    y = "Precision (Positive Predictive Value)"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    # plot.title = element_text(size = 16, face = "bold", margin = margin(b = 15)),
    # plot.subtitle = element_text(size = 15, color = "gray60")
  ) +
  theme_prism()

print(pr_scatter)

ggsave(file.path(save_dir, paste0("precision_recall_scatter_", dataset_name, ".png")), pr_scatter, 
       width = 10, height = 5, dpi = 300)
saveRDS(pr_scatter,
        here("results", "panels", paste0("pr_scatter_", dataset_name, ".rds")))

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

# PR curves by exploratory ss ---------------------------------------------

# --- combine effect sizes: one point per (criterion_clean x sample_size) ---
pr_summary_agg <- pr_summary[, .(
  precision = mean(precision, na.rm = TRUE),
  recall    = mean(recall, na.rm = TRUE),
  f1_score  = mean(f1_score, na.rm = TRUE),
  auc_pr    = mean(auc_pr, na.rm = TRUE),
  n_obs     = sum(n_obs, na.rm = TRUE),
  n_positive = sum(n_positive, na.rm = TRUE),
  n_negative = sum(n_negative, na.rm = TRUE)
), by = .(criterion, criterion_clean, sample_size)]

# --- settings ---
frames_dir <- file.path(save_dir, "pr_frames")
dir.create(frames_dir, showWarnings = FALSE, recursive = TRUE)

n_values <- sort(unique(pr_summary_agg$sample_size))
frame_files <- character(length(n_values))

# --- make frames (one plot per sample_size) ---
for(i in seq_along(n_values)) {
  n <- n_values[i]
  
  plot_dt <- pr_summary_agg[sample_size == n]
  
  p <- ggplot(
    plot_dt,
    aes(x = recall, y = precision, color = criterion_clean)
  ) +
    geom_point(size = 3, alpha = 0.85) +
    ggrepel::geom_text_repel(
      aes(label = criterion_clean),
      size = 3, max.overlaps = Inf, show.legend = FALSE
    ) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent) +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "Precision vs Recall (binary criteria)",
      subtitle = paste0(
        "exploratory_n = ", n,
        " (effect classes averaged; frame ", i, " / ", length(n_values), ")"
      ),
      x = "Recall",
      y = "Precision",
      color = "Criterion"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "bottom")
  
  frame_path <- file.path(frames_dir, sprintf("pr_points_n_%03d.png", i))
  ggsave(frame_path, p, width = 12, height = 8, dpi = 200, bg = "white")
  frame_files[i] <- frame_path
}

# --- stitch frames into a gif ---
gif_path <- file.path(save_dir, "pr_precision_recall_by_n.gif")

gifski::gifski(
  png_files = frame_files,
  gif_file = gif_path,
  width = 1200, height = 800,
  delay = 0.8,
  loop = TRUE
)

gif_path

# Plot PR trajectory across exploratory sample sizes

pr_summary_agg[, sample_size := factor(sample_size, levels = c(5, 10, 15, 20))]
setorder(pr_summary_agg, criterion_clean, sample_size)

p_pr_trajectories <- ggplot(
  pr_summary_agg,
  aes(
    x = recall,
    y = precision,
    group = criterion_clean,
    color = criterion_clean
  )
) +
  geom_path(linewidth = 1.0, alpha = 0.85) +
  geom_point(aes(shape = sample_size), size = 3.2, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = pr_summary_agg[sample_size == "20"],
    aes(label = criterion_clean),
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE,
    box.padding = 1.3,
    point.padding = 1.9
  ) +
  scale_x_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
  scale_y_continuous(limits = c(0.75, 1), labels = percent_format(accuracy = 1)) +
  scale_shape_discrete(name = "Exploratory n") +
  labs(
    title = NULL, # "Precision–Recall trajectories across exploratory sample size",
    subtitle = NULL, # paste("Points are n =", paste(n_keep, collapse = ", "), " (effect classes averaged)"),
    x = "Recall",
    y = "Precision",
    color = "Criterion"
  ) +
  theme_prism(base_size = 14) +
  theme(
    legend.position = "right",
    legend.box = "vertical"
  )

print(p_pr_trajectories)

ggsave(file.path(save_dir, paste0("p_pr_trajectories_", dataset_name, ".png")), p_pr_trajectories, 
       width = 10, height = 5, dpi = 300)

saveRDS(p_pr_trajectories,
        here("results", "panels", paste0("p_pr_trajectories_", dataset_name, ".rds")))

# Criteria landscape: FPR vs PR -------------------------------------------

setDT(overall_pr)
setDT(pr_pooled_data)

# Merge PR with FPR 
criteria_landscape <- merge(
  overall_pr,
  pr_pooled_data,
  by.x = "criterion",
  by.y = "Criterion",
  all.x = TRUE
)

setnames(criteria_landscape, "mean_FPR", "fpr")

criteria_landscape_plot <- ggplot(
  criteria_landscape,
  aes(x = fpr, y = mean_recall)
) +
  geom_point(
    aes(size = mean_precision),
    shape = 21,
    fill = "#6C6C6C",
    # colour = "black",
    # stroke = 0.2,
    alpha = 0.9
  ) +
  ggrepel::geom_text_repel(
    aes(label = criterion_clean),
    size = 4,
    max.overlaps = Inf,
    box.padding = 0,
    point.padding = 17,
    show.legend = FALSE
  ) +
  scale_x_continuous(
    limits = c(0, 0.75),
    labels = percent_format(accuracy = 1)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = percent_format(accuracy = 1)
  ) +
  scale_size_continuous(
    range = c(3, 10),
    labels = percent_format(accuracy = 1),
    name = "Precision"
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),
    name = NULL
  ) +
  labs(
    title = NULL, # "Replication criteria by operating characteristics",
    # subtitle = "x = false positive propensity (99% shrinkage), y = sensitivity (recall at 0–20% shrinkage)",
    x = "False positive propensity (99% shrinkage)",
    y = "Sensitivity (recall at 0–20% shrinkage)"
  ) +
  theme_prism() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(face = "bold")
  )
  

print(criteria_landscape_plot)

ggsave(
  filename = file.path(save_dir, paste0("criteria_landscape_", dataset_name, ".png")),
  plot = criteria_landscape_plot,
  width = 11,
  height = 7,
  dpi = 300,
  bg = "white"
)


