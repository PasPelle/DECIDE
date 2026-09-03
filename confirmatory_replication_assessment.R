# INTRO: Analysis of the DECIDE confirmatory projects including effect size comparison, variance analysis, replication assessment with the replication flags function and effect size deconstruction
library(here)
source(here("packages.R"))

# Import and pre-processing______________####
data_dir        <- here("data") 
save_dir_decide <- here("results")
source(here("all_functions.R"))

# Panels folder for figure building
panels_dir <- file.path(save_dir_decide, "panels")          
if (!dir.exists(panels_dir)) dir.create(panels_dir, recursive = TRUE) 

df <- fread(here("data", "MASTERSHEET_DECIDE_anonymized.csv"), na.strings = c("NA", ""))
all_completed <- copy(df)

# Take only relevant comparisons
all_completed_2 <- all_completed[all_completed$Group_ES %in% c("Ctrl", "Treated")]
all_completed_2[, Value := as.numeric(Value)]

# Remove rows with NA in the observed value (Value column): this can't be done as Project F has pre-computed means so the raw values are NAs
# all_completed_2 <- all_completed_2[!is.na(Value)]

# 1st batch: groups with precomputed means, SD, n  (i.e. Project F exploratory and Project G were extracted from graphs with precomputed means)
batch_1_precomputed_stats <- all_completed_2[
  !is.na(Mean) & !is.na(SD) & !is.na(n),
  .(mean = Mean, 
    sd = SD, 
    n = n),
  by = .(Project_ID, Project_Name, project_letter, Stage, Center, Group_ES)
]

## Process batch 2: groups with raw values and subgroups (Project G, Project H and Project C)

# Check Subgroup Sizes. Filter for Project G and Project C with valid Value
# batch_2_subgroup_sizes <- all_completed_2[
#   !is.na(Subgroup),   # Select all those rows that have subgroups
#   .N,
#   by = .(Project_Name, Stage, Center, Group_ES, Subgroup)
# ][order(Project_Name, Stage, Group_ES, -N)]

## Batch 2: First compute means of subgroups by Project and treatment . 
subgroup_stats <- all_completed_2[
  !is.na(Subgroup),                # select only those rows with subgroup labels
  .(
    mean = mean(Value),            # subgroup means
    sd   = fifelse(.N == 1, 0, sd(Value)),             #  SD, imputation of 0 for those subgroups with n=1
    n    = .N                      # total sample size
  ),
  by = .(Project_ID, Project_Name, project_letter, Stage, Center, Group_ES, Subgroup)
]
# Batch 2: subgroup weighted mean and pooled SD
batch_2_meta_stats <- subgroup_stats[,
                                     .(
                                       mean = weighted.mean(mean, w = n),
                                       sd = sqrt(sum((n - 1) * sd^2)/ sum(n-1)),
                                       n = sum(n)
                                       ),
                                     by = .(Project_ID, Project_Name, project_letter, Stage, Center, Group_ES)
                                     ]
# Check the n of each subgroup
subgroup_n <- df[!is.na(Subgroup), 
                      .(n = .N), 
                      by = .(Project_ID, Stage, Subgroup)]
subgroup_n
# Batch 3: all the other groups with raw data
batch_3_raw_stats <- all_completed_2[
  !is.na(Value) & Value != "" & is.na(Subgroup) & is.na(Mean), # Get only those with raw data (no precomputed mean) and no subgroups
  .(
    mean = mean(Value),
    sd = sd(Value),
    n = .N
    ),
  by = .(Project_ID, Project_Name, project_letter, Stage, Center, Group_ES)
                             ]


# Combine the 3 batches
all_stats <- rbind(batch_1_precomputed_stats, 
                   batch_2_meta_stats, 
                   batch_3_raw_stats)

# Reshape to wide format to compute Hedge's g
all_stats_wide <- all_stats[
  order(Project_ID, Stage, Center, Group_ES)
  ][
    , if (.N == 2) {
      m1 <- mean[Group_ES == "Ctrl"]
      m2 <- mean[Group_ES == "Treated"]
      sd1 <- sd[Group_ES == "Ctrl"]
      sd2 <- sd[Group_ES == "Treated"]
      n1 <- n[Group_ES == "Ctrl"]
      n2 <- n[Group_ES == "Treated"]
      s_pooled <- sqrt(((n1 - 1)*sd1^2 + (n2 - 1)*sd2^2) / (n1 + n2 - 2))
      d <- (m1 - m2) / s_pooled # Cohen's D
      J <- 1 - (3 / (4*(n1 + n2) - 9))  # correction for Hedges' g
      g <- d * J
      se_g <- sqrt((n1 + n2) / (n1 * n2) + (g^2 / (2 * (n1 + n2))))
      
      .(hedges_g = g, se_g = se_g, n1 = n1, n2 = n2, m1 = m1, m2 = m2, sd1 = sd1, sd2 = sd2)
      
    },
    by = .(Project_ID, Project_Name, project_letter, Stage, Center)
  ]

# Add letters column for anonymous plots
# all_stats_wide[, Project_ID := as.numeric(Project_ID)]
# sorted_ids <- sort(unique(all_stats_wide$Project_ID))

write.xlsx(all_stats_wide, file.path(save_dir_decide, "all_stats_wide.xlsx"))


# Create the project_letters lookup table
project_letters <- unique(all_stats_wide[, .(Project_ID, Project_Name, project_letter)])

print("Project Name to Letter Mapping:")
print(project_letters[order(project_letter)])

# Summary treatment arms
df_treatment_arms <- unique(df, by = c("Stage", "Treatment", "Center","project_letter"))

write.xlsx(df_treatment_arms, here("results", "treatment_arms.xlsx"))

# Experimental Units_________________________####

all_stats_wide[, total_n := n1 + n2]

# All EU in exploratory vs confirmatory
all_stats_wide$Stage <- factor(all_stats_wide$Stage, levels = c("exploratory", "confirmatory"))

eu_plot <- ggplot(all_stats_wide, aes(x = Stage, y = total_n, fill = Stage)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, alpha = 0.5, size = 2) +
  labs(
    title = "Experimental Units (EU) by Stage",
    x = NULL,
    y = "EU"
  ) +
  scale_fill_brewer(palette="Set2") +
  theme_minimal(base_size = 10)+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

eu_plot

ggsave(
  filename = file.path(save_dir_decide, "eu_plot.png"),
  plot = eu_plot,
  width = 5,
  height = 6,
  dpi = 50
)

# EU ratio by project (selecting only ctrl and treated groups)

# Summarize confirmatory stats
confirmatory_eu_decide <- all_stats_wide[Stage == "confirmatory", .(
  confirmatory_n = total_n
), by = project_letter]



# Get exploratory total_n 
exploratory_eu_decide <- all_stats_wide[Stage == "exploratory", .(
  exploratory_n = total_n
), by = project_letter]

# Merge and compute ratio
ratio_eu <- merge(exploratory_eu_decide, confirmatory_eu_decide, by = "project_letter")

ratio_eu[, 
  ratio_eu := confirmatory_n / exploratory_n
  ]

# Plot
ratio_eu_plot <- ggplot(ratio_eu, aes(x = project_letter, y = ratio_eu, colour = project_letter)) +
  geom_jitter(width = 0.33, height = 0, size = 4) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(RColorBrewer::brewer.pal(9, "Set1"), "#8B4513")) +

  labs(
    x = "Project",
    y = "Ratio: Confirmatory EU / Exploratory EU"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22)
  )
  

ratio_eu_plot

ggsave(
  filename = file.path(save_dir_decide, "ratio_eu_plot.png"),
  plot = ratio_eu_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Plot planned, executed EU and exploratory EU (not in manuscript)

# average confirmatoy eu by project 
confirmatory_avg_dt <- confirmatory_eu_decide[
  , .(mean_confirmatory_n = mean(confirmatory_n, na.rm = TRUE)),
  by = project_letter
]

# extract planned eu per group and multiply by 2 (ctrl and Treated)
planned_eu <- df[!is.na(planned_EU_total),
                 .(project_letter, planned_EU_per_group)]

planned_eu <- planned_eu[, planned_eu_x2 := planned_EU_per_group * 2]

# replace project G value as the ss is planned disbalanced in the proposal
planned_eu[project_letter == "G", planned_eu_x2 := 28]
  
merged_eu <- merge(exploratory_eu_decide, confirmatory_avg_dt, by = "project_letter")
merged_eu <- merge(merged_eu, planned_eu[, .(planned_eu_x2, project_letter)], 
                   by = "project_letter",
                   all.x = TRUE)

merged_eu_long <- melt(merged_eu, 
                       id.vars = "project_letter",
                       measure.vars = c("exploratory_n", "mean_confirmatory_n", "planned_eu_x2"),
                       variable.name = "Type",
                       value.name = "Count")
setDT(merged_eu_long)

## Plot

# Reorder
merged_eu_long$Type <- factor(merged_eu_long$Type,
                              levels = c("exploratory_n",
                                         "planned_eu_x2",
                                         "mean_confirmatory_n"))

planned_executed_eu_plot <- ggplot(merged_eu_long, 
                                   aes(x = project_letter, 
                                       y = Count, 
                                       color = Type,
                                       shape = Type,
                                       alpha = Type)) +
  geom_point(size = 5) +
  geom_line(aes(group = project_letter), color = "gray70", linewidth = 1) +
  scale_shape_manual(
    name = "Stage",  # Same name as color scale
    values = c("exploratory_n" = 17, 
               "planned_eu_x2" = 16, 
               "mean_confirmatory_n" = 16),
    labels = c("exploratory_n" = "exploratory",
               "planned_eu_x2" = "planned\nconfirmatory",
               "mean_confirmatory_n" = "executed\nconfirmatory")
  ) +
  scale_color_manual(
    name = "Stage",
    values = c("exploratory_n" = "#b3cd7a", 
               "planned_eu_x2" = "#78adcd", 
               "mean_confirmatory_n" = "#4a8db3"),
    labels = c("exploratory_n" = "exploratory",
               "planned_eu_x2" = "planned\nconfirmatory",
               "mean_confirmatory_n" = "executed\nconfirmatory")
  ) +
  scale_alpha_manual(
    values = c("exploratory_n" = 1,
               "planned_eu_x2" = 0.4,
               "mean_confirmatory_n" = 1),
    guide = "none"  # Hide alpha from legend
  ) +
  # scale_y_continuous(breaks = seq(0, 1200, 100),
  #                    limits = c(0, 100)) +
  # scale_y_break(c(300, 1000)) +
  
  labs(
    x = NULL,
    y = "Sample Size"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(alpha = c("exploratory_n" = 1,          # Make the transparency in legend the same as graph
                                    "planned_eu_x2" = 0.4,
                                    "total_executed_EU" = 1))
    )) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(hjust = 1),
    axis.text.y = element_text(hjust = 1),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),  # Distance between legend items
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank()
  ) 

planned_executed_eu_plot

ggsave(
  filename = file.path(save_dir_decide, "planned_executed_eu_plot.png"),
  plot = planned_executed_eu_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Effect Size _________________________####

#### Plot ES per project and stage ####
# # Ensure correct factor levels so projects are ordered by ID
# all_stats_wide[, Project_Label := factor(
#   paste0(Project_ID, "_", Project_Name),
#   levels = unique(paste0(Project_ID, "_", Project_Name))
# )]

# Ensure Stage is a factor with correct left-to-right order
all_stats_wide[, Stage := factor(Stage, levels = c("exploratory", "confirmatory"))]

# Force base x position per project
all_stats_wide[, base_x := as.numeric(factor(project_letter))]

# Count per group to differentiate within confirmatory
all_stats_wide[, stage_index := seq_len(.N), by = .(Project_ID, Stage)]

# add 95% CI
all_stats_wide[,
               `:=`(ci_lower = hedges_g - 1.96 * se_g,
                 ci_upper = hedges_g + 1.96 * se_g
                 )]


# Enhanced spacing logic
all_stats_wide[, jittered_x := fifelse(
  Stage == "exploratory",
  base_x - 0.25,  # left offset for exploratory
  base_x + fifelse(stage_index == 1, 0.12, 0.37)  # separate confirmatory replicates
)]



hedges_effects_by_project <- ggplot(all_stats_wide, aes(x = jittered_x, y = hedges_g, color = Stage, shape = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.1
  ) +
  geom_vline(xintercept = seq(1.5, max(all_stats_wide$base_x) - 0.5, by = 1), 
             linetype = "solid", color = "gray80", linewidth = 0.5) +
  scale_color_manual(values = c("exploratory" = scales::alpha("#003366", 0.9), 
                                "confirmatory" = "#00AFBB"),
                     labels = c(
                       "exploratory" = "Exploratory",
                       "confirmatory" = "Confirmatory"
                     )
                     ) +
  scale_shape_manual(
    values = c(
      "exploratory"  = 17,
      "confirmatory" = 16
    ),
    labels = c(
      "exploratory"  = "Exploratory",
      "confirmatory" = "Confirmatory"
    )
  ) +
  scale_x_continuous(
    breaks = unique(all_stats_wide$base_x),
    labels = unique(all_stats_wide$project_letter)
  ) +
  scale_y_continuous(limits = c(-3, 6)) +
  labs(
    x = "pCS",
    y = "Effect Size (Hedge's g)",
    color = "Stage",
    shape = "Stage"
  ) +
  theme_prism() +
  theme(
    axis.text.x = element_text(hjust = 1),
    axis.text.y = element_text(hjust = 1),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )


hedges_effects_by_project

saveRDS(hedges_effects_by_project, file.path(panels_dir, "hedges_effects_by_project.rds"))

ggsave(
  filename = file.path(save_dir_decide, "hedges_effects_by_project.png"),
  plot = hedges_effects_by_project,
  width = 10,
  height = 6,
  dpi = 300
)

# How many multi-lab projects are null?
all_stats_wide[
  Stage == "confirmatory",
  .(
    total_multi_lab = .N,
    ci_spanning_0   = sum(ci_lower <= 0 & ci_upper >= 0),
    proportion      = mean(ci_lower <= 0 & ci_upper >= 0)
  )
]

# Meta-analysis with fixed-effect model
meta_results <- all_stats_wide[Stage == "confirmatory", {
  res <- run_meta_fixed(hedges_g, se_g, Center)
  extract_meta_fixed(res)
}, by = .(Project_ID, Project_Name, Stage, project_letter)]

# Add back the exploratory data

exploratory_decide <- all_stats_wide[Stage == "exploratory"]
setnames(exploratory_decide, "hedges_g", "g_pooled") # harmonize name columns even though expl g is not pooled (not meta-analysed being 1 lab)
setnames(exploratory_decide, "se_g", "se")

# add the exploratory p-valuesa, same inferential scale as confirmatory (normal z-based)

exploratory_decide[, pval := 2 * pnorm(-abs(g_pooled / se))]

exploratory_rows <- exploratory_decide[, .(
  Project_ID,
  Project_Name,
  Stage,
  project_letter,
  g_pooled,
  se,
  ci_lower, 
  ci_upper,
  pval,
  I2 = NA_real_          # exploratory has no heterogeneity
)]

meta_results <- rbind(meta_results, exploratory_rows, use.names = TRUE, fill = TRUE)


#### ES per project and stage with pooled confirmatory ####

# Prepare meta_results for plotting — pooled confirmatory + exploratory
meta_wide_plot <- meta_results[, .(project_letter, Stage, g_pooled, se, ci_lower, ci_upper)]

# Factor and base x position
meta_wide_plot[, project_letter := factor(project_letter, levels = sort(unique(project_letter)))]
meta_wide_plot[, Stage := factor(Stage, levels = c("exploratory", "confirmatory"))]
meta_wide_plot[, base_x := as.numeric(project_letter)]

# x offset: exploratory left, confirmatory right
meta_wide_plot[, jittered_x := fifelse(
  Stage == "exploratory",
  base_x - 0.15,
  base_x + 0.15
)]

hedges_pooled_by_project <- ggplot(
  meta_wide_plot,
  aes(x = jittered_x, y = g_pooled, color = Stage, shape = Stage)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.08
  ) +
  geom_vline(
    xintercept = seq(1.5, max(meta_wide_plot$base_x) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "exploratory"  = scales::alpha("#003366", 0.9),
      "confirmatory" = "#00AFBB"
    ),
    labels = c(
      "exploratory"  = "Exploratory",
      "confirmatory" = "Confirmatory"
    )
  ) +
  scale_shape_manual(
    values = c(
      "exploratory"  = 17,
      "confirmatory" = 16
    ),
    labels = c(
      "exploratory"  = "Exploratory",
      "confirmatory" = "Confirmatory"
    )
  ) +
  scale_x_continuous(
    breaks = 1:length(levels(meta_wide_plot$project_letter)),
    labels = levels(meta_wide_plot$project_letter)
    # limits = c(0.5, 10.5)
  ) +
  scale_y_continuous(breaks = seq(-4, 6, by = 2)) +
  labs(
    x     = "pCS",
    y     = "Effect Size (Hedges' g)",
    color = "Stage",
    shape = "Stage"
  ) +
  theme_prism() +
  theme(
    axis.text.x        = element_text(hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    legend.position = "bottom"
  )

hedges_pooled_by_project

saveRDS(hedges_pooled_by_project, file.path(panels_dir, "hedges_pooled_by_project.rds"))

ggsave(
  filename = file.path(save_dir_decide, "hedges_pooled_by_project.png"),
  plot     = hedges_pooled_by_project,
  width    = 10,
  height   = 6,
  dpi      = 300
)


meta_long <- melt(
  meta_results,
  id.vars = c("Project_ID", "Project_Name", "Stage", "project_letter"),
  measure.vars = c("g_pooled", "se", "ci_lower", "ci_upper", "pval"),
  variable.name = "metric",
  value.name = "value"
)
# Achieve a wide format table with multiple values per Stage
meta_wide <- dcast(
  meta_long,
  Project_ID + Project_Name + project_letter ~ Stage + metric,
  value.var = "value"
)
meta_wide <- as.data.table(meta_wide)

# Exploratory g is not a pooled g
setnames(meta_wide, "exploratory_g_pooled", "exploratory_g")

# # Sceptical p-value
# # Compute z-scores and variance ratio
# meta_wide[, `:=`(                                   # := operator has to be in backticks here because it's used as a fx to create multiple columns
#   exploratory_z = exploratory_g / exploratory_se,
#   confirmatory_z = confirmatory_g_pooled / confirmatory_se,
#   variance_ratio = (exploratory_se / confirmatory_se)^2
# )]
# 
# 
# # Compute sceptical p-values
# meta_wide[, sceptical_p := pSceptical(
#   zo = exploratory_z,
#   zr = confirmatory_z,
#   c = variance_ratio,
#   alternative = "two.sided",
#   type = "golden"
# )]
# 
# 
# # Flag significant results
# meta_wide[, sceptical_sig := sceptical_p < 0.05]

# Calculate max value with some padding
max_g <- max(abs(meta_wide$exploratory_g), abs(meta_wide$confirmatory_g_pooled), na.rm = TRUE)
max_limit <- max_g * 1.1  # Add 10% padding

XY_scatterplot_ES <- ggplot(meta_wide, aes(x = abs(exploratory_g), y = abs(confirmatory_g_pooled))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(size = 3, color ="darkgreen", alpha = 0.8) +
  geom_text_repel(aes(label = project_letter),
                  size = 4,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  segment.color = "grey50"
                  ) +
  # scale_color_manual(
  #   values = c("TRUE" = "darkcyan", "FALSE" = "deeppink"),
  #   labels = c("TRUE" = "Sceptical-Significant", "FALSE" = "Not Significant"),
  #   name = "Confirmatory\nSceptical p < 0.05"
  # ) +
  coord_fixed(ratio = 1, xlim = c(0, max_limit), ylim = c(0, max_limit), clip = "off") +
  labs(
    x = "|Exploratory Hedges' g|",
    y = "|Confirmatory Hedges' g|"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    # legend.text = element_text(size = 22),
    # legend.title = element_text(size = 22),
    plot.margin = margin(20, 20, 20, 20)  # Add margins around the plot
  )

XY_scatterplot_ES

pdf(here("results", "XY_scatterplot_ES.pdf"),
    width = 8, height = 8)
XY_scatterplot_ES # Plot goes into PDF 
dev.off()  # Close PDF

ggsave("XY_scatterplot_ES.png", plot = XY_scatterplot_ES,
       path = here("results"),
       width = 5, height = 6, dpi = 200)



# Variance Analysis_________________________####

#### Heterogeneity indicators ####
# Confirmatory lab heterogeneity
# Basic stats
confirmatory <- all_stats_wide[Stage == "confirmatory",]

project_stats <- confirmatory[, {
  m <- rma(yi = hedges_g, sei = se_g, method = "REML", data = .SD)
  list(tau2 = m$tau2, I2 = m$I2, QEp = m$QEp)
}, by = .(Project_ID, Project_Name, project_letter), .SDcols = c("hedges_g", "se_g")]


project_stats[, `:=`(
  tau2 = formatC(tau2, format = "e", digits = 2),    # Format e is the Exp notation
  I2   = formatC(I2,   format = "f", digits = 1),    # Format f is fixed point notation (good for percentages, ES)
  QEp  = formatC(QEp,  format = "e", digits = 2)
)]

fwrite(project_stats, here("results", "project_stats.csv"))


# # Plot I2 heterogeneity in confirmatory projects
project_stats[, I2_numeric := as.numeric(I2)]
project_stats[, project_letter := factor(project_letter, levels = rev(sort(unique(project_letter))))]

heterogeneity_plot <- ggplot(project_stats, aes(x = I2_numeric, y = project_letter)) +
  geom_bar(
    stat = "identity",
    fill = "#3182bd",
    width = 0.7,
    color = "black",
    linewidth = 0.3
  ) +
  geom_vline(xintercept = c(25, 50, 75), linetype = "dashed", color = "gray50", linewidth = 0.3) +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    x = expression(paste("Heterogeneity (", I^2, " %)")),
    y = NULL,
    title = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

heterogeneity_plot

# Calculate range and spread metrics
# center_variability <- confirmatory[, {
#   m <- rma(yi = hedges_g, sei = se_g, method = "REML", data = .SD)
#   
#   list(
#     n_centers = .N,
#     pooled_g = m$b[1],
#     min_g = min(hedges_g),
#     max_g = max(hedges_g),
#     range_g = max(hedges_g) - min(hedges_g),  # Simple spread
#     sd_g = sd(hedges_g),  # SD of observed effects
#     cv_g = sd(hedges_g) / abs(mean(hedges_g))  # Coefficient of variation
#   )
# }, by = .(project_letter, Project_Name), .SDcols = c("hedges_g", "se_g")]
# 
# print(center_variability[order(project_letter)])
# 
# cv_plot <- ggplot(center_variability, 
#                   aes(x = cv_g, y = project_letter)) +
#   geom_bar(stat = "identity", fill = "#d62728", width = 0.7,
#            color = "black", linewidth = 0.3) +
#   geom_vline(xintercept = c(20, 40, 60), linetype = "dashed", 
#              color = "gray50", linewidth = 0.3) +
#   geom_text(aes(label = paste0("n=", n_centers)), 
#             hjust = -0.2, size = 3) +
#   scale_x_continuous(
#     expand = expansion(mult = c(0, 0.08))) +
#   coord_cartesian(xlim = c(0, 15)) +
#   labs(
#     x = "Coefficient of Variation (%)",
#     y = NULL
#   ) +
#   theme_classic(base_size = 12) +
#   theme(
#     panel.grid = element_blank(),
#     axis.line = element_line(color = "black", linewidth = 0.5),
#     axis.ticks = element_line(color = "black", linewidth = 0.5),
#     axis.ticks.length = unit(0.15, "cm"),
#     axis.text = element_text(color = "black", size = 10),
#     axis.title.x = element_text(margin = margin(t = 10), size = 11),
#     plot.margin = margin(10, 10, 10, 10)
#   )
# 
# cv_plot
# 
# ggsave("cv_plot.pdf", plot = cv_plot,
#        path = here("results"),
#        width = 7, height = 5, device = cairo_pdf)

# ggsave("cv_plot.png", plot = cv_plot,
#        path = here("results"),
#        width = 5, height = 6, dpi = 200)
# 

#### CV controls exploratory vs pooled confirmatory #####

# Compare variance of exploratory controls vs pooled confirmatory controls

# Exploratory: single-lab cv for controls
ctrl_exploratory <- all_stats[Group_ES == "Ctrl" & Stage == "exploratory",
                              .(mean_expl = mean,
                                sd_expl = sd,
                                n_expl = n,
                                cv_exploratory = sd / mean),
                              by = .(project_letter, Project_Name)]

# Confirmatory: pooled observation across labs (Column Center) for controls and then cv
ctrl_confirmatory <- all_stats[Group_ES == "Ctrl" & Stage == "confirmatory",
                               {
                                 pooled_mean <- mean(mean)
                                 pooled_sd   <- sqrt(sum((n - 1) * sd^2 + n * (mean - pooled_mean)^2) / (sum(n) - 1))
                                 .(mean_pooled = pooled_mean,
                                   sd_pooled   = pooled_sd,
                                   n_total     = sum(n))
                               },
                               by = .(project_letter, Project_Name)]

ctrl_confirmatory[, cv_confirmatory := sd_pooled / abs(mean_pooled)]

merged_ctrl <- ctrl_exploratory[ctrl_confirmatory, on = .(project_letter, Project_Name)]
cv_merged_ctrl <- merged_ctrl[, .(project_letter, Project_Name, cv_exploratory, cv_confirmatory)]


var_compare_long <- melt(cv_merged_ctrl,
                         measure.vars = c("cv_exploratory", "cv_confirmatory"),
                         variable.name = "Stage",
                         value.name = "cv")
setDT(var_compare_long)
var_compare_long[, Stage := fifelse(Stage == "cv_exploratory", "exploratory", "confirmatory")]


cv_decide_pooled_labs_ctrl_plot <- ggplot(
  var_compare_long,
  aes(x = project_letter, y = cv, color = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, length(unique(cv_merged_ctrl$project_letter)) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("exploratory" = scales::alpha("#003366", 0.9), "confirmatory" = "#00AFBB"),
    labels = c("exploratory" = "Exploratory", "confirmatory" = "Confirmatory")
  ) +
  labs(title = "Control", x = "Project", y = "Variance (CV)", color = "Stage") +
  theme_prism() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

cv_decide_pooled_labs_ctrl_plot


### CV difference expl vs conf pooled labs ###

# Compute CV difference (confirmatory - exploratory) per project, center, group
cv_diff_pooled_dt <- cv_merged_ctrl[
  ,
  .(cv_diff = cv_confirmatory - cv_exploratory),
  by = .(project_letter)
]

# Test if cv_diff is significantly different from 0 per Group_ES
cv_diff_pooled_test <- t.test(cv_diff_pooled_dt$cv_diff, mu = 0)
cv_diff_pooled_test
wilcox.test(cv_diff_pooled_dt$cv_diff, mu = 0)  #

# Add significance label for annotation
p_val <- cv_diff_pooled_test$p.value
sig_label <- fcase(
  p_val < 0.001, "p < 0.001",
  p_val < 0.01,  "p < 0.01",
  p_val < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_val, 3))
)

# BoxPlot
cv_diff_ctrl_pooled_plot <- ggplot(cv_diff_pooled_dt, aes(x = "", y = cv_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#003366", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#003366", width = 0.15, size = 2) +
  geom_text(
    data = cv_diff_pooled_dt,
    aes(x = Inf, y = Inf, label = sig_label),
    hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE
  )+
  labs(x = NULL, y = "CV difference \n(pooled confirmatory - exploratory)", title = "Control") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

cv_diff_ctrl_pooled_plot


# Per-project test: is CV_confirmatory > CV_exploratory?
cv_test_per_project <- merged_ctrl[,
                                   {
                                     # Approximate SE of each CV
                                     se_cv_expl  <- cv_exploratory / sqrt(2 * n_expl)  * sqrt(1 + 2 * cv_exploratory^2)
                                     se_cv_conf  <- cv_confirmatory / sqrt(2 * n_total) * sqrt(1 + 2 * cv_confirmatory^2)

                                     # SE of the difference
                                     se_diff <- sqrt(se_cv_expl^2 + se_cv_conf^2)

                                     # t-statistic for the difference (confirmatory - exploratory)
                                     cv_diff <- cv_confirmatory - cv_exploratory
                                     t_stat  <- cv_diff / se_diff

                                     # Degrees of freedom (Welch-Satterthwaite)
                                     df <- (se_cv_expl^2 + se_cv_conf^2)^2 /
                                       (se_cv_expl^4 / (n_expl - 1) + se_cv_conf^4 / (n_total - 1))

                                     # One-tailed p-value: is confirmatory CV *larger* than exploratory?
                                     p_value <- pt(-abs(t_stat), df = df)  # one-tailed

                                     .(cv_diff, t_stat, df, p_value)
                                   },
                                   by = .(project_letter, Project_Name)
]

# Significance label per project
cv_test_per_project[, sig_label := fcase(
  p_value < 0.001, "p < 0.001",
  p_value < 0.01,  "p < 0.01",
  p_value < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_value, 3))
)]

# Overall test across projects
overall_cv_test <- t.test(cv_test_per_project$cv_diff, mu = 0, alternative = "greater")
p_val_overall <- overall_cv_test$p.value
sig_label_overall <- fcase(
  p_val_overall < 0.001, "p < 0.001",
  p_val_overall < 0.01,  "p < 0.01",
  p_val_overall < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_val_overall, 3))
)

# Plot
cv_diff_per_project_plot <- ggplot(cv_test_per_project, aes(x = "", y = cv_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#003366", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#003366", width = 0.15, size = 2) +
  labs(x = NULL, y = "CV difference\n(confirmatory - exploratory)", title = "Control") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
cv_diff_per_project_plot

cv_decide_pooled_labs_ctrl_plot | cv_diff_per_project_plot


#### CV by initiating lab #####

# Create a within-lab replication dataset

# subset, for each project (project_letter), the ones where the exploratory and confirmatory Stage have the same Center
in_lab_repl_dt <- all_stats[
  all_stats[,
            .(Center_ref = intersect(
              Center[Stage == "exploratory"],
              Center[Stage == "confirmatory"])),
              by = project_letter],
            on = .(project_letter, Center = Center_ref)
  ]


# compute coefficient of variation
in_lab_repl_dt <- in_lab_repl_dt[, cv := sd / abs(mean),
  by = project_letter]


cv_decide_in_lab_plot_ctrl <- ggplot(
  in_lab_repl_dt[Group_ES == "Ctrl"],
  aes(x = project_letter, y = cv, color = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, length(unique(in_lab_repl_dt$project_letter)) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("exploratory" = scales::alpha("#003366", 0.9), "confirmatory" = "#00AFBB"),
    labels = c("exploratory" = "Exploratory", "confirmatory" = "Confirmatory")
  ) +
  labs(title = "Control", x = "Project", y = "Variance (CV)", color = "Stage") +
  theme_prism() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
cv_decide_in_lab_plot_ctrl

cv_decide_in_lab_plot_treated <- ggplot(
  in_lab_repl_dt[Group_ES == "Treated"],
  aes(x = project_letter, y = cv, color = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, length(unique(in_lab_repl_dt$project_letter)) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("exploratory" = scales::alpha("#003366", 0.9), "confirmatory" = "#00AFBB"),
    labels = c("exploratory" = "Exploratory", "confirmatory" = "Confirmatory")
  ) +
  labs(title = "Treated", x = "Project", y = "Variance (CV)", color = "Stage") +
  theme_prism() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
cv_decide_in_lab_plot_treated


### CV difference initiating labs ###

# Compute CV difference (confirmatory - exploratory) per project, center, group
cv_diff_dt <- in_lab_repl_dt[
  ,
  .(cv_diff = cv[Stage == "confirmatory"] - cv[Stage == "exploratory"]),
  by = .(project_letter, Center, Group_ES)
]

# Test if cv_diff is significantly different from 0 per Group_ES
cv_diff_test <- cv_diff_dt[
  ,
  {
    test <- t.test(cv_diff, mu = 0)
    .(p_value = test$p.value,
      mean_diff = test$estimate,
      ci_lower = test$conf.int[1],
      ci_upper = test$conf.int[2])
  },
  by = Group_ES
]

# Add significance label for annotation
cv_diff_test[, sig_label := fcase(
  p_value < 0.001, "p < 0.001",
  p_value < 0.01,  "p < 0.01",
  p_value < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_value, 3))
)]


# BoxPlot
cv_diff_plot_ctrl <- ggplot(cv_diff_dt[Group_ES == "Ctrl"],
                            aes(x = project_letter, y = cv_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#003366", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#003366", width = 0.15, size = 2) +
  geom_text(
    data = cv_diff_test[Group_ES == "Ctrl"],
    aes(x = Inf, y = Inf, label = sig_label),
    hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE
  ) +
  labs(title = "Control", x = "Project", y = "CV difference (confirmatory - exploratory)") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(hjust = 1))

cv_diff_plot_ctrl


# Boxplot difference initiating labs in ctrl only
cv_diff_plot_pooled_ctrl <- ggplot(cv_diff_dt[Group_ES == "Ctrl"],
                                   aes(x = "", y = cv_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#003366", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#003366", width = 0.15, size = 2) +
  geom_text(
    data = cv_diff_test[Group_ES == "Ctrl"],
    aes(x = Inf, y = Inf, label = sig_label),
    hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE
  )+
  labs(x = NULL, y = "CV difference \n(confirmatory - exploratory)", title = "Control") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

cv_diff_plot_pooled_ctrl

# Side by side plot (with patchwork library)
# Note: | is side by side, / is stacked. You can combine them too, e.g. (A | B) / C.
cv_ctrl_plot <- cv_decide_in_lab_plot_ctrl | cv_diff_plot_pooled_ctrl
cv_ctrl_plot

ggsave(
  filename = file.path(save_dir_decide, "cv_decide_initiating_labs_ctrl.png"),
  plot = cv_ctrl_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Boxplot difference initiating labs in ctrl and treated together
cv_diff_plot_pooled <- ggplot(cv_diff_dt, aes(x = Group_ES, y = cv_diff, fill = Group_ES)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Group_ES), width = 0.15, size = 2) +
  geom_text(
    data = cv_diff_test,
    aes(x = Group_ES, y = Inf, label = sig_label),
    vjust = 1.5, size = 4, inherit.aes = FALSE
  ) +
  scale_fill_manual(values  = c("Ctrl" = "#003366", "Treated" = "#00AFBB")) +
  scale_color_manual(values = c("Ctrl" = "#003366", "Treated" = "#00AFBB")) +
  labs(x = NULL, y = "CV difference (confirmatory - exploratory)") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

cv_diff_plot_pooled

#### Mean difference exploratory vs pooled confirmatory controls ####

# merged_ctrl <- ctrl_exploratory[ctrl_confirmatory, on = .(project_letter, Project_Name)]

merged_ctrl_long <-  melt(merged_ctrl,
                          measure.vars = c("mean_expl", "mean_pooled"),
                          variable.name = "Stage",
                          value.name = "mean")
setDT(merged_ctrl_long)
# Rename stage values
merged_ctrl_long[, Stage := fifelse(Stage == "mean_expl", "exploratory", "confirmatory")]

# Plot exploratory mean vs confirmatory pooled means

mean_decide_pooled_labs_ctrl_plot <- ggplot(
  merged_ctrl_long,
  aes(x = project_letter, y = mean, color = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, length(unique(cv_merged_ctrl$project_letter)) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("exploratory" = scales::alpha("#003366", 0.9), "confirmatory" = "#00AFBB"),
    labels = c("exploratory" = "Exploratory", "confirmatory" = "Confirmatory")
  ) +
  labs(title = "Control", x = "Project", y = "Mean", color = "Stage") +
  theme_prism() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

mean_decide_pooled_labs_ctrl_plot

# Mean difference
# merged_ctrl[, mean_diff := mean_pooled - mean_expl] #subtract confirmatory pooled mean

# Hedges' g per project
merged_ctrl_hedges <- merged_ctrl[,
                           {
                             # Mean difference and SD ratio
                             mean_diff <- mean_pooled - mean_expl
                             sd_ratio  <- sd_pooled / sd_expl

                             # Pooled SD for Hedges' g
                             sd_pooled_hedges <- sqrt(((n_expl - 1) * sd_expl^2 + (n_total - 1) * sd_pooled^2) /
                                                        (n_expl + n_total - 2))

                             # Cohen's d
                             d <- (mean_pooled - mean_expl) / sd_pooled_hedges

                             # Hedges' g correction factor
                             df  <- n_expl + n_total - 2
                             J   <- 1 - (3 / (4 * df - 1))
                             g   <- d * J

                             # SE of Hedges' g
                             se_g <- sqrt((n_expl + n_total) / (n_expl * n_total) + g^2 / (2 * (n_expl + n_total)))

                             # t-statistic and two-tailed p-value
                             t_stat  <- g / se_g
                             p_value <- 2 * pt(-abs(t_stat), df = df)

                             .(mean_diff, sd_ratio, sd_pooled, sd_expl, g, se_g, t_stat, p_value)
                           },
                           by = .(project_letter, Project_Name)
]

# Significance label for overall one-sample t-test (is g different from 0 across projects?)
overall_test <- t.test(merged_ctrl_hedges$g, mu = 0)
p_val <- overall_test$p.value

merged_ctrl_hedges[, sig_label := fcase(
  p_value < 0.001, "p < 0.001",
  p_value < 0.01,  "p < 0.01",
  p_value < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_value, 3))
)]

hedges_g_plot_ctrl <- ggplot(merged_ctrl_hedges, aes(x = "", y = g)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#003366", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#003366", width = 0.15, size = 2) +
  labs(x = NULL, y = "Hedges' g\n(confirmatory - exploratory)", title = "Control") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())
hedges_g_plot_ctrl

mean_decide_pooled_labs_ctrl_plot | hedges_g_plot_ctrl

#### Plot CV difference and mean difference per project ####

# Merge CV per-project test and Hedges' g into one long table for forest plot
cv_forest <- cv_test_per_project[, .(project_letter, Project_Name,
                                     estimate = cv_diff,
                                     se = (cv_diff / t_stat),  # back-calculate SE
                                     metric = "CV difference")]

g_forest <- merged_ctrl_hedges[, .(project_letter, Project_Name,
                            estimate = g,
                            se = se_g,
                            metric = "Hedges' g")]

forest_dt <- rbind(cv_forest, g_forest)

# Compute 95% CI
forest_dt[, ci_lower := estimate - 1.96 * se]
forest_dt[, ci_upper := estimate + 1.96 * se]

# Order projects consistently
forest_dt[, project_letter := factor(project_letter,
                                     levels = sort(unique(project_letter)))]

# Forest plot
forest_plot <- ggplot(forest_dt,
                      aes(x = estimate, y = project_letter, color = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_pointrange(
    aes(xmin = ci_lower, xmax = ci_upper),
    position = position_dodge(width = 0.5),
    size = 0.6
  ) +
  scale_color_manual(
    values = c("CV difference" = scales::alpha("#003366", 0.9),
               "Hedges' g"    = "#00AFBB"),
    name = NULL
  ) +
  labs(
    x = "Effect size (with 95% CI)",
    y = "Project",
    title = "CV difference vs Hedges' g per project"
  ) +
  theme_prism() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_blank()
  )

forest_plot



# CV difference of the initiating lab (from cv_diff_dt, Ctrl only)
cv_inlab <- cv_diff_dt[Group_ES == "Ctrl",
                       .(cv_diff_inlab = round(cv_diff, 3)),
                       by = project_letter]

# CV difference expl vs pooled confirmatory (from cv_test_per_project)
cv_pooled <- cv_test_per_project[,
                                 .(project_letter,
                                   cv_diff_pooled = round(cv_diff, 3),
                                   p_cv = round(p_value, 3))]

# Hedges' g (from merged_ctrl_hedges)
hedges <- merged_ctrl_hedges[,
                      .(project_letter,
                        g = round(g, 3),
                        p_g = round(p_value, 3))]

# Merge all
summary_table <- cv_inlab[cv_pooled, on = "project_letter"]
summary_table <- summary_table[hedges, on = "project_letter"]

summary_table[, g_label := paste0(g, " (", fcase(
  p_g < 0.001, "p < 0.001",
  p_g < 0.01,  "p < 0.01",
  p_g < 0.05,  "p < 0.05",
  default = paste0("p = ", p_g)
), ")")]

# Format p-values inline
summary_table[, cv_diff_pooled_label := paste0(cv_diff_pooled, " (", fcase(
  p_cv < 0.001, "p < 0.001",
  p_cv < 0.01,  "p < 0.01",
  p_cv < 0.05,  "p < 0.05",
  default = paste0("p = ", p_cv)
), ")")]

summary_table[, conclusion := fcase(
  p_cv < 0.05 & p_g < 0.05,  "Both significant",
  p_cv < 0.05 & p_g >= 0.05, "CV diff only",
  p_cv >= 0.05 & p_g < 0.05, "Hedges' g only",
  default = "Neither significant"
)]


display_table <- summary_table[, .(
  Project                         = project_letter,
  `CV diff (initiating lab)`      = cv_diff_inlab,
  `CV diff (expl vs pooled conf)` = cv_diff_pooled_label,
  `Hedges' g`                     = g_label,
  `Conclusion`                    = conclusion
)]

setorder(display_table, Project)

display_table |>
  gt() |>
  tab_header(title = "Variance summary — Controls") |>
  cols_align(align = "center") |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Conclusion",
      rows = Conclusion == "Both significant"
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#fff3cd"),  # yellow
    locations = cells_body(
      columns = "Conclusion",
      rows = Conclusion %in% c("CV diff only", "Hedges' g only")
    )
  ) |>
  tab_style(
    style = cell_fill(color = "#d4edda"),
    locations = cells_body(
      columns = "Conclusion",
      rows = Conclusion == "Neither significant"
    )
  )


# # Compute mean difference and SD ratio per project
# merged_ctrl[, `:=`(
#   mean_diff = mean_pooled - mean_expl,
#   sd_ratio  = sd_pooled / sd_expl
# )]

# Reshape to long for forest-style plot
decomp_long <- melt(
  merged_ctrl_hedges[, .(project_letter, Project_Name, mean_diff, sd_ratio, g)],
  id.vars       = c("project_letter", "Project_Name"),
  variable.name = "metric",
  value.name    = "value"
)
setDT(decomp_long)

decomp_long[, metric := fcase(
  metric == "mean_diff", "Mean difference (numerator)",
  metric == "sd_ratio",  "SD ratio (denominator)",
  metric == "g",         "Hedges' g (combined)"
)]

decomp_long[, project_letter := factor(project_letter, levels = rev(sort(unique(project_letter))))]
decomp_long[, ref := fifelse(metric == "SD ratio (denominator)", 1, 0)]

decomp_plot <- ggplot(decomp_long, aes(x = value, y = project_letter, color = metric)) +
  geom_vline(aes(xintercept = ref), linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Mean difference (numerator)" = scales::alpha("#003366", 0.9),
    "SD ratio (denominator)"      = "#CC3300",
    "Hedges' g (combined)"        = "#00AFBB"
  )) +
  facet_wrap(~ metric, scales = "free_x", ncol = 3) +
  labs(x = NULL, y = "Project", title = "Decomposition of effect size — Controls") +
  theme_prism() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        panel.grid.major.y = element_line(color = "gray90"))

decomp_plot

# Replication assessment________________####

# Add more columns for the assessment: exploratory n1 and n2 to meta_wide
exploratory_n <- all_stats_wide[Stage == "exploratory", .(n1 = unique(n1), n2 = unique(n2)), by = Project_Name]
setnames(exploratory_n, c("n1", "n2"), c("n1_exploratory", "n2_exploratory"))
meta_wide <- merge(exploratory_n[,.(
  Project_Name,
  n1_exploratory,
  n2_exploratory
)
],
meta_wide,
by = "Project_Name")

# # Add more columns for the assessment: inter-study variance column (tau2)
# meta_wide <- merge(
#   meta_wide,
#   project_stats[, .(Project_Name, tau2_confirmatory = as.numeric(tau2))],
#   by = "Project_Name",
#   all = TRUE
# )


# Add more columns for the assessment: pooled Smallest Effec Size SDE
# Compute SDE
confirmatory[, sde := NA_real_]  # initialize column

confirmatory[n1 >= 2 & n2 >= 2, sde := mapply(function(n1, n2) {
  pwr::pwr.t2n.test(n1 = n1, n2 = n2, sig.level = 0.05, power = 0.8)$d
}, n1, n2)]


# Convert d to g
confirmatory[, sde_g := sde *(1- (3 / (4*(n1 + n2) - 9)))]

# Compute SE for each SDE
confirmatory[, se_sde_g := sqrt((n1 + n2)/(n1 * n2) + (sde_g^2 / (2 * (n1 + n2))))]

# Meta-analyze SDEs across centers to get pooled SDE threshold
sde_meta <- confirmatory[
  ,
  {
    res <- run_meta_fixed(
      te      = sde_g,
      se      = se_sde_g,
      studlab = Center
    )
    list(
      SDE_pooled = res$TE.common,
      SDE_se     = res$seTE.common,
      SDE_lower  = res$lower.common,
      SDE_upper  = res$upper.common
    )
  },
  by = Project_Name,
  .SDcols = c("sde_g", "se_sde_g", "Center")
]
# Add the SDE column
meta_wide <- merge(
  meta_wide,
  sde_meta[, .(Project_Name, SDE_pooled)],
  by = "Project_Name",
  all.x = TRUE
)

# Add more columns for the assessment:  import SESOI extracted from BMBF proposals

SESOI <- fread(here("data", "SESOI_DECIDE_anonymized.csv"), sep = ";", dec = ",", encoding = "Latin-1")

SESOI[, SESOI_g := as.numeric(SESOI_g)]
# merge with confirmatory
meta_wide <- merge(
  meta_wide,
  SESOI[, .(Project_Name, SESOI_g)],
  by = "Project_Name",
  all.x = TRUE
)

# === Run replication assessment function ===
meta_wide <- compute_replication_flags(
  meta_wide,
  cfg = list(alpha = 0.05, z_pi = 1.96, power_d33 = 0.33, power_sde = 0.8,
             sceptical_type = "golden", sceptical_alternative = "two.sided"),
  cols = list(
    id = "Project_Name",
    g_exp = "exploratory_g", se_exp = "exploratory_se",
    ci_lo_exp = "exploratory_ci_lower", ci_hi_exp = "exploratory_ci_upper", p_exp = "exploratory_pval",
    n1_exp = "n1_exploratory", n2_exp = "n2_exploratory",
    g_conf = "confirmatory_g_pooled", se_conf = "confirmatory_se",
    ci_lo_conf = "confirmatory_ci_lower", ci_hi_conf = "confirmatory_ci_upper", p_conf = "confirmatory_pval",
    tau2_conf = "tau2_confirmatory",
    SDE_pooled = "SDE_pooled",
    SESOI_g = "SESOI_g"
  ),
  compute_thresholds = TRUE
)


#### Replication matrix Confirmatory ####
replication_matrix_decide <- meta_wide[, .(
  Project_ID,
  Project_Name,
  project_letter,
  `Exploratory g in confirmatory CI` = ci_agreement,
  # `Exploratory g in confirmatory PI` = exploratory_within_confirmatory_PI,
  `Exploratory and confirmatory in the same direction` = direction_agreement,
  `Sceptical p < 0.05` = sceptical_sig,
  `Both p < 0.05 & Same Direction` = ttest_sig,
  `Small Telescopes` = small_telescope_confirmed,
  `Confirmatory larger than its SDE` = sde_confirmed,
  `CI above SESOI` = ci_above_sesoi
  # `lnRR confirmatory same direction and >= exploratory` = lnRR_confirmation
)]

replication_long_decide <- melt(replication_matrix_decide, 
                                id.vars = c("Project_ID", "Project_Name", "project_letter"),
                                variable.name = "criterion", 
                                value.name = "success")
replication_long_decide <- as.data.table(replication_long_decide)

# Count TRUEs per project to rank them
replication_long_decide[, success_count := sum(success == TRUE, na.rm = TRUE), by = Project_ID]

#  Use project letter as factor and order by success count
replication_long_decide[, project_letter := factor(as.character(project_letter), levels = sort(unique(as.character(project_letter))))]

replication_long_decide[, success_char := as.character(success)]

# # Order projects by ID
# replication_long_decide[, Project_Name := factor(Project_Name,
#                                           levels = unique(meta_wide$Project_Name))]

# Store replication results
saveRDS(replication_long_decide, file.path(save_dir_decide, "replication_long_decide.rds"))


# Heatmap
heatmap_plot_decide <- ggplot(replication_long_decide, aes(x = project_letter, y = criterion, fill = success_char)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.85) +  # Added width/height < 1 to shrink tiles
  
  # More professional color scheme (Nature/Cell style)
  scale_fill_manual(
    values = c("TRUE" = "#1F6F70", "FALSE" = "#D7D7D7"),
    labels = c("TRUE" = "Met", "FALSE" = "Not met"),
    name = "Criterion status"
  ) +
  
  # Cleaner, more descriptive labels
  scale_y_discrete(
    labels = c(
      "CI above SESOI" = "CI above SESOI",
      "Confirmatory larger than its SDE" = "Confirmatory ES > mDES",
      "Small Telescopes" = "Small Telescopes",
      "Both p < 0.05 & Same Direction" = "Significant & Same Direction",
      "Sceptical p < 0.05" = "Sceptical p-value",
      "Exploratory and confirmatory in the same direction" = "Same Direction",
      # "Exploratory g in confirmatory PI" = "Exploratory ES within c-PI ",
      "Exploratory g in confirmatory CI" = "Exploratory ES within c-CI"
    )
  ) +
  
  labs(
    x = "Confirmatory Studies",
    y = NULL,
    title = NULL
  ) +
  
  # Professional theme with LARGER text
  theme_minimal(base_size = 14) +  
  theme(
    # Axis text - 
    axis.text.x = element_text(size = 20, face = "bold", color = "black"),  
    axis.text.y = element_text(size = 20, hjust = 1, color = "black"),      
    
    # Axis titles
    axis.title.x = element_text(size = 17, face = "bold", margin = margin(t = 10)),
    
    # Legend - also bigger
    legend.position = "bottom",
    legend.title = element_text(size = 21, face = "bold"),
    legend.text = element_text(size = 20),
    legend.key.size = unit(0.6, "cm"),
    legend.box.spacing = unit(0.3, "cm"),
    
    # Panel
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    # Plot margins - increased to accommodate larger text
    plot.margin = margin(10, 10, 10, 10, "pt"),
    
    # Background
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  theme_prism()

heatmap_plot_decide

# Save with adjusted dimensions for larger labels
ggsave(
  filename = file.path(save_dir_decide, "replication_heatmap_large_labels.png"),
  plot = heatmap_plot_decide,
  width = 16,   
  height = 11,  
  dpi = 300,
  bg = "white"
)


# Effect Size deconstruction__________________####

#### Shrinkage indicators - normalized mean difference (not included in the manuscript) ####
# Indicators of Hedges g, mean diff, sd, ctrl mean diff, ctrl sd
# with formula |exploratory| - |confirmatory| / |exploratory| + |confirmatory|
# < 0: confirmatory inflated
# 0: same magnitude
# 0 - 1: confirmatory smaller than exploratory
# 1: total shrinkage
# >1: sign error

mean_diff_dt <- all_stats[,
                          .(mean_diff = mean[Group_ES == "Ctrl"] - mean[Group_ES == "Treated"],
                            # variance sum law for independent variables to get the SD of thedifference:
                            # SD of the difference is done by summing the variances (squared SD of the single SD) and then taking the square root to return to the SD
                            sd_diff   = sqrt(sd[Group_ES == "Ctrl"]^2 + sd[Group_ES == "Treated"]^2),
                            n_ctrl    = n[Group_ES == "Ctrl"],
                            n_treated = n[Group_ES == "Treated"]),
                          by = .(project_letter, Project_Name, Stage, Center)
]

mean_diff_confirmatory <- mean_diff_dt[Stage == "confirmatory", {
  pooled_mean <- mean(mean_diff)
  pooled_sd   <- sqrt(sum((n_ctrl + n_treated - 2) * sd_diff^2 +
                            (n_ctrl + n_treated) * (mean_diff - pooled_mean)^2) /
                        (sum(n_ctrl + n_treated) - 1))
  .(mean_diff = pooled_mean, sd_pooled = pooled_sd,
    n_total = sum(n_ctrl + n_treated), Stage = "confirmatory")
}, by = .(project_letter, Project_Name)]

mean_diff_exploratory <- mean_diff_dt[Stage == "exploratory",
                                      .(project_letter, Project_Name, Stage, mean_diff,
                                        sd_pooled = sd_diff, n_total = n_ctrl + n_treated)]

mean_diff_bind <- rbind(mean_diff_exploratory, mean_diff_confirmatory, fill = TRUE)

g_to_join <- meta_results[, .(project_letter,
                              Stage = as.character(Stage),
                              g_pooled, se, ci_lower, ci_upper, pval)]

mean_diff_bind <- g_to_join[mean_diff_bind, on = .(project_letter, Stage)]

shrinkage_dt <- mean_diff_bind[,
                               .(
                                 g_shrinkage    = (abs(g_pooled[Stage == "exploratory"])  - abs(g_pooled[Stage == "confirmatory"])) /
                                   (abs(g_pooled[Stage == "exploratory"])  + abs(g_pooled[Stage == "confirmatory"])),
                                 mean_shrinkage = (abs(mean_diff[Stage == "exploratory"]) - abs(mean_diff[Stage == "confirmatory"])) /
                                   (abs(mean_diff[Stage == "exploratory"]) + abs(mean_diff[Stage == "confirmatory"])),
                                 # Note: this is the SD of the difference ctrl - treated
                                 sd_shrinkage   = (abs(sd_pooled[Stage == "exploratory"]) - abs(sd_pooled[Stage == "confirmatory"])) /
                                   (abs(sd_pooled[Stage == "exploratory"]) + abs(sd_pooled[Stage == "confirmatory"]))
                               ),
                               by = .(project_letter, Project_Name)
]

merged_ctrl_shrinkage <- merged_ctrl[
  ctrl_exploratory[, .(project_letter, mean_expl)], on = "project_letter"
][
  ctrl_confirmatory[, .(project_letter, mean_pooled)], on = "project_letter"
]

merged_ctrl_shrinkage[, `:=`(
  ctrl_mean_shrinkage = (abs(mean_expl) - abs(mean_pooled)) / (abs(mean_expl) + abs(mean_pooled)),
  ctrl_sd_shrinkage   = (abs(sd_expl)   - abs(sd_pooled))   / (abs(sd_expl)   + abs(sd_pooled))
)]

# Project factor order (shared across all panels)
project_levels <- rev(sort(unique(shrinkage_dt$project_letter)))


# ── Panel 1: Hedges g (ctrl vs treated) ──────────────────────────────────────
p1_dt <- shrinkage_dt[, .(project_letter, value = g_shrinkage,
                          metric = "Effect Site Shrinkage")]
p1_dt[, project_letter := factor(project_letter, levels = project_levels)]

p1 <- ggplot(p1_dt, aes(x = value, y = project_letter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = "#00AFBB") +
  facet_wrap(~ metric) +
  labs(x = NULL, y = "Confirmatory Studies") +
  shrinkage_theme()
p1
# ── Panel 2: Mean diff (ctrl vs treated) ─────────────────────────────────────
p2_dt <- shrinkage_dt[, .(project_letter, value = mean_shrinkage,
                          metric = "Mean Shrinkage")]
p2_dt[, project_letter := factor(project_letter, levels = project_levels)]

p2 <- ggplot(p2_dt, aes(x = value, y = project_letter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = scales::alpha("#003366", 0.9)) +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  shrinkage_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# ── Panel 3: SD (ctrl vs treated) ────────────────────────────────────────────
p3_dt <- shrinkage_dt[, .(project_letter, value = sd_shrinkage,
                          metric = "Variance Inflation")]
p3_dt[, project_letter := factor(project_letter, levels = project_levels)]

p3 <- ggplot(p3_dt, aes(x = value, y = project_letter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = "#CC3300") +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  shrinkage_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# ── Panel 4: Mean diff (ctrl expl vs conf) ───────────────────────────────────
p4_dt <- merged_ctrl_shrinkage[, .(project_letter, value = ctrl_mean_shrinkage,
                                   metric = "Control Stability")]
p4_dt[, project_letter := factor(project_letter, levels = project_levels)]

p4 <- ggplot(p4_dt, aes(x = value, y = project_letter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = "#E7B800") +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  shrinkage_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# ── Panel 5: Distance method (Niangoran 2023) — mean difference as variable Y ──
# Variable Y = mean difference (ctrl - treated) per confirmatory center
# Already computed in mean_diff_dt: mean_diff and sd_diff per center
# SS_i = (n_total_i - 1)*sd_diff_i^2 + n_total_i*(mean_diff_i - grand_mean_diff)^2
# Di = (SS_i / df1_i) / (SS_total / df2) ~ F(df1_i, df2)

alpha_dist <- 0.05

distance_decide_dt <- mean_diff_dt[Stage == "confirmatory"][,
                                                            {
                                                              # Total n per center: ctrl + treated
                                                              n_total_i <- n_ctrl + n_treated
                                                              
                                                              # Grand mean difference: weighted by total n per center
                                                              y_grand   <- sum(n_total_i * mean_diff) / sum(n_total_i)
                                                              
                                                              # Total sample size across all centers in this project
                                                              N_total   <- sum(n_total_i)
                                                              
                                                              # Numerator df per center: total n minus 1
                                                              df1       <- n_total_i - 1
                                                              
                                                              # Denominator df: total observations minus 1
                                                              df2       <- N_total - 1
                                                              
                                                              # Per-center SS: within-center component (from sd_diff) + between-center component
                                                              SS_i      <- (n_total_i - 1) * sd_diff^2 +
                                                                n_total_i * (mean_diff - y_grand)^2
                                                              
                                                              # Total SS across all centers
                                                              SS_total  <- sum(SS_i)
                                                              
                                                              # Di: center's mean squared deviation relative to overall mean squared deviation
                                                              Di        <- (SS_i / df1) / (SS_total / df2)
                                                              
                                                              # Critical F value per center
                                                              f_crit    <- qf(1 - alpha_dist, df1 = df1, df2 = df2)
                                                              
                                                              # Flag atypical centers
                                                              atypical  <- Di > f_crit
                                                              
                                                              .(Center, n_ctrl, n_treated, n_total_i, mean_diff, sd_diff,
                                                                y_grand, Di, f_crit, atypical)
                                                            },
                                                            by = .(project_letter, Project_Name)
]

# Align factor levels with other panels
distance_decide_dt[, project_letter := factor(project_letter, levels = project_levels)]

# Strip label to match other panels
distance_decide_dt[, metric := "Atypical Labs"]

distance_decide_dt[, atypical := factor(atypical, levels = c("FALSE", "TRUE"))]
distance_decide_dt[, atypical := as.logical(as.character(atypical))]

p5_distance_decide <- ggplot(
  distance_decide_dt,
  aes(x = Di, y = project_letter, color = atypical)
) +
  # Di = 1: center exactly at the grand mean difference
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("FALSE" = "#333333", "TRUE" = "#B8B8B8"),
    labels = c("FALSE" = "Typical Lab", "TRUE" = "Atypical Lab (p < 0.05)")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous(
    limits = c(0, 2),         
    breaks = seq(0, 2, by = 0.5)
  ) +
  # facet_wrap(~ metric) +
  labs(x = "Distance Statistic (Dᵢ)", y = "pCS", title = NULL) +
  shrinkage_theme() +
  theme(
    # axis.text.y     = element_blank(),
    # axis.ticks.y    = element_blank(),
    # axis.line.y     = element_blank(),
    # strip.text = element_blank(),
    legend.position = "bottom",
    legend.title    = element_blank(),
    plot.title = element_blank(),
    axis.ticks.length = unit(2.75, "pt"), 
    axis.title.x = element_text(margin = margin(t = 4))
  )

p5_distance_decide

saveRDS(p5_distance_decide, file.path(panels_dir, "p5_distance_decide.rds"))   

ggsave(
  filename = file.path(save_dir_decide, "atypical_labs_p_decide.png"),
  plot     = p5_distance_decide,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# # Exploratory sd_diff for each project
# mean_diff_dt[Stage == "exploratory", 
#              .(project_letter, Center, sd_diff_expl = sd_diff)]
# 
# # Confirmatory sd_diff per center for each project
# mean_diff_dt[Stage == "confirmatory",
#              .(project_letter, Center, sd_diff_conf = sd_diff)]
# 
# # Side by side for project A specifically
# mean_diff_dt[project_letter == "A",
#              .(Stage, Center, mean_diff, sd_diff, n_ctrl, n_treated)]


# ── Combine all 5 ────────────────────────────────────────────────────────────
shrinkage_decomp_decide_plot <- p1 + p2 + p3 + p4 + p5_distance_decide +
  plot_layout(nrow = 1, widths = c(1.4, 1, 1, 1, 1), guides = "collect") &
  theme(legend.position = "none")

shrinkage_decomp_decide_plot

ggsave(
  filename = file.path(save_dir_decide, "shrinkage_decomp_decide_plot.png"),
  plot     = shrinkage_decomp_decide_plot,
  width    = 14,
  height   = 6,
  dpi      = 300
)

## Summary table of all shrinkage indicators
# pull directly from the source tables
shrinkage_summary_decide <- shrinkage_dt[,
                                         .(project_letter,
                                           Project_Name,
                                           g_shrinkage    = round(g_shrinkage, 3),
                                           mean_shrinkage = round(mean_shrinkage, 3),
                                           sd_shrinkage   = round(sd_shrinkage, 3))
]

# Add raw expl/conf values from mean_diff_bind
expl_conf_decide <- mean_diff_bind[,
                                   .(project_letter,
                                     Stage,
                                     g_pooled  = round(g_pooled, 3),
                                     mean_diff = round(mean_diff, 3),
                                     sd_pooled = round(sd_pooled, 3))
]

expl_decide <- expl_conf_decide[Stage == "exploratory",
                                .(project_letter,
                                  g_expl        = g_pooled,
                                  mean_diff_expl = mean_diff,
                                  sd_expl       = sd_pooled)
]

conf_decide <- expl_conf_decide[Stage == "confirmatory",
                                .(project_letter,
                                  g_conf        = g_pooled,
                                  mean_diff_conf = mean_diff,
                                  sd_conf       = sd_pooled)
]

shrinkage_summary_decide <- shrinkage_summary_decide[
  expl_decide, on = "project_letter"
][
  conf_decide, on = "project_letter"
]

# Ctrl mean shrinkage (from merged_ctrl_shrinkage)
ctrl_mean_summary_decide <- merged_ctrl_shrinkage[,
                                                  .(project_letter,
                                                    ctrl_mean_expl      = round(mean_expl, 3),
                                                    ctrl_mean_conf      = round(mean_pooled, 3),
                                                    ctrl_mean_shrinkage = round(ctrl_mean_shrinkage, 3))
]

# Distance Di: one row per project
distance_summary_decide <- distance_decide_dt[,
                                              .(n_centers    = .N,
                                                n_atypical   = sum(atypical == "TRUE"),
                                                median_Di    = round(median(Di), 3)),
                                              by = project_letter
]

# Merge all
summary_table_decide <- shrinkage_summary_decide[
  ctrl_mean_summary_decide, on = "project_letter"
][
  distance_summary_decide, on = "project_letter"
]

setorder(summary_table_decide, project_letter)

# Display with gt
summary_table_decide_dt <- summary_table_decide[, .(
  Project                  = project_letter,
  `g (expl)`               = g_expl,
  `g (conf)`               = g_conf,
  `Hedges g shrinkage`     = g_shrinkage,
  `Mean diff (expl)`       = mean_diff_expl,
  `Mean diff (conf)`       = mean_diff_conf,
  `Mean diff shrinkage`    = mean_shrinkage,
  `SD (expl)`              = sd_expl,
  `SD (conf)`              = sd_conf,
  `SD shrinkage`           = sd_shrinkage,
  `Ctrl mean (expl)`       = ctrl_mean_expl,
  `Ctrl mean (conf)`       = ctrl_mean_conf,
  `Ctrl mean shrinkage`    = ctrl_mean_shrinkage,
  `N centers`              = n_centers,
  `Median Di`              = median_Di,
  `N atypical centers`     = ifelse(n_centers <= 2, "—", as.character(n_atypical))
)]

summary_table_decide_dt |>
  gt() |>
  tab_header(title = "Shrinkage indicators — DECIDE dataset") |>
  tab_spanner(
    label   = "Hedges g",
    columns = c("g (expl)", "g (conf)", "Hedges g shrinkage")
  ) |>
  tab_spanner(
    label   = "Mean difference (ctrl vs treated)",
    columns = c("Mean diff (expl)", "Mean diff (conf)", "Mean diff shrinkage")
  ) |>
  tab_spanner(
    label   = "SD (ctrl vs treated)",
    columns = c("SD (expl)", "SD (conf)", "SD shrinkage")
  ) |>
  tab_spanner(
    label   = "Ctrl mean (expl vs conf)",
    columns = c("Ctrl mean (expl)", "Ctrl mean (conf)", "Ctrl mean shrinkage")
  ) |>
  tab_spanner(
    label   = "Distance method",
    columns = c("N centers", "Median Di", "N atypical centers")
  ) |>
  cols_align(align = "center") |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Hedges g shrinkage",
      rows    = `Hedges g shrinkage` > 0.5 | `Hedges g shrinkage` < 0
    )
  ) |>
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Mean diff shrinkage",
      rows    = `Mean diff shrinkage` > 0.5 | `Mean diff shrinkage` < 0
    )
  ) |>
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "SD shrinkage",
      rows    = `SD shrinkage` < -0.1
    )
  ) |>
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Ctrl mean shrinkage",
      rows    = `Ctrl mean shrinkage` > 0.2 | `Ctrl mean shrinkage` < -0.2
    )
  ) |>
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "N atypical centers",
      rows    = `N atypical centers` != "—" & `N atypical centers` != "0"
    )
  ) |>
  tab_footnote(
    footnote  = "— indicates only 2 confirmatory centers: Distance method not informative with < 3 centers.",
    locations = cells_column_labels(columns = "N atypical centers")
  )

# Save as Excel
write.xlsx(
  summary_table_decide_dt,
  file.path(save_dir_decide, "shrinkage_summary_decide.xlsx"),
  rowNames = FALSE
)

plot_dt <- shrinkage_summary_decide[ctrl_mean_summary_decide, on = "project_letter"]


#### Factor contribution: deltas (not included in the manuscript)####

# Isolated counterfactual source contributions to Hedges' g shrinkage
# Three independent substitutions, each starting from g_exp:
#
#   delta_ctrl    : replace μ_ctrl_exp    with μ_ctrl_conf    (control stability)
#   delta_treated : replace μ_treated_exp with μ_treated_conf (treatment response)
#   delta_sd      : replace SD_exp        with SD_conf        (variance inflation)
#
# Each delta is on the raw Δg scale (not normalized).
# Negative delta = that source reduced g → contributing to shrinkage
# Positive delta = that source increased g → protective against shrinkage
#
# NOTE: the three deltas are independent and do NOT sum to total shrinkage.

# ── Exploratory components (single lab) ──────────────────────────────────────

expl_raw <- all_stats[Stage == "exploratory" & Group_ES %in% c("Ctrl", "Treated"),
                      .(mean = mean, sd = sd, n = n, Group_ES),
                      by = .(project_letter, Project_Name)]

expl_ctrl <- expl_raw[Group_ES == "Ctrl",
                      .(project_letter, Project_Name,
                        mu_ctrl_exp = mean, sd_ctrl_exp = sd, n_ctrl_exp = n)]

expl_treated <- expl_raw[Group_ES == "Treated",
                         .(project_letter, Project_Name,
                           mu_treated_exp = mean, sd_treated_exp = sd, n_treated_exp = n)]

expl_components <- expl_ctrl[expl_treated, on = .(project_letter, Project_Name)]

expl_components[, sd_exp_pooled := sqrt(
  ((n_ctrl_exp - 1) * sd_ctrl_exp^2 + (n_treated_exp - 1) * sd_treated_exp^2) /
    (n_ctrl_exp + n_treated_exp - 2)
)]

expl_components[, g_exp := hedges_g_from_parts(
  mu_ctrl_exp, mu_treated_exp, sd_exp_pooled, n_ctrl_exp, n_treated_exp
)]

# ── Confirmatory components (pooled across labs) ──────────────────────────────

# conf_ctrl <- all_stats[Stage == "confirmatory" & Group_ES == "Ctrl", {
#   pm  <- mean(mean)
#   psd <- sqrt(sum((n - 1) * sd^2 + n * (mean - pm)^2) / (sum(n) - 1))
#   .(mu_ctrl_conf = pm, sd_ctrl_conf = psd, n_ctrl_conf = sum(n))
# }, by = .(project_letter, Project_Name)]
# 
# conf_treated <- all_stats[Stage == "confirmatory" & Group_ES == "Treated", {
#   pm  <- mean(mean)
#   psd <- sqrt(sum((n - 1) * sd^2 + n * (mean - pm)^2) / (sum(n) - 1))
#   .(mu_treated_conf = pm, sd_treated_conf = psd, n_treated_conf = sum(n))
# }, by = .(project_letter, Project_Name)]

# Weighted means
conf_ctrl <- all_stats[Stage == "confirmatory" & Group_ES == "Ctrl", {
  w   <- n / sum(n)
  pm  <- sum(w * mean)
  psd <- sqrt(sum((n - 1) * sd^2 + n * (mean - pm)^2) / (sum(n) - 1))
  .(mu_ctrl_conf = pm, sd_ctrl_conf = psd, n_ctrl_conf = sum(n))
}, by = .(project_letter, Project_Name)]

conf_treated <- all_stats[Stage == "confirmatory" & Group_ES == "Treated", {
  w   <- n / sum(n)
  pm  <- sum(w * mean)
  psd <- sqrt(sum((n - 1) * sd^2 + n * (mean - pm)^2) / (sum(n) - 1))
  .(mu_treated_conf = pm, sd_treated_conf = psd, n_treated_conf = sum(n))
}, by = .(project_letter, Project_Name)]

conf_components <- conf_ctrl[conf_treated, on = .(project_letter, Project_Name)]

conf_components[, sd_conf_pooled := sqrt(
  ((n_ctrl_conf - 1) * sd_ctrl_conf^2 + (n_treated_conf - 1) * sd_treated_conf^2) /
    (n_ctrl_conf + n_treated_conf - 2)
)]

conf_components[, g_conf := hedges_g_from_parts(
  mu_ctrl_conf, mu_treated_conf, sd_conf_pooled, n_ctrl_conf, n_treated_conf
)]

# ── Single decomp_dt — built once, columns added progressively ───────────────

decomp_dt <- expl_components[conf_components, on = .(project_letter, Project_Name)]

# Total shrinkage on |g| scale — used throughout
decomp_dt[, total_shrinkage := abs(g_exp) - abs(g_conf)]

# Shared plot settings
component_colors <- c(
  "Control stability"  = "#8E7CC3",   
  "Treatment response" = "#6FA287",
  "Variance inflation" = "#B5563C"
)

project_ord <- rev(sort(unique(decomp_dt$project_letter)))
# ── Independent counterfactual substitutions ─────────────────────────


# Each source swapped independently, always starting from g_exp.
# The other two sources remain at exploratory values.
# Normalized by |g_exp| so projects with different effect sizes are comparable.
# Δ / |g_exp| = -1 means that source alone would have eliminated the entire effect.
#
# Problem: deltas do NOT sum to total_shrinkage.
# The gap is the interaction term — the joint effect of sources changing
# simultaneously, invisible when looking at each source in isolation.

decomp_dt[, `:=`(
  g_c_ind = hedges_g_from_parts(
    mu_ctrl_conf, mu_treated_exp, sd_exp_pooled, n_ctrl_exp, n_treated_exp
  ),
  g_t_ind = hedges_g_from_parts(
    mu_ctrl_exp, mu_treated_conf, sd_exp_pooled, n_ctrl_exp, n_treated_exp
  ),
  g_s_ind = hedges_g_from_parts(
    mu_ctrl_exp, mu_treated_exp, sd_conf_pooled, n_ctrl_conf, n_treated_conf
  )
)]

# Positive = contributing to shrinkage (reduced g), negative = protective
decomp_dt[, `:=`(
  delta_ctrl_ind    = -(g_c_ind - g_exp) / abs(g_exp),
  delta_treated_ind = -(g_t_ind - g_exp) / abs(g_exp),
  delta_sd_ind      = -(g_s_ind - g_exp) / abs(g_exp),
  interaction_ind   = total_shrinkage / abs(g_exp) -
    (-(g_c_ind - g_exp) / abs(g_exp)) -
    (-(g_t_ind - g_exp) / abs(g_exp)) -
    (-(g_s_ind - g_exp) / abs(g_exp))
)]

message("Independent deltas (normalized by |g_exp|):")
print(decomp_dt[, .(
  project_letter,
  total_shrinkage   = round(total_shrinkage, 3),
  delta_ctrl        = round(delta_ctrl_ind, 3),
  delta_treated     = round(delta_treated_ind, 3),
  delta_sd          = round(delta_sd_ind, 3),
  sum_deltas        = round(delta_ctrl_ind + delta_treated_ind + delta_sd_ind, 3),
  interaction       = round(interaction_ind, 3)
)])

# Plot
deltas_long <- melt(
  decomp_dt[, .(
    project_letter,
    `Control stability`      = delta_ctrl_ind,
    `Treatment response` = delta_treated_ind,
    `Variance inflation` = delta_sd_ind
  )],
  id.vars = "project_letter", variable.name = "component", value.name = "delta"
)
setDT(deltas_long)
set(deltas_long, j = "project_letter", value = factor(deltas_long$project_letter, levels = project_ord))
set(deltas_long, j = "component", value = factor(deltas_long$component,
                                                levels = c("Control stability", "Treatment response", "Variance inflation")))

p_deltas <- ggplot(deltas_long, aes(x = delta, y = project_letter, fill = component)) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = component_colors) +
  scale_x_continuous(
    sec.axis = dup_axis(
      breaks = c(-0.5, 0.5),
      labels = c("← protective", "contributing →"),
      name   = NULL
    )
  ) +
  labs(
    x       = "Δ|g| / |g_exp|",
    y       = NULL,
    fill    = NULL,
    title   = "Step 1: Independent counterfactual substitutions",
    caption = paste(
      "Each bar: isolated effect of one source, holding the other two at exploratory values.",
      "Normalized by |g_exp|. Bars do NOT sum to total shrinkage — interaction term not shown.",
      sep = "\n"
    )
  ) +
  theme_prism() +
  theme(
    legend.position  = "bottom",
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.line.x.top  = element_blank(),
    panel.border     = element_blank()
  )

p_deltas

ggsave(
  filename = file.path(save_dir_decide, "shrinkage_deltas_independent.png"),
  plot = p_deltas, width = 10, height = 7, dpi = 300
)
#### Factor contribution: sequential deltas (not included in the manuscript)####
# Fixed substitution order: ctrl mean → treated mean → SD
# Each step is the MARGINAL contribution given all previous swaps:
#
#   delta_ctrl    = g_ctrl_swap     - g_exp          (swap ctrl mean only)
#   delta_treated = g_both_means    - g_ctrl_swap    (swap treated mean, ctrl already swapped)
#   delta_sd      = g_conf          - g_both_means   (swap SD, both means already swapped)
#
# By construction: delta_ctrl + delta_treated + delta_sd = g_conf - g_exp
# i.e. the three terms fill the total shrinkage bar exactly.
#
# Sign convention:
#   total_shrinkage = g_exp - g_conf  (positive = shrinkage, negative = inflation)
#   each delta is negated so positive segment = contributing to shrinkage
#
# Ordering rationale: baseline stability first (ctrl), then treatment response,
# then variance — matching the biological logic of the three shrinkage indicators.

# Ordering 1: ctrl → SD → treated
decomp_dt[, `:=`(
  g_c_seq1  = hedges_g_from_parts(
    mu_ctrl_conf, mu_treated_exp, sd_exp_pooled,  n_ctrl_exp,  n_treated_exp   # ctrl only
  ),
  g_cs_seq1 = hedges_g_from_parts(
    mu_ctrl_conf, mu_treated_exp, sd_conf_pooled, n_ctrl_conf, n_treated_conf  # ctrl + SD
  )
)]

decomp_dt[, `:=`(
  delta_ctrl_seq1    = -(g_c_seq1  - g_exp),        # ctrl first
  delta_sd_seq1      = -(g_cs_seq1 - g_c_seq1),     # SD second, ctrl already in
  delta_treated_seq1 = -(g_conf    - g_cs_seq1)      # treated last
)]

# Ordering 2: SD → ctrl → treated
decomp_dt[, `:=`(
  g_s_seq2  = hedges_g_from_parts(
    mu_ctrl_exp,  mu_treated_exp, sd_conf_pooled, n_ctrl_conf, n_treated_conf  # SD only
  ),
  g_sc_seq2 = hedges_g_from_parts(
    mu_ctrl_conf, mu_treated_exp, sd_conf_pooled, n_ctrl_conf, n_treated_conf  # SD + ctrl
  )
)]

decomp_dt[, `:=`(
  delta_sd_seq2      = -(g_s_seq2  - g_exp),        # SD first
  delta_ctrl_seq2    = -(g_sc_seq2 - g_s_seq2),     # ctrl second, SD already in
  delta_treated_seq2 = -(g_conf    - g_sc_seq2)      # treated last
)]

# Sanity checks
decomp_dt[, `:=`(
  check_seq1 = delta_ctrl_seq1 + delta_sd_seq1 + delta_treated_seq1,
  check_seq2 = delta_ctrl_seq2 + delta_sd_seq2 + delta_treated_seq2
)]

message("STEP 2 — Sequential deltas (check1 and check2 should both equal g_exp - g_conf):")
print(decomp_dt[, .(project_letter,
                    total  = round(g_exp - g_conf, 3),
                    check1 = round(check_seq1, 3),
                    check2 = round(check_seq2, 3),
                    delta_ctrl1 = round(delta_ctrl_seq1, 3),
                    delta_ctrl2 = round(delta_ctrl_seq2, 3),
                    delta_sd1   = round(delta_sd_seq1, 3),
                    delta_sd2   = round(delta_sd_seq2, 3))])

# Combine both orderings for plot
ord1 <- decomp_dt[, .(
  project_letter,
  ordering             = "ctrl → SD → treated",
  `Control stability`      = delta_ctrl_seq1,
  `Treatment response` = delta_treated_seq1,
  `Variance inflation` = delta_sd_seq1
)]

ord2 <- decomp_dt[, .(
  project_letter,
  ordering             = "SD → ctrl → treated",
  `Control stability`      = delta_ctrl_seq2,
  `Treatment response` = delta_treated_seq2,
  `Variance inflation` = delta_sd_seq2
)]

ord_long <- melt(
  rbind(ord1, ord2),
  id.vars       = c("project_letter", "ordering"),
  variable.name = "component",
  value.name    = "delta"
)
setDT(ord_long)

set(ord_long, j = "project_letter", value = factor(
  ord_long$project_letter, levels = project_ord
))

set(ord_long, j = "component", value = factor(
  ord_long$component,
  levels = c("Control stability", "Treatment response", "Variance inflation")
))

set(ord_long, j = "ordering", value = factor(
  ord_long$ordering,
  levels = c("ctrl → SD → treated", "SD → ctrl → treated")
))

p_sequential <- ggplot(ord_long, aes(x = delta, y = project_letter, fill = component)) +
  geom_col(position = "stack", width = 0.7) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  facet_wrap(~ ordering) +
  scale_fill_manual(values = component_colors) +
  labs(
    x       = "Hedges' g",
    y       = NULL,
    fill    = NULL,
    title   = "Step 2: Sequential substitution — attribution changes with ordering",
    caption = paste(
      "Total bar length identical in both panels — telescoping sum guarantees this.",
      "Colour composition changes because ctrl and SD enter g multiplicatively (SD is the denominator).",
      "There is no biologically motivated correct ordering — this arbitrariness motivates Step 3.",
      sep = "\n"
    )
  ) +
  theme_prism() +
  theme(
    legend.position    = "bottom",
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank(),
    axis.ticks.x.top   = element_blank(),
    axis.line.x.top    = element_blank(),
    panel.border       = element_blank(),
    strip.text         = element_text(size = 10, face = "bold")
  )

p_sequential

ggsave(
  filename = file.path(save_dir_decide, "shrinkage_sequential_sequential.png"),
  plot     = p_sequential,
  width    = 12,
  height   = 7,
  dpi      = 300
)


### Factor contribution: Shapley values ####

# Total shrinkage = |g_exp| - |g_conf| (reduction in absolute effect size)
# Positive = effect shrank, negative = effect inflated
#
# Three sources: ctrl mean drift, treatment response, variance inflation
# Shapley values computed on |g| so all projects on same interpretable scale
#
# Guarantees:
#   phi_ctrl + phi_treated + phi_sd = |g_exp| - |g_conf| (exact, no residual)
#   No arbitrary ordering assumption
#   Each source gets fair attribution
#
# Sign convention:
#   positive phi = that source contributed to shrinkage (reduced |g|)
#   negative phi = that source was protective (increased |g|)


# All 7 coalition |g| values
decomp_dt[, `:=`(
  g0    = abs(g_exp),
  
  # Single swaps
  g_c   = abs(hedges_g_from_parts(mu_ctrl_conf,  mu_treated_exp,  sd_exp_pooled,
                                  n_ctrl_exp,  n_treated_exp)),
  g_t   = abs(hedges_g_from_parts(mu_ctrl_exp,   mu_treated_conf, sd_exp_pooled,
                                  n_ctrl_exp,  n_treated_exp)),
  g_s   = abs(hedges_g_from_parts(mu_ctrl_exp,   mu_treated_exp,  sd_conf_pooled,
                                  n_ctrl_conf, n_treated_conf)),
  
  # Double swaps
  g_ct  = abs(hedges_g_from_parts(mu_ctrl_conf,  mu_treated_conf, sd_exp_pooled,
                                  n_ctrl_exp,  n_treated_exp)),
  g_cs  = abs(hedges_g_from_parts(mu_ctrl_conf,  mu_treated_exp,  sd_conf_pooled,
                                  n_ctrl_conf, n_treated_conf)),
  g_ts  = abs(hedges_g_from_parts(mu_ctrl_exp,   mu_treated_conf, sd_conf_pooled,
                                  n_ctrl_conf, n_treated_conf)),
  
  # Triple swap
  g_cts = abs(g_conf)
)]

# Shapley values (before negation: negative = reduced |g| = shrinkage)
decomp_dt[, `:=`(
  phi_ctrl = (2/6) * (g_c   - g0)  +
    (1/6) * (g_ct  - g_t) +
    (1/6) * (g_cs  - g_s) +
    (2/6) * (g_cts - g_ts),
  
  phi_treated = (2/6) * (g_t   - g0)  +
    (1/6) * (g_ct  - g_c) +
    (1/6) * (g_ts  - g_s) +
    (2/6) * (g_cts - g_cs),
  
  phi_sd = (2/6) * (g_s   - g0)  +
    (1/6) * (g_cs  - g_c) +
    (1/6) * (g_ts  - g_t) +
    (2/6) * (g_cts - g_ct)
)]

# Negate: positive = contributed to |g| reduction = shrinkage
decomp_dt[, `:=`(
  phi_ctrl    = -phi_ctrl,
  phi_treated = -phi_treated,
  phi_sd      = -phi_sd
)]

# Sanity check
decomp_dt[, check_shapley := phi_ctrl + phi_treated + phi_sd]
message("STEP 3 — Shapley values (diff should be 0 for all rows):")
print(decomp_dt[, .(
  project_letter,
  total_shrinkage = round(total_shrinkage, 3),
  phi_ctrl        = round(phi_ctrl, 3),
  phi_treated     = round(phi_treated, 3),
  phi_sd          = round(phi_sd, 3),
  check           = round(check_shapley, 3),
  diff            = round(total_shrinkage - check_shapley, 8)
)])

# Plot
shapley_long <- melt(
  decomp_dt[, .(
    project_letter,
    `Control stability`      = phi_ctrl,
    `Treatment response` = phi_treated,
    `Variance inflation` = phi_sd
  )],
  id.vars = "project_letter", variable.name = "component", value.name = "phi"
)
setDT(shapley_long)
set(shapley_long, j = "project_letter", value = factor(shapley_long$project_letter, levels = project_ord))
set(shapley_long, j = "component", value = factor(shapley_long$component,
                                                levels = c("Control stability", "Treatment response", "Variance inflation")))

p_shapley <- ggplot(shapley_long, aes(x = phi, y = project_letter, fill = component)) +
  geom_col(position = "stack", width = 0.7) +
  geom_point(
    data = decomp_dt,
    aes(x = total_shrinkage, 
        y = factor(as.character(project_letter), levels = levels(shapley_long$project_letter)), 
        fill = "Total shrinkage", shape = "Total shrinkage"),
    inherit.aes = FALSE,
    size = 2.5,
    color = "black"
  ) +
  scale_shape_manual(name = NULL, 
                     values = c("Total shrinkage" = 23),
                     labels = c("Total shrinkage" = "Total Shrinkage")
                     ) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  scale_fill_manual(
    values = c(component_colors, "Total shrinkage" = "black"),
    breaks = c("Total shrinkage", names(component_colors)),
    labels = c(
      "Total shrinkage"     = "Total shrinkage",
      "Control stability"   = "Control Stability",
      "Treatment response"  = "Treatment Response",
      "Variance inflation"  = "Variance Inflation"
    )
  ) +
  scale_x_continuous(
    breaks = seq(-2.5, 7.5, 2.5),
    sec.axis = dup_axis(
      breaks = c(-1.5, 3),
      labels = c("← reducing", "contributing to shrinkage→"),
      name   = NULL
    )
  ) +
  guides(
    fill  = guide_legend(nrow = 2, ncol = 2, byrow = TRUE, override.aes = list(shape = c(23, NA, NA, NA))),
    shape = "none"
  ) +
  labs(
    x     = "Hedges' g",
    y     = NULL,
    fill  = NULL,
    title = NULL
  ) +
  theme_prism() +
  theme(
    legend.position      = "bottom",
    panel.grid.major.x   = element_line(color = "gray90"),
    panel.grid.major.y   = element_blank(),
    axis.ticks.x.top     = element_blank(),
    axis.ticks.y         = element_blank(),
    panel.border         = element_blank(),
    axis.line.x.top      = element_blank(),
    axis.line.y.right    = element_blank(),
    axis.text.x.top      = element_text(size = 9, margin = margin(b = 4))
  )

p_shapley

saveRDS(p_shapley, file.path(panels_dir, "p_shapley.rds"))
        
ggsave(
  filename = file.path(save_dir_decide, "shrinkage_shapley.png"),
  plot = p_shapley, width = 10, height = 7, dpi = 300
)

# Summary table (Shapley)
decomp_summary <- decomp_dt[, .(
  Project             = project_letter,
  `|g| (exp)`         = round(abs(g_exp), 3),
  `|g| (conf)`        = round(abs(g_conf), 3),
  `Total shrinkage`   = round(total_shrinkage, 3),
  `φ ctrl stability`      = round(phi_ctrl, 3),
  `φ treatment`       = round(phi_treated, 3),
  `φ variance`        = round(phi_sd, 3)
)]

setorder(decomp_summary, Project)
print(decomp_summary)

# Plot total shrinkage along with Shapley values
p_total <- ggplot(decomp_dt, aes(x = total_shrinkage, 
                                 y = factor(project_letter, levels = project_ord))) +
  geom_col(fill = "#4C9BE8", width = 0.7) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  labs(
    x     = "|exploratory g| − |confirmatory g|",
    y     = NULL,
    title = NULL
  ) +
  theme_prism() +
  theme(
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank()
  )

p_comb_decide <- p_total + p_shapley + plot_layout(ncol = 2)
p_comb_decide

ggsave(
  filename = file.path(save_dir_decide, "shrinkage_shapley.png"),
  plot = p_comb_decide, width = 10, height = 7, dpi = 300
)

##### Compare absolute vs. sign-projected Shapley decomposition ####

# Recompute using decomp_dt signed coalition g-values
# decomp_dt currently stores g0, g_c, etc as absolute values
# Take the signed versions and rebuild them the same way they were originally derived, before the abs() was applied.

decomp_compare <- decomp_dt[, .(
  project_letter,
  g_exp, g_conf,
  
  # signed coalition values
  g0_signed    = g_exp,
  g_c_signed   = hedges_g_from_parts(mu_ctrl_conf,  mu_treated_exp,  sd_exp_pooled,  n_ctrl_exp,  n_treated_exp),
  g_t_signed   = hedges_g_from_parts(mu_ctrl_exp,   mu_treated_conf, sd_exp_pooled,  n_ctrl_exp,  n_treated_exp),
  g_s_signed   = hedges_g_from_parts(mu_ctrl_exp,   mu_treated_exp,  sd_conf_pooled, n_ctrl_conf, n_treated_conf),
  g_ct_signed  = hedges_g_from_parts(mu_ctrl_conf,  mu_treated_conf, sd_exp_pooled,  n_ctrl_exp,  n_treated_exp),
  g_cs_signed  = hedges_g_from_parts(mu_ctrl_conf,  mu_treated_exp,  sd_conf_pooled, n_ctrl_conf, n_treated_conf),
  g_ts_signed  = hedges_g_from_parts(mu_ctrl_exp,   mu_treated_conf, sd_conf_pooled, n_ctrl_conf, n_treated_conf),
  g_cts_signed = g_conf
)]

# Compute both versions row by row
results_list <- lapply(seq_len(nrow(decomp_compare)), function(i) {
  row <- decomp_compare[i]
  s <- sign(row$g0_signed)  # sign anchor = exploratory direction
  
  # ORIGINAL: v(S) = |g(S)|
  orig <- compute_shapley(
    abs(row$g0_signed), abs(row$g_c_signed), abs(row$g_t_signed), abs(row$g_s_signed),
    abs(row$g_ct_signed), abs(row$g_cs_signed), abs(row$g_ts_signed), abs(row$g_cts_signed)
  )
  total_orig <- abs(row$g0_signed) - abs(row$g_cts_signed)
  
  # FIX: v(S) = sign(g_exp) * g(S)
  fixed <- compute_shapley(
    s * row$g0_signed, s * row$g_c_signed, s * row$g_t_signed, s * row$g_s_signed,
    s * row$g_ct_signed, s * row$g_cs_signed, s * row$g_ts_signed, s * row$g_cts_signed
  )
  # Sign-projected shrinkage: |g_exp| - sign(g_exp)*g_conf
  # Equals |g_exp| - |g_conf| when signs agree, but correctly detects full
  # reversals as maximal shrinkage instead of scoring them as zero/negative like the older |g_exp| - |g_conf|
  total_fixed <- s * row$g0_signed - s * row$g_cts_signed
  
  data.table(
    project_letter    = row$project_letter,
    sign_reversal      = sign(row$g0_signed) != sign(row$g_cts_signed),
    
    total_shrinkage_orig  = round(total_orig, 3),
    total_shrinkage_fixed = round(total_fixed, 3),
    
    phi_ctrl_orig     = round(orig$phi_ctrl, 3),
    phi_ctrl_fixed    = round(fixed$phi_ctrl, 3),
    
    phi_treated_orig  = round(orig$phi_treated, 3),
    phi_treated_fixed = round(fixed$phi_treated, 3),
    
    phi_sd_orig       = round(orig$phi_sd, 3),
    phi_sd_fixed      = round(fixed$phi_sd, 3)
  )
})

comparison_dt <- rbindlist(results_list)

# Check that additivity holds under both formulas
comparison_dt[, check_orig  := round(phi_ctrl_orig  + phi_treated_orig + phi_sd_orig - total_shrinkage_orig, 6)]
comparison_dt[, check_fixed := round(phi_ctrl_fixed + phi_treated_fixed + phi_sd_fixed - total_shrinkage_fixed, 6)]

print(comparison_dt)

# Flag which projects change meaningfully (not just full reversals)
comparison_dt[, total_changed := total_shrinkage_orig != total_shrinkage_fixed]
comparison_dt[, attribution_changed := 
                (phi_ctrl_orig != phi_ctrl_fixed) | 
                (phi_treated_orig != phi_treated_fixed) | 
                (phi_sd_orig != phi_sd_fixed)
]

comparison_dt[, .(project_letter, sign_reversal, total_changed, attribution_changed)]

##### Signed Shapley plots ####

shapley_long_signed <- melt(
  comparison_dt[, .(
    project_letter,
    `Control stability`  = phi_ctrl_fixed,
    `Treatment response` = phi_treated_fixed,
    `Variance inflation` = phi_sd_fixed
  )],
  id.vars = "project_letter", variable.name = "component", value.name = "phi"
)
setDT(shapley_long_signed)
set(shapley_long_signed, j = "project_letter", value = factor(shapley_long_signed$project_letter, levels = project_ord))
set(shapley_long_signed, j = "component", value = factor(shapley_long_signed$component,
                                                         levels = c("Control stability", "Treatment response", "Variance inflation")))

p_shapley_signed <- ggplot(shapley_long_signed, aes(x = phi, y = project_letter, fill = component)) +
  geom_col(position = "stack", width = 0.7) +
  geom_point(
    data = comparison_dt,
    aes(x = total_shrinkage_fixed, 
        y = factor(as.character(project_letter), levels = levels(shapley_long_signed$project_letter)), 
        fill = "Total shrinkage", shape = "Total shrinkage"),
    inherit.aes = FALSE,
    size = 2.5,
    color = "black"
  ) +
  scale_shape_manual(name = NULL, 
                     values = c("Total shrinkage" = 23),
                     labels = c("Total shrinkage" = "Total Shrinkage")
  ) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  scale_fill_manual(
    values = c(component_colors, "Total shrinkage" = "black"),
    breaks = c("Total shrinkage", names(component_colors)),
    labels = c(
      "Total shrinkage"     = "Total shrinkage",
      "Control stability"   = "Control Stability",
      "Treatment response"  = "Treatment Response",
      "Variance inflation"  = "Variance Inflation"
    )
  ) +
  scale_x_continuous(
    sec.axis = dup_axis(
      breaks = c(-2, 4),
      labels = c("← reducing", "contributing to shrinkage→"),
      name   = NULL
    )
  ) +
  guides(
    fill  = guide_legend(nrow = 2, ncol = 2, byrow = TRUE, override.aes = list(shape = c(23, NA, NA, NA))),
    shape = "none"
  ) +
  labs(
    x     = "Hedges' g",
    y     = NULL,
    fill  = NULL,
    title = NULL
  ) +
  theme_prism() +
  theme(
    legend.position      = "bottom",
    panel.grid.major.x   = element_line(color = "gray90"),
    panel.grid.major.y   = element_blank(),
    axis.ticks.x.top     = element_blank(),
    axis.ticks.y         = element_blank(),
    panel.border         = element_blank(),
    axis.line.x.top      = element_blank(),
    axis.line.y.right    = element_blank(),
    axis.text.x.top      = element_text(size = 9, margin = margin(b = 4))
  )

p_shapley_signed

p_total_signed <- ggplot(comparison_dt, aes(x = total_shrinkage_fixed, 
                                            y = factor(project_letter, levels = project_ord))) +
  geom_col(fill = "#4C9BE8", width = 0.7) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  labs(
    x     = "|exploratory g| − sign(exploratory g)·confirmatory g",
    y     = NULL,
    title = NULL
  ) +
  theme_prism() +
  theme(
    panel.grid.major.x = element_line(color = "gray90"),
    panel.grid.major.y = element_blank()
  )

p_comb_decide_signed <- p_total_signed + p_shapley_signed + plot_layout(ncol = 2)
p_comb_decide_signed

# Save under DIFFERENT names — does not overwrite original p_shapley/p_total
saveRDS(p_shapley_signed, file.path(panels_dir, "p_shapley_signed.rds"))

ggsave(
  filename = file.path(save_dir_decide, "shrinkage_shapley_signed.png"),
  plot = p_comb_decide_signed, width = 10, height = 7, dpi = 300
)
