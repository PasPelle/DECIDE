# INTRO Replication assessment of the retrospective dataset and comparison plots 
# with the confirmatory dataset
library(here)
source(here("packages.R"))

# Run functions script
source(here("all_functions.R"))

# Import DECIDE, retrospective multi-lab
decide_df        <- read_excel(here("results", "all_stats_wide.xlsx"))
retrospective_df <- fread(here("data", "mastersheet_external_studies.csv"), dec = ",")
decide_iv        <- fread(here("data", "decide_mIV_anonymized.csv")) # This file contains the minimal IV score used for comparison with retrospective dt

save_dir_external <- here("results", "retrospective_dt_and_comparisons")
save_dir_decide <- here("results")
panels_dir <- file.path(save_dir_decide, "panels")             
if (!dir.exists(panels_dir)) dir.create(panels_dir, recursive = TRUE)

# Load
replication_long_decide <- setDT(readRDS(here("results", "replication_long_decide.rds")))

#### Data preparation #####
# Reshape datasets for next analyses

## retrospective multi-lab
retrospective_dt <- as.data.table(retrospective_df)

# Leave only experimental groups of the Hedge's g 
retrospective_castable <- retrospective_dt[!is.na(group_es)]

# Remove species that don't match with the exploratory (project 9)
retrospective_dt <- retrospective_dt[!(id == 9 & subgroup %in% c("Mouse", "Rabbit"))]

# Exclude projects #4 and #19 (not eligible after internal discussions)
retrospective_dt <- retrospective_dt[!id %in% c(4, 19)]

# create unique lab IDs
retrospective_castable[, unique_lab := paste(id, center, sep = "_")]

retrospective_castable <- retrospective_castable[group_es %in% c("ctrl", "treated")]

# dcast to get only 1 row per lab
retrospective_dt_wide <- dcast(retrospective_castable,
                          unique_lab + id + stage + center + subgroup ~  group_es,
                          value.var = "sample_size"
                          )
setDT(retrospective_dt_wide)
# add the numeric columns back to the wide dt
ext_num_dt <- unique(
  retrospective_castable[
    ,
    .(
      unique_lab, stage, center, subgroup,
      mean_difference, pooled_sd, J, hedges_g, se_g, EU, iv_score
    )
  ],
  by = c("unique_lab", "stage", "center", "subgroup")
)

retrospective_dt_wide <- ext_num_dt[
  retrospective_dt_wide,
  on = .(unique_lab, stage, center, subgroup)
]

# Variance Analysis 1 _________________________####

##### Coefficient of Variation ####
multi_lab_dt <- retrospective_dt_wide[stage == "Multi_lab"]

## Extract CV from retrospective studies
# multi_lab_dt[, hedges_g := as.numeric(as.character("hedges_g"))]


cv_ext_by_project <- multi_lab_dt[
  ,
  .(mean_g = mean(hedges_g),
    sd_g = sd(hedges_g),
    cv_g = sd(hedges_g) / abs(mean(hedges_g)),
    dataset = "retrospective"
    ),
  by = id
]

# Extract CV from DECIDE confirmatory studies
decide_dt <- as.data.table(decide_df)
confirmatory_decide <- decide_dt[Stage == "confirmatory",]

cv_decide_by_project <- confirmatory_decide[
  ,
  .(mean_g = mean(hedges_g),
    sd_g = sd(hedges_g),
    cv_g = sd(hedges_g) / abs(mean(hedges_g)),
    dataset = "DECIDE"
  ),
  by = project_letter]


# Bind the 2 datasets for plotting
setnames(cv_decide_by_project, "project_letter", "id")
ext_decide <- rbind(cv_ext_by_project, cv_decide_by_project)

cv_decide_vs_retrospective <- ggplot(ext_decide, 
                                     aes(x = factor(dataset, levels = c("retrospective", "DECIDE")), 
                                         y = cv_g, fill = dataset)) +
  geom_boxplot(alpha = 0.85, width = 0.5, outlier.shape = NA, linewidth = 0.6) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) + 
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(
    title = NULL,
    x = NULL,
    y = "Coefficient of Variation (CV)"
  ) +
  scale_x_discrete(
    labels = c(
      retrospective = "Retrospective\nMulti-lab",
      DECIDE = "Confirmatory\nMulti-lab"
    )
  ) +
  scale_y_log10(
    breaks = c(0.01, 0.05, 0.1, 0.5, 1, 10),
    labels = c("0.01", "0.05", "0.1", "0.5", "1", "10")
  ) +
  theme_minimal() +
  theme_prism(
    axis_text_angle = 45) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5)
  )

cv_decide_vs_retrospective

saveRDS(cv_decide_vs_retrospective, file.path(panels_dir, "cv_decide_vs_retrospective.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "cv_decide_vs_retrospective.png"),
  plot = cv_decide_vs_retrospective,
  width = 10,
  height = 6,
  dpi = 300
)

stat_table_cv   <- make_pw_wilcox_stat_tbl(data = ext_decide, 
                                           y = "cv_g",  
                                           group = "dataset")


##### Root Mean Square Error (RMSE) ####
rmse_ext <- multi_lab_dt[,
                           {
                             g_bar <- mean(hedges_g)
                             list(
                                 k = .N,
                                 g_bar = g_bar,
                                 rmse = sqrt(mean((hedges_g - g_bar)^2)),
                                 dataset = "retrospective"
                             )
                           },
                           by = id
                         ]
                         

rmse_decide <- confirmatory_decide[,
                         {
                           g_bar <- mean(hedges_g)
                           list(
                             k = .N,
                             g_bar = g_bar,
                             rmse = sqrt(mean((hedges_g - g_bar)^2)),
                             dataset = "DECIDE"
                           )
                         },
                         by = project_letter
                         ]


# Bind the 2 datasets for plotting
setnames(rmse_decide, "project_letter", "id")
rmse_ext_decide <- rbind(rmse_ext, rmse_decide)

rmse_ext_decide_plot <- ggplot(rmse_ext_decide, 
                               aes(x = factor(dataset, levels = c("retrospective", "DECIDE")), 
                                   y = rmse, fill = dataset)) +
  geom_boxplot(alpha = 0.85, width = 0.5, outlier.shape = NA, linewidth = 0.6) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  labs(
    title = NULL,
    x = NULL,
    y = "Root Mean Squared Error (RMSE)"
  ) +
  scale_x_discrete(
    labels = c(
      retrospective = "Multi-lab \nRetrospective",
      DECIDE = "Multi-lab \nConfirmatory"
    )
  ) +
  # scale_y_log10(
  #   breaks = c(0.01, 0.05, 0.1, 0.5, 1, 10),
  #   labels = c("0.01", "0.05", "0.1", "0.5", "1", "10")
  # ) +
  theme_minimal() +
  theme_prism(
    axis_text_angle = 45
  ) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
rmse_ext_decide_plot

saveRDS(rmse_ext_decide_plot, file.path(panels_dir, "rmse_ext_decide_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "rmse_ext_decide_plot.png"),
  plot = rmse_ext_decide_plot,
  width = 10,
  height = 6,
  dpi = 300
)

stat_table_rmse <- make_pw_wilcox_stat_tbl(data = rmse_ext_decide,  y = "rmse",  group = "dataset")


# Meta-analysis_______________ ####


## Confirmatory DECIDE ##
meta_results_decide <- confirmatory_decide[
  ,
  {
    res <- run_meta_fixed(hedges_g, se_g, Center)
    c(extract_meta_fixed(res),list(dataset = "DECIDE"))
  },
  by = .(project_letter)
]

## retrospective multi-lab
meta_results_ext <- multi_lab_dt[
  ,
  {
    res <- run_meta_fixed(hedges_g, se_g, center)
    c(extract_meta_fixed(res),list(dataset = "retrospective"))
  },
  by = id
]

# Mapping id for plotting
id_plotid_map <- data.table(
  id      = c(1, 6, 9,  11, 12, 13, 16, 17, 18),
  plot_id = c(1, 2, 3,  4,  5,  6,  7,  8,  9)
)
id_plotid_map[, plot_id := as.factor(plot_id)]

# Add plot_id to meta_results_ext
meta_results_ext[, id := as.numeric(as.character(id))]
meta_results_ext <- merge(meta_results_ext, id_plotid_map, by = "id", all.x = TRUE)

#### Meta-analysis summary tables

k_decide <- confirmatory_decide[, .N, by = project_letter]
setnames(k_decide, "N", "k")
setnames(k_decide, "project_letter", "id")

meta_summary_decide <- meta_results_decide[, id := project_letter][
  project_stats[, .(id = project_letter, I2, tau2)], on = "id"
][
  k_decide, on = "id"
]

meta_summary_decide[, .(
  Project    = id,
  k,
  g          = round(g_pooled, 3),
  SE         = round(se, 3),
  `CI lower` = round(ci_lower, 3),
  `CI upper` = round(ci_upper, 3),
  p          = pval,
  `I² (%)`   = i.I2,
  `τ²`       = tau2
)][order(Project)] |>
  gt() |>
  tab_header(title = "Meta-analysis summary — DECIDE confirmatory dataset") |>
  tab_spanner(
    label   = "Pooled estimate",
    columns = c("g", "SE", "CI lower", "CI upper", "p")
  ) |>
  tab_spanner(
    label   = "Between-lab variability (descriptive)",
    columns = c("I² (%)", "τ²")
  ) |>
  fmt_scientific(columns = "p", decimals = 2) |>
  cols_align(align = "center") |>
  tab_style(
    style     = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) |>
  tab_footnote(
    footnote = paste(
      "Pooled estimates from fixed-effects meta-analysis.",
      "I² and τ² are descriptive measures of between-lab variability",
      "and should not be interpreted as formal heterogeneity tests",
      "given the small number of laboratories per project (mostly k = 2)."
    ),
    locations = cells_column_spanners(spanners = "Between-lab variability (descriptive)")
  )
meta_summary_decide_export <- meta_summary_decide[, .(
  Project    = id,
  k,
  g          = round(g_pooled, 3),
  SE         = round(se, 3),
  `CI lower` = round(ci_lower, 3),
  `CI upper` = round(ci_upper, 3),
  p          = pval,
  `I2 (%)`   = i.I2,
  `tau2`     = tau2
)][order(Project)]

# Get exploratory g and se per project for DECIDE
exploratory_decide_es <- decide_dt[Stage == "exploratory", .(
  project_letter,
  g_exploratory  = round(hedges_g, 3),
  se_exploratory = round(se_g, 3)
)]

meta_summary_decide_export <- meta_summary_decide[
  exploratory_decide_es[, .(id = project_letter, g_exploratory, se_exploratory)],
  on = "id"
][, .(
  Project        = id,
  k,
  `g exploratory`  = g_exploratory,
  `SE exploratory`  = se_exploratory,
  `g multi-lab`  = round(g_pooled, 3),
  `SE multi-lab` = round(se, 3),
  `95% CI`       = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
  p              = pval,
  `I2 (%)`       = i.I2,
  `tau2`         = tau2
)][order(Project)]

write.xlsx(
  meta_summary_decide_export,
  file.path(save_dir_decide, "meta_summary_decide.xlsx"),
  rowNames = FALSE
)

# Fixed-effect meta-analysis summary table - retrospective
k_ext <- multi_lab_dt[, .N, by = id]
setnames(k_ext, "N", "k")
k_ext <- id_plotid_map[k_ext, on = "id"]

meta_summary_ext <- meta_results_ext[
  k_ext[, .(plot_id, k)]
  , on = "plot_id"
]

meta_summary_ext_export <- meta_summary_ext[!is.na(plot_id), .(
  Project         = plot_id,
  k,
  `g multi-lab`   = round(g_pooled, 3),
  `SE multi-lab`  = round(se, 3),
  `95% CI`        = paste0("[", round(ci_lower, 3), ", ", round(ci_upper, 3), "]"),
  p               = pval,
  `I2 (%)`        = round(as.numeric(I2) * 100, 1),
  `tau2`          = round(tau2, 4)
)][order(Project)]


# Get exploratory g and se per project for retrospective
retrospective_exploratory <- retrospective_dt_wide[stage == "Exploratory", ]

exploratory_ext <- retrospective_exploratory[, .(
  id = id,
  `g exploratory` = round(hedges_g, 3),
  `SE exploratory` = round(se_g, 3)
)]

# Join plot_id
exploratory_ext <- id_plotid_map[exploratory_ext, on = "id"]
setnames(exploratory_ext, "plot_id", "Project")
exploratory_ext[, Project := as.factor(Project)]

meta_summary_ext_export <- meta_summary_ext_export[
  exploratory_ext[, .(Project, `g exploratory`, `SE exploratory`)], 
  on = "Project"
]

# Reorder columns
setcolorder(meta_summary_ext_export, c(
  "Project", "k",
  "g exploratory", "SE exploratory",
  "g multi-lab", "SE multi-lab", "95% CI", "p",
  "I2 (%)", "tau2"
))

print(meta_summary_ext_export)

write.xlsx(
  meta_summary_ext_export,
  file.path(save_dir_external, "meta_summary_retrospective.xlsx"),
  rowNames = FALSE
)

##### Heterogeneity plots - decide confirmatory vs retrospective multi-lab ####

# Prepare numeric I2 and in % and ordered factor using plot_id
meta_results_ext[, I2_numeric := as.numeric(I2 *100)]
meta_results_ext[, plot_id := factor(plot_id, levels = rev(sort(unique(plot_id))))]

# Plot
heterogeneity_plot_ext <- ggplot(meta_results_ext, aes(x = I2_numeric, y = plot_id)) +
  geom_bar(
    stat = "identity", 
    fill = "#E7B800",
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
    title = "retrospective"
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )

heterogeneity_plot_ext

# Plot I2 heterogeneity in % in confirmatory projects
meta_results_decide[, I2_numeric := as.numeric(I2 * 100)]
meta_results_decide[, project_letter := factor(project_letter, levels = rev(sort(unique(project_letter))))]

heterogeneity_plot_decide <- ggplot(meta_results_decide, aes(x = I2_numeric, y = project_letter)) +
  geom_bar(
    stat = "identity",
    fill = "#00AFBB",
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
    title = "confirmatory"
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 10),
    axis.title.x = element_text(margin = margin(t = 10), size = 11),
    plot.margin = margin(10, 10, 10, 10)
  )
heterogeneity_plot_decide

heterogeneity_plot_ext / heterogeneity_plot_decide

# IV scores_______________________________####

# Merge decide IV scores (s100a has no confirmatory data so no letter was assigned)
decide_iv <- decide_iv[!is.na(IV_final)]
decide_iv[, dataset := "DECIDE"]

setnames(decide_iv, "Project_ID", "unique_lab")
setnames(decide_iv, "IV_final", "iv_score")
setnames(decide_iv, "Study_Type", "stage")
decide_iv[trimws(unique_lab) == "", unique_lab := "K"]

# retrospective: removed project 12 bc it was used for 2 mutli-lab outcomes coming from the same study (Llovera 2015) so I have deleted it to not double count the same IV
retrospective_iv <- retrospective_dt[!is.na(iv_score) & id != 12,
                           .(id, 
                             stage, 
                             iv_score, 
                             dataset = "retrospective")]
retrospective_iv <- id_plotid_map[retrospective_iv, on = "id"]

# Bind and plot IV scores
combined_iv <- rbindlist(
  list(
    decide_iv[, .(project_id = unique_lab, stage, dataset, iv_score)],
    retrospective_iv[, .(project_id = as.character(plot_id), stage, dataset, iv_score)]  
    ),
  use.names = TRUE,
  fill = TRUE
)
combined_iv[, iv_score := as.numeric(iv_score)]

# Build combined group variable
combined_iv$group <- interaction(
  combined_iv$dataset,
  combined_iv$stage,
  sep = " - "
)

# Set the order of levels exactly as they appear in the data
combined_iv$group <- factor(
  combined_iv$group,
  levels = c(
    "retrospective - Exploratory",
    "retrospective - Multi_lab",
    "DECIDE - Exploratory",
    "DECIDE - Confirmatory"
  )
)

# Plot only retrospective and DECIDE
combined_iv_decide_retrospective <- combined_iv[dataset != "pct"]

# Stats
kruskal.test(iv_score ~ group, data = combined_iv_decide_retrospective)

# Pairwise Wilcoxon with BH correction
# pw <- pairwise.wilcox.test(combined_iv_decide_retrospective$iv_score, 
#                            combined_iv_decide_retrospective$group, 
#                            p.adjust.method = "BH")
# 
# # Convert p-value matrix to long table of comparisons
# p_long <- as.data.table(as.table(pw$p.value))
# setnames(p_long, c("group2", "group1", "p_adj"))  # note: matrix is lower triangle
# p_long <- p_long[!is.na(p_adj)]
# 
# # Keep only significant comparisons (BH-adjusted)
# p_sig <- p_long[p_adj < 0.05]
# 
# # Nice label format for plotting
# p_sig[, p_label := ifelse(p_adj < 0.001, "P < 0.001",
#                           paste0("P = ", formatC(p_adj, format = "f", digits = 3)))]
# 
# # order comparisons so  "shorter" comparisons are lower 
# # This uses the factor order of group on the x-axis.
# lvl <- levels(combined_iv_decide_retrospective$group)
# p_sig[, i1 := match(group1, lvl)]
# p_sig[, i2 := match(group2, lvl)]
# p_sig[, span := abs(i2 - i1)]
# setorder(p_sig, span, i1, i2)   # shorter brackets first
# 
# # y positions with sufficient headroom 
# y_rng  <- range(combined_iv_decide_retrospective$iv_score, na.rm = TRUE)
# y_span <- diff(y_rng)
# ymax   <- y_rng[2]
# 
# step <- 0.08 * y_span   # key change: larger step
# p_sig[, y.position := ymax + step * seq_len(.N)]
# 
# # Table for ggpubr (same name as before)
# stat_table_iv <- p_sig[, .(group1, group2, p.adj = p_adj, p = p_adj,
#                       p_label, y.position)]

combined_iv_decide_retrospective[stage == "Confirmatory", stage := "Multi_lab"]

write.csv(combined_iv_decide_retrospective, 
          here("results", "combined_iv_decide_retrospective.csv"),
          row.names = FALSE)

combined_iv_decide_retrospective[
  ,
  stage := factor(stage, levels = c("Exploratory", "Multi_lab"))
]
combined_iv_decide_retrospective[
  ,
  dataset := factor(dataset, levels = c("DECIDE", "retrospective"))
]

stat_table_iv <- make_pw_wilcox_stat_tbl(
  data = combined_iv_decide_retrospective,
  y = "iv_score",
  group = "group"
)

stat_table_iv_selected <- stat_table_iv[
  group1 == "retrospective - Exploratory" & group2 == "retrospective - Multi_lab" |
  group1 == "DECIDE - Exploratory" & group2 == "DECIDE - Confirmatory" |
  group1 == "retrospective - Multi_lab" & group2 == "DECIDE - Confirmatory"
  
]

# 2x2 ANOVA
anova_result_mIV <- aov(iv_score ~ stage * dataset , data = combined_iv_decide_retrospective)
summary(anova_result_mIV)



# ART requires factors as factor(), and no missing data in the model vars
combined_iv_decide_retrospective[, `:=`(
  stage = factor(stage),
  dataset = factor(dataset)
)]

m_art <- art(iv_score ~ stage * dataset, data = combined_iv_decide_retrospective)

# Check ART alignment worked correctly (should show F ≈ 0 for the "wrong" effects per row)
summary(m_art)

# Omnibus ANOVA on aligned ranks — main effects + interaction
anova(m_art)

combined_iv_decide_retrospective[, iv_score_ord := factor(iv_score, ordered = TRUE)]

m_clm <- clm(iv_score_ord ~ stage * dataset, data = combined_iv_decide_retrospective)
summary(m_clm)

# Type II/III test for main effects + interaction (like an ANOVA table)

Anova(m_clm, type = "III")


cols_dataset <- c("DECIDE" = "#00AFBB", "retrospective" = "#E7B800")

combined_iv_decide_retrospective_plot <-
  ggplot(
    combined_iv_decide_retrospective,
    aes(
      x = group,
      y = iv_score,
      fill = dataset
    )
  ) +
  geom_boxplot(
    aes(alpha = stage),
    outlier.shape = NA,
    width = 0.55,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) +
  scale_fill_manual(values = cols_dataset) +
  scale_alpha_manual(
    values = c("Exploratory" = 0.35, "Multi_lab" = 0.85)
  ) +
  guides(fill = "none", alpha = "none") +
  scale_x_discrete(
    labels = c(
      "retrospective - Exploratory" = "Exploratory",
      "retrospective - Multi_lab"   = "Multi-lab",
      "DECIDE - Exploratory"        = "Exploratory",
      "DECIDE - Confirmatory"       = "Multi-lab"
    )
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "minimal Internal Validity (mIV)"
  ) +
  ggpubr::stat_pvalue_manual(
    stat_table_iv_selected,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(clip = "off") +
  annotate(
    "text", x = 1.5, y = 0, label = "Retrospective",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  annotate(
    "text", x = 3.5, y = 0, label = "Confirmatory",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  theme_prism() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(vjust = 0.5),
    axis.title.y = element_text(size = 16),
    plot.margin = margin(10, 40, 50, 10)
  )

combined_iv_decide_retrospective_plot

saveRDS(combined_iv_decide_retrospective_plot, file.path(panels_dir, "combined_iv_decide_retrospective_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "combined_iv_decide_retrospective.png"),
  plot = combined_iv_decide_retrospective_plot,
  width = 10,
  height = 6,
  dpi = 300
)

#### Bayesian analysis of iv data from confirmatory studies

# read in data for confirmatory studies 
combined_iv_decide_retrospective$iv_score <- ordered(combined_iv_decide_retrospective$iv_score)
my_priors <- c(
  prior(normal(0, 1.5), class = "b"),
  prior(normal(0, 1.5), class = "Intercept")
)

m.1 <- brm(
  formula = iv_score ~  stage*dataset+(1|project_id),
  data = combi,
  ned_iv_decide_retrospective,
  prior = my_priors,
  family = cumulative ("probit"),
  chains = 4,          
  cores = 4,           
  iter = 2000,         
  seed = 42,
  save_pars = save_pars(all = TRUE)
)

summary(m.1)

# 1. Create a conditions dataframe for all levels of 'dataset'
conds <- make_conditions(m.1, vars = "dataset")

# 2. Generate the conditional effects plot
conditional_effects(
  m.1, 
  effects = "stage", 
  conditions = conds, 
  categorical = TRUE
)

# Plots actual data (y) vs. simulated data from the model (yrep)
pp_check(m.1, type = "bars", ndraws = 100) +
  ggplot2::labs(title = "Posterior Predictive Check: Actual vs. Predicted")

# Add LOO criterion to the model
m.1 <- add_criterion(m.1, "loo", moment_match = TRUE, overwrite = TRUE)

# Print the LOO results
loo(m.1)


# Experimental units ______________________####
# EU in exploratory vs multi-lab in DECIDE vs retrospective

retrospective_eu <- multi_lab_dt[,.(EU, dataset = "retrospective")]
decide_eu <- confirmatory_decide[,
                                 .(EU = n1 + n2,
                                   dataset = "DECIDE")]

decide_eu   <- decide_eu[!is.na(EU)]
retrospective_eu <- retrospective_eu[!is.na(EU)]



# bind DECIDE df
comb_decide_exernal_eu <- rbindlist(list(
  exploratory_eu_decide[,
                 .(
                    group = "exploratory DECIDE",
                    dataset = "DECIDE",
                    stage = "exploratory",
                    eu = exploratory_n 
                  )
                 ],
  confirmatory_eu_decide[,
                  .(
                    group = "confirmatory DECIDE",
                    dataset = "DECIDE",
                    stage = "Multi_lab",
                    eu = confirmatory_n
                  )
                  ],
  retrospective_exploratory[,
              .(
                group = "exploratory retrospective",
                dataset = "retrospective",
                stage = "exploratory",
                eu = EU
              )
              ],
  multi_lab_dt[,
               .(
                 group = "multi_lab retrospective",
                 dataset = "retrospective",
                 stage = "Multi_lab",
                 eu = EU
               )
               ]
  )
  )

comb_decide_exernal_eu$group <- factor(
  comb_decide_exernal_eu$group,
  levels = c(
    "exploratory retrospective",
    "multi_lab retrospective",
    "exploratory DECIDE",
    "confirmatory DECIDE"
  )
)

stat_tbl_eu <- make_pw_wilcox_stat_tbl(
  data = comb_decide_exernal_eu,
  y = "eu",
  group = "group"
)

stat_tbl_eu_selected <- stat_tbl_eu[
  group1 == "exploratory retrospective" & group2 == "multi_lab retrospective" |
    group1 == "exploratory DECIDE" & group2 == "confirmatory DECIDE" |
    group1 == "multi_lab retrospective" & group2 == "confirmatory DECIDE"
  
]


comb_decide_exernal_eu_plot <-
  ggplot(
    comb_decide_exernal_eu,
    aes(
      x = group,
      y = eu,
      fill = dataset
    )
  ) +
  geom_boxplot(
    aes(alpha = stage),
    outlier.shape = NA,
    width = 0.55,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  )+
  scale_fill_manual(values = cols_dataset) +
  scale_alpha_manual(
    values = c("exploratory" = 0.35, "Multi_lab" = 0.85)
  ) +
  scale_x_discrete(
    labels = c(
      "exploratory retrospective" = "Exploratory",
      "multi_lab retrospective"   = "Multi-lab",
      "exploratory DECIDE"        = "Exploratory",
      "confirmatory DECIDE"       = "Multi-lab"
    )
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Experimental Units",
    fill = NULL
  ) +
  ggpubr::stat_pvalue_manual(
    stat_tbl_eu_selected,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(clip = "off") +
  annotate(
    "text", x = 1.5, y = 0, label = "Retrospective",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  annotate(
    "text", x = 3.5, y = 0, label = "Confirmatory",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  theme_prism() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    plot.margin = margin(10, 40, 50, 10)
  )

comb_decide_exernal_eu_plot

saveRDS(comb_decide_exernal_eu_plot, file.path(panels_dir, "comb_decide_exernal_eu_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "eu_plot_retrospective_decide.png"),
  plot = comb_decide_exernal_eu_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Print out total n of EU summing all projects in Expl vs multi-lab
total_eu <- comb_decide_exernal_eu[,
                                   .(total_eu = sum(eu)),
                                   by = group]
total_eu


# Replication assessment_________________________####

##### Replication assessment for the retrospective dataset ####
## Prepare dataset with exploratory and confirmatory
setnames(meta_results_ext,"g_pooled", "g_pooled_confirmatory")
setnames(meta_results_ext,"se", "se_g_confirmatory")
setnames(meta_results_ext,"ci_lower", "ci_lower_confirmatory")
setnames(meta_results_ext,"ci_upper", "ci_upper_confirmatory")
setnames(meta_results_ext,"pval", "pval_confirmatory")
setnames(meta_results_ext,"I2", "I2_confirmatory")

## Exploratory
retrospective_exploratory <- retrospective_dt_wide[stage == "Exploratory"]

# Add pvalue from Hedge's g
retrospective_exploratory <- retrospective_exploratory[,
                                             z_expl := hedges_g / se_g
                                             ]
retrospective_exploratory <- retrospective_exploratory[,
                                             p_val_exploratory := 2 * stats::pnorm(-abs(z_expl))
                                             ]

setnames(retrospective_exploratory, old = c("hedges_g", "se_g"), new = c("hedges_g_exploratory","se_g_exploratory"))
setnames(retrospective_exploratory, old = c("ctrl", "treated"), new = c("n1_exploratory", "n2_exploratory"))

# add 95% CI
retrospective_exploratory[,
               `:=`(ci_lower_exploratory = hedges_g_exploratory - 1.96 * se_g_exploratory,
                    ci_upper_exploratory = hedges_g_exploratory + 1.96 * se_g_exploratory
               )]


# Merge exploratory and confirmatory
retrospective_all <- merge(
  retrospective_exploratory[, .(id, 
                           hedges_g_exploratory,
                           se_g_exploratory,
                           p_val_exploratory,
                           n1_exploratory,
                           n2_exploratory,
                           ci_lower_exploratory,
                           ci_upper_exploratory)],
  meta_results_ext,
  by = "id",
  all = TRUE
)


# Add more columns for the assessment: pooled Smallest Effect Size SDE
setnames(multi_lab_dt, old = c("ctrl", "treated"), new = c("n1_confirmatory", "n2_confirmatory"))
multi_lab_dt[
  n1_confirmatory >= 2 & n2_confirmatory >= 2,
  sde := mapply(
    function(n1, n2) {
      pwr::pwr.t2n.test(
        n1 = n1,
        n2 = n2,
        sig.level = 0.05,
        power = 0.8
      )$d
    },
    n1_confirmatory,
    n2_confirmatory
  )
]


# Convert d to g
multi_lab_dt[, sde_g := sde *(1- (3 / (4*(n1_confirmatory + n2_confirmatory) - 9)))]

# Compute SE for each SDE
multi_lab_dt[, se_sde_g := sqrt((n1_confirmatory + n2_confirmatory)/(n1_confirmatory * n2_confirmatory) + (sde_g^2 / (2 * (n1_confirmatory + n2_confirmatory))))]

# Meta-analyze SDEs across centers to get pooled SDE threshold
sde_meta <- multi_lab_dt[
  ,
  {
    m <- rma(
      yi = sde_g,
      sei = se_sde_g,
      method = "REML",
      data = .SD
    )
    list(
      SDE_pooled = m$beta[1],
      SDE_se = m$se,
      SDE_lower = m$ci.lb,
      SDE_upper = m$ci.ub
    )
  },
  by = id,
  .SDcols = c("sde_g", "se_sde_g")
]

# Add the SDE column
retrospective_all <- merge(
  retrospective_all,
  sde_meta[, .(id, SDE_pooled)],
  by = "id",
  all.x = TRUE
)

# === Run replication assessment function ===
retrospective_all <- compute_replication_flags(
  retrospective_all,
  cfg = list(alpha = 0.05, z_pi = 1.96, power_d33 = 0.33, power_sde = 0.8,
             sceptical_type = "golden", sceptical_alternative = "two.sided"),
  cols = list(
    id = "id",
    g_exp = "hedges_g_exploratory", se_exp = "se_g_exploratory",
    ci_lo_exp = "ci_lower_exploratory", ci_hi_exp = "ci_upper_exploratory", p_exp = "p_val_exploratory",
    n1_exp = "n1_exploratory", n2_exp = "n2_exploratory",
    g_conf = "g_pooled_confirmatory", se_conf = "se_g_confirmatory",
    ci_lo_conf = "ci_lower_confirmatory", ci_hi_conf = "ci_upper_confirmatory", p_conf = "pval_confirmatory",
    tau2 = "tau2",
    SDE_pooled = "SDE_pooled"
  ),
  compute_thresholds = TRUE
)


##### Replication matrix ####
replication_matrix_retrospective <- retrospective_all[, .(
  id,
  `Exploratory g in confirmatory CI` = ci_agreement,
  # `Exploratory g in confirmatory PI` = exploratory_within_confirmatory_PI,
  `Exploratory and confirmatory in the same direction` = direction_agreement,
  `Skeptical p < 0.05` = sceptical_sig,
  `Both p < 0.05 & Same Direction` = ttest_sig,
  `Small Telescopes` = small_telescope_confirmed,
  `Confirmatory larger than its SDE` = sde_confirmed
)]

replication_long_retrospective <- melt(replication_matrix_retrospective, 
                         id.vars = "id",
                         variable.name = "criterion", 
                         value.name = "success")
replication_long_retrospective <- as.data.table(replication_long_retrospective)

# Count TRUEs per project to rank them
replication_long_retrospective[, success_count := sum(success == TRUE, na.rm = TRUE), by = id]

#  Use project letter as factor and order by success count
# replication_long_retrospective[, project_letter := factor(as.character(project_letter), levels = sort(unique(as.character(project_letter))))]

replication_long_retrospective[, success_char := as.character(success)]

# # # Order projects by ID
# replication_long_retrospective[, id := factor(id)]
# 
# replication_long_retrospective[
#   order(id),
#   new_id := as.factor(.GRP),
#   by = id
# ]
replication_long_retrospective[, id := as.numeric(as.character(id))]
replication_long_retrospective <- merge(replication_long_retrospective, id_plotid_map, by = "id", all.x = TRUE)

# Heatmap
heatmap_plot_ext <- ggplot(replication_long_retrospective, aes(x = plot_id, y = criterion, fill = success_char)) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.85) +  # Added width/height < 1 to shrink tiles
  
  # More professional color scheme (Nature/Cell style)
  scale_fill_manual(
    values = c("TRUE" = "#1F6F70", "FALSE" = "#D7D7D7"),
    labels = c("TRUE" = "Met", "FALSE" = "Not met"),
    name = "Criterion status"
  ) +
    scale_y_discrete(
    drop = FALSE,
    limits = rev(c(
      "CI above SESOI",
      "Confirmatory larger than its SDE",
      "Small Telescopes",
      "Both p < 0.05 & Same Direction",
      "Skeptical p < 0.05",
      "Exploratory and confirmatory in the same direction",
      # "Exploratory g in confirmatory PI",
      "Exploratory g in confirmatory CI"
    )),
    labels = c(
      "CI above SESOI" = "CI above SESOI",
      "Confirmatory larger than its SDE" = "Confirmatory ES > mDES",
      "Small Telescopes" = "Small Telescopes",
      "Both p < 0.05 & Same Direction" = "Significant & Same Direction",
      "Skeptical p < 0.05" = "Skeptical p-value",
      "Exploratory and confirmatory in the same direction" = "Same Direction",
      # "Exploratory g in confirmatory PI" = "Exploratory ES within c-PI ",
      "Exploratory g in confirmatory CI" = "Exploratory ES within c-CI"
    )
  )+
  labs(
    x = "Retrospective Dataset",
    y = NULL,
    title = NULL
  ) +
  # Professional theme with LARGER text
  theme_minimal(base_family = "Arial", base_size = 14) +  
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


heatmap_plot_ext


# Save with adjusted dimensions for larger labels
ggsave(
  filename = file.path(save_dir_external, "replication_heatmap_retrospective.png"),
  plot = heatmap_plot_ext,
  width = 16,   
  height = 11,  
  dpi = 300,
  bg = "white"
)

# Plot together the 2 dataset heatmaps

# Shared y-axis order & labels
y_order <- rev(c(
  "CI above SESOI",
  "Confirmatory larger than its SDE",
  "Small Telescopes",
  "Both p < 0.05 & Same Direction",
  "Skeptical p < 0.05",
  "Exploratory and confirmatory in the same direction",
  # "Exploratory g in confirmatory PI",
  "Exploratory g in confirmatory CI"
))

# Relabel to improve readability
y_labels <- c(
  "CI above SESOI"                                       = "CI above SESOI",
  "Confirmatory larger than its SDE"                     = "Confirmatory ES > mDES",
  "Small Telescopes"                                     = "Small Telescopes",
  "Both p < 0.05 & Same Direction"                       = "Significant & Same Direction",
  "Skeptical p < 0.05"                                   = "Skeptical p-value",
  "Exploratory and confirmatory in the same direction"   = "Same Direction",
  # "Exploratory g in confirmatory PI"                     = "Exploratory ES within c-PI",
  "Exploratory g in confirmatory CI"                     = "Exploratory ES within c-CI"
)

# Shared theme base
shared_theme <- list(
  theme_prism(base_family = "Arial", base_size = 14),
  theme(
    axis.text.x     = element_text(size = 16, face = "bold", color = "black"),
    axis.title.x    = element_text(size = 15, face = "bold", margin = margin(b = 8)),
    legend.position = "none",        # handled once via patchwork
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background= element_rect(fill = "white", color = NA)
  )
)
replication_long_decide[, success_char := as.character(as.logical(success))]

# Right panel: Confirmatory
p_decide <- ggplot(
  replication_long_decide,
  aes(x = project_letter, y = criterion, fill = success_char)
) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.85) +
  scale_fill_manual(
    values   = c("TRUE" = "#00AFBB", "FALSE" = "#D7D7D7"),
    labels   = c("TRUE" = "Met (Confirmatory)", "FALSE" = "Not met"),
    name     = "Criterion status",
    na.value = "#D7D7D7"
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(drop = FALSE, limits = y_order, labels = y_labels) +
  labs(x = NULL, y = NULL) +
  ggtitle("Confirmatory Studies") +
  shared_theme +
  theme(
    axis.text.x  = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title   = element_text(size = 30, face = "bold", hjust = 0.5,
                                margin = margin(b = 8)),
    plot.margin  = margin(10, 10, 10, 5, "pt")
  )

replication_long_retrospective[, success_char := as.character(as.logical(success))]
# Left panel: Retrospective
p_retro <- ggplot(
  replication_long_retrospective,
  aes(x = factor(id), y = criterion, fill = success_char)
) +
  geom_tile(color = "white", linewidth = 0.8, width = 0.85, height = 0.85) +
  scale_fill_manual(
    values   = c("TRUE" = "#E7B800", "FALSE" = "#D7D7D7"),
    labels   = c("TRUE" = "Met (Retrospective)", "FALSE" = "Not met"),
    name     = "Criterion status",
    na.value = "#D7D7D7"
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(drop = FALSE, limits = y_order, labels = y_labels) +
  labs(x = NULL, y = NULL) +           # no bottom title
  ggtitle("Retrospective Projects") +  # title at top instead
  shared_theme +
  theme(
    axis.text.x  = element_text(size = 16, face = "bold", color = "black"),
    axis.text.y  = element_text(size = 16, hjust = 1, color = "black"),
    plot.title   = element_text(size = 30, face = "bold", hjust = 0.5,
                                margin = margin(b = 8)),
    plot.margin  = margin(10, 5, 10, 10, "pt")
  )


# Combine with patchwork
# Width ratio: confirmatory has 4 cols (A–D), retrospective has 9 → ratio ~ 4:9
combined_heatmap_plot <- p_retro + p_decide +
  plot_layout(
    widths = c(9, 10),   # adjust to match column counts
    guides = "collect"
  ) &
  theme(
    legend.position      = "bottom",
    legend.title         = element_text(size = 18, face = "bold"),
    legend.text          = element_text(size = 17),
    legend.key.size      = unit(0.6, "cm"),
    legend.box.spacing   = unit(0.3, "cm"),
    plot.title           = element_text(size = 32, face = "bold", hjust = 0.5,
                                        margin = margin(b = 8))
  )

combined_heatmap_plot

saveRDS(combined_heatmap_plot, file.path(panels_dir, "combined_heatmap_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "replication_heatmap_combined.png"),
  plot     = combined_heatmap_plot,
  width    = 22,
  height   = 10,
  dpi      = 300,
  bg       = "white"
)


# Effect Size_________________________________________####
##### SDE of pooled labs ####

# Compare the SDE of exploratory vs pooled confirmatory labs

## retrospective
# Sum EUs per project and compute multi-lab SDE
retrospective_multilab_sde_sum_n <- multi_lab_dt[,
                               .(
                                 n1_total = sum(n1_confirmatory),
                                 n2_total = sum(n2_confirmatory)
                               ),
                               by = id
                                ][,
                                  sde := mapply(function(n1, n2) {
                                    pwr::pwr.t2n.test(n1 = n1,
                                                      n2 = n2,
                                                      power = 0.8)$d
                                    },
                                    n1_total,
                                    n2_total
                                    )
                                  ][,
                                    stage := "multi lab"
                                    ]
# retrospective_all <- sde_by_project[retrospective_all, on = "id"]

# retrospective exploratory SDE 


retrospective_exploratory_sde <- retrospective_all[, .(
  id,
  n1_exploratory,
  n2_exploratory
  )
  ][,
    sde := mapply(function(n1, n2) {
              pwr::pwr.t2n.test(n1 = n1,
                                n2 = n2,
                                power = 0.8)$d
              },
              n1_exploratory,
              n2_exploratory)
              ][,
                stage :="exploratory"
              ]

setnames(retrospective_exploratory_sde, 
         old = c("n1_exploratory", "n2_exploratory"),
         new = c("n1_total","n2_total"))

retrospective_sde_sum_n <- rbind(
  retrospective_multilab_sde_sum_n,
  retrospective_exploratory_sde
  )
retrospective_sde_sum_n[,
                 dataset := "retrospective"]


## DECIDE exploratory and confirmatory SDE

# DECIDE sum EU and compute confirmatory SDE

decide_confirmatory_sde <- confirmatory_decide[,
                               .(
                                 n1_total = sum(n1),
                                 n2_total = sum(n2)
                               ),
                               by = project_letter
][,
  sde := mapply(function(n1, n2) {
    pwr::pwr.t2n.test(n1 = n1,
                      n2 = n2,
                      power = 0.8)$d
  },
  n1_total,
  n2_total
  )
][,
  stage := "confirmatory"         # adding back columns I need fir rbind
]

# DECIDE exploratory sde

exploratory_decide <- decide_dt[Stage == "exploratory"
                                ][,
                                  sde := mapply(function(n1, n2){
                                    pwr::pwr.t2n.test(n1 = n1,
                                                      n2 = n2,
                                                      power = 0.8)$d
                                  },
                                  n1,
                                  n2
                                  )
                                ]

setnames(exploratory_decide, 
         old = c("n1","n2", "Stage"),
         new = c("n1_total", "n2_total", "stage")
         )

decide_sde_sum_n <- rbind(
  decide_confirmatory_sde,
  exploratory_decide[,
                     .(project_letter,
                       stage,
                       n1_total,
                       n2_total,
                       sde)]
)
decide_sde_sum_n[,
                 dataset := "DECIDE"]
setnames(decide_sde_sum_n, "project_letter", "id")

combined_sde <- rbind(decide_sde_sum_n, retrospective_sde_sum_n)

# Build combined group variable
combined_sde$group <- interaction(
  combined_sde$dataset,
  combined_sde$stage,
  sep = " - "
)

# Set the order of levels exactly as they appear in the data
combined_sde$group <- factor(
  combined_sde$group,
  levels = c(
    "retrospective - exploratory",
    "retrospective - multi lab",
    "DECIDE - exploratory",
    "DECIDE - confirmatory"
  )
)

combined_sde[stage == "confirmatory", stage := "Multi-lab"]
combined_sde[stage == "exploratory", stage := "Exploratory"]
combined_sde[stage == "multi lab", stage := "Multi-lab"]


stat_table_sde <- make_pw_wilcox_stat_tbl(
  data = combined_sde,
  y = "sde",
  group = "group"
)

stat_table_sde_selected <- stat_table_sde[
  group1 == "retrospective - exploratory" & group2 == "retrospective - multi lab" |
    group1 == "DECIDE - exploratory" & group2 == "DECIDE - confirmatory" |
    group1 == "retrospective - multi lab" & group2 == "DECIDE - confirmatory"
]

combined_sde_decide_retrospective_plot <-
  ggplot(
    combined_sde,
    aes(
      x = group,
      y = sde,
      fill = dataset
    )
  ) +
  geom_boxplot(
    aes(alpha = stage),
    outlier.shape = NA,
    width = 0.55,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) +
  scale_fill_manual(values = cols_dataset) +
  scale_alpha_manual(
    values = c("Exploratory" = 0.35, "Multi_lab" = 0.85)
  ) +
  scale_x_discrete(
    labels = c(
      "retrospective - exploratory" = "Exploratory",
      "retrospective - multi lab"   = "Multi-lab",
      "DECIDE - exploratory"        = "Exploratory",
      "DECIDE - confirmatory"       = "Multi-lab"
    )
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Smallest Detectable Effect Size (SDE)",
    fill = NULL
  ) +
  ggpubr::stat_pvalue_manual(
    stat_table_sde_selected,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4
  ) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  coord_cartesian(clip = "off") +
  annotate(
    "text", x = 1.5, y = 0, label = "Retrospective",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  annotate(
    "text", x = 3.5, y = 0, label = "Confirmatory",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  theme_prism() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    plot.margin = margin(10, 40, 50, 10)
  )

combined_sde_decide_retrospective_plot

saveRDS(combined_sde_decide_retrospective_plot, file.path(panels_dir, "combined_sde_decide_retrospective_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "combined_sde_decide_retrospective_plot.png"),
  plot = combined_sde_decide_retrospective_plot,
  width = 10,
  height = 6,
  dpi = 300
)

kruskal.test(iv_score ~ group, data = combined_iv_decide_retrospective)

pairwise.wilcox.test(
  combined_iv_decide_retrospective$iv_score,
  combined_iv_decide_retrospective$group,
  p.adjust.method = "BH"   # Benjamini–Hochberg correction
)



##### ES by project retrospective multi-lab  ---------------------------------------

# Ensure correct factor levels so projects are ordered by ID
# retrospective_dt_wide[, Project_Label := factor(
#   paste0(Project_ID, "_", Project_Name),
#   levels = unique(paste0(Project_ID, "_", Project_Name))
# )]

# Ensure Stage is a factor with correct left-to-right order
retrospective_dt_wide[, stage := factor(stage, levels = c("Exploratory", "Multi_lab"))]

# Add 95% CI
retrospective_dt_wide[,
                 `:=`(
                   ci_lower = hedges_g - 1.96 * se_g,
                   ci_upper = hedges_g + 1.96 * se_g
                 )]

# Force base x position per project
retrospective_dt_wide[, base_x := as.numeric(factor(id))]

# Count per group to differentiate within Multi_lab
retrospective_dt_wide[, stage_index := seq_len(.N), by = .(id, stage)]
retrospective_dt_wide[, n_in_group := .N, by = .(id, stage)]
dx <- 0.12  # define horizontal spacing between labs

# Enhanced spacing logic
retrospective_dt_wide[, jittered_x := fifelse(
  stage == "Exploratory",
  base_x - 0.25,
  {
    # Multi_lab: centre the offsets so they don't drift too far right
    center_offset <- (stage_index - (n_in_group + 1) / 2) * dx
    base_x + 0.25 + center_offset
  }
)]


# retrospective_dt_wide[, id_f := factor(id, levels = sort(unique(id)))]

# Relabel for plotting
# retrospective_dt_wide[
#   order(id),
#   new_id := as.factor(.GRP),
#   by = id
# ]
retrospective_dt_wide[, id := as.numeric(as.character(id))]
retrospective_dt_wide <- merge(retrospective_dt_wide, id_plotid_map, by = "id", all.x = TRUE)


x_breaks <- unique(retrospective_dt_wide[, .(base_x, plot_id)])[order(base_x)]


hedges_by_project_retrospective_plot <- ggplot(
  retrospective_dt_wide,
  aes(x = jittered_x, y = hedges_g, color = stage, shape = stage)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.08
  ) +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, max(retrospective_dt_wide$base_x) - 0.5, by = 1),
    linetype = "solid", color = "gray85", linewidth = 0.4
  ) +
  scale_x_continuous(
    breaks = x_breaks$base_x,
    labels = x_breaks$plot_id,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_color_manual(
    values = c(
      Exploratory = scales::alpha("#CC3300", 0.9),
      Multi_lab   = "#E7B800"
    ),
    labels = c(
      Exploratory = "Exploratory",
      Multi_lab   = "Multi-lab"
    )
  ) +
  scale_shape_manual(
    values = c(
      Exploratory = 17,
      Multi_lab   = 16
    ),
    labels = c(
      Exploratory = "Exploratory",
      Multi_lab   = "Multi-lab"
    )
  ) +
  labs(
    x = "Retrospective Dataset",
    y = "Effect Size (Hedge's g)",
    color = "Stage",
    shape = "Stage"
  ) +
  coord_cartesian(clip = "off") +
  # theme_minimal(base_family = "Arial", base_size = 22) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  theme_prism()

hedges_by_project_retrospective_plot

saveRDS(hedges_by_project_retrospective_plot, file.path(panels_dir, "hedges_effects_by_project_retrospective.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "hedges_effects_by_project_retrospective.png"),
  plot = hedges_by_project_retrospective_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# How many multi-lab projects are null?
retrospective_dt_wide[
  stage == "Multi_lab",
  .(
    total_multi_lab = .N,
    ci_spanning_0   = sum(ci_lower <= 0 & ci_upper >= 0),
    proportion      = mean(ci_lower <= 0 & ci_upper >= 0)
  )
]

##### ES per project and stage with pooled multi-lab — Retrospective ####

# Prepare data: one exploratory row + one pooled multi-lab row per project
retrospective_meta_plot <- rbindlist(list(
  retrospective_exploratory[, .(
    id,
    stage    = "Exploratory",
    g_pooled = hedges_g_exploratory,
    se       = se_g_exploratory,
    ci_lower = ci_lower_exploratory,
    ci_upper = ci_upper_exploratory
  )],
  meta_results_ext[, .(
    id,
    stage    = "Multi_lab",
    g_pooled = g_pooled_confirmatory,
    se       = se_g_confirmatory,
    ci_lower = ci_lower_confirmatory,
    ci_upper = ci_upper_confirmatory
  )]
), use.names = TRUE)

# Join plot_id for anonymous x axis
retrospective_meta_plot <- id_plotid_map[retrospective_meta_plot, on = "id"]

# Alphabetical factor order on plot_id
retrospective_meta_plot[, plot_id := factor(
  as.character(plot_id),
  levels = sort(unique(as.character(plot_id)))
)]
retrospective_meta_plot[, stage := factor(stage, levels = c("Exploratory", "Multi_lab"))]
retrospective_meta_plot[, base_x := as.numeric(plot_id)]

# x offset: exploratory left, pooled multi-lab right
retrospective_meta_plot[, jittered_x := fifelse(
  stage == "Exploratory",
  base_x - 0.15,
  base_x + 0.15
)]

hedges_pooled_by_project_ext <- ggplot(
  retrospective_meta_plot,
  aes(x = jittered_x, y = g_pooled, color = stage, shape = stage)
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.08
  ) +
  geom_vline(
    xintercept = seq(1.5, max(retrospective_meta_plot$base_x) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c(
      "Exploratory" = scales::alpha("#CC3300", 0.9),
      "Multi_lab"   = "#E7B800"
    ),
    labels = c(
      "Exploratory" = "Exploratory",
      "Multi_lab"   = "Multi-lab (pooled)"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Exploratory" = 17,
      "Multi_lab"   = 16
    ),
    labels = c(
      "Exploratory" = "Exploratory",
      "Multi_lab"   = "Multi-lab (pooled)"
    )
  ) +
  scale_x_continuous(
    breaks = sort(unique(retrospective_meta_plot$base_x)),
    labels = levels(retrospective_meta_plot$plot_id)
  ) +
  labs(
    x     = "Retrospective Dataset",
    y     = "Effect Size (Hedges' g)",
    color = "Stage",
    shape = "Stage"
  ) +
  theme_prism() +
  theme(
    axis.text.x        = element_text(hjust = 1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

hedges_pooled_by_project_ext

saveRDS(hedges_pooled_by_project_ext, file.path(panels_dir, "hedges_pooled_by_project_ext.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "hedges_pooled_by_project_retrospective.png"),
  plot     = hedges_pooled_by_project_ext,
  width    = 10,
  height   = 6,
  dpi      = 300
)

##### Combined Scatterplot Effect Sizes --------------------------------
# Goal: Plot ES shrinkage of pooled g (REML) of DECIDE vs retrospective

# Prepare dataset for plotting

retrospective_exploratory[, dataset := "retrospective"]
exploratory_decide[, dataset := "DECIDE"]



combined_retrospective_decide_g <- rbindlist(list(
                                          retrospective_exploratory[,
                                                               .(
                                                                 dataset,
                                                                 id,
                                                                 stage = "exploratory",
                                                                 effect_size = hedges_g_exploratory
                                                                 )
                                                               ],
                                          exploratory_decide[,
                                                             .(
                                                               dataset,
                                                               id = project_letter,
                                                               stage = "exploratory",
                                                               effect_size = hedges_g 
                                                                )
                                                             ],
                                          meta_results_decide[,
                                                              .(
                                                                dataset,
                                                                id = project_letter,
                                                                stage = "confirmatory/multi-lab",
                                                                effect_size = g_pooled
                                                                )
                                                              ],
                                          meta_results_ext[,
                                                           .(
                                                             dataset,
                                                             id,
                                                             stage = "confirmatory/multi-lab",
                                                             effect_size = g_pooled_confirmatory
                                                              )
                                                           ]
                                            ),
                                        use.names = TRUE
                                        )
combined_retrospective_decide_g_casted <- dcast(combined_retrospective_decide_g, 
                                           id + dataset ~ stage,
                                           value.var = "effect_size")
setDT(combined_retrospective_decide_g_casted)
setnames(combined_retrospective_decide_g_casted, "confirmatory/multi-lab", "multilab")
# S errors
combined_retrospective_decide_g_casted[sign(multilab) != sign(exploratory), 
                                       .(id, exploratory, multilab)]

# add a significant p-value column from the meta-analysis

setnames(meta_results_decide,"pval", "pval_confirmatory")
setnames(meta_results_decide,"project_letter", "id")

meta_p_values <- rbindlist(list(
  meta_results_ext[,.(dataset, id, pval_confirmatory)],
  meta_results_decide[,.(dataset, id, pval_confirmatory)]
))
combined_retrospective_decide_g_casted <- merge(combined_retrospective_decide_g_casted, 
      meta_p_values[, .(id, pval_confirmatory)],
      by = "id",
      all = TRUE)
combined_retrospective_decide_g_casted <- combined_retrospective_decide_g_casted[,
  p_value_significant := pval_confirmatory < 0.05
]


max_limit <- max(
  abs(combined_retrospective_decide_g_casted$exploratory),
  abs(combined_retrospective_decide_g_casted$multilab),
  na.rm = TRUE
)

max_limit <- max_limit * 1.05     # add 5% headroom


combined_scatterplot_es <- ggplot(combined_retrospective_decide_g_casted, 
                                  aes(x = abs(exploratory), 
                                      y = abs(multilab), 
                                      colour = p_value_significant,
                                      shape = dataset)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
  geom_point(size = 3, alpha = 0.8) +
  # geom_text_repel(aes(label = project_letter),
  #                 size = 4,
  #                 box.padding = 0.5,
  #                 point.padding = 0.3,
  #                 segment.color = "grey50"
  # ) +
  scale_color_manual(
    values = c("TRUE" = "#1F6F70", "FALSE" = "#B3B3B3"),
    labels = c("TRUE" = "significant", "FALSE" = "not significant"),
    name = "Meta-analytic p < 0.05"
  ) +
  scale_shape_manual(
    values = c("DECIDE" = 15, "retrospective" = 17),
    labels = c("DECIDE" = "Confirmatory", 
               "retrospective" = "Retrospective"),
    name = "Dataset"
  ) +
  # coord_fixed(ratio = 1, xlim = c(0, max_limit), ylim = c(0, max_limit), clip = "off") +
  scale_x_continuous(limits = c(0, max_limit)) +
  scale_y_continuous(limits = c(0, max_limit)) +
  labs(
    x = "|Exploratory Hedges'g|",
    y = "|Multi-lab Hedges'g|"
  ) +
  # scale_color_manual(values = c("retrospective" ="#E7B800", "DECIDE" = "#00AFBB"),
  #                   labels = c("retrospective" ="retrospective", "DECIDE" = "DECIDE"),
  #                   name = "framework"
  #                   ) +
  theme_prism() +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 22),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.position = "right",
    legend.box = "vertical",
    # legend.text = element_text(size = 22),
    # legend.title = element_text(size = 22),
    plot.margin = margin(20, 20, 20, 20),
    aspect.ratio = 1
  )

combined_scatterplot_es

## Move the legend outside the density plots: first remove the legend, then combine with cowplot:
# extract legend from the scatter plot
legend <- get_legend(
  combined_scatterplot_es + 
    theme(legend.position = "right",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14))
)

# make scatter without legend for ggMarginal
scatter_no_legend <- combined_scatterplot_es + theme(legend.position = "none")

marginal_plot <- ggMarginal(
  scatter_no_legend,
  type = "density",
  groupColour = TRUE,
  groupFill = TRUE,
  margins = "both",
  size = 5
)

# combine marginal plot with legend
final_scatter_plot <- plot_grid(
  marginal_plot,
  legend,
  nrow = 1,
  rel_widths = c(1, 0.3, 0.1)
)

final_scatter_plot

saveRDS(final_scatter_plot, file.path(panels_dir, "final_scatter_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "combined_scatterplot_es.png"),
  plot = final_scatter_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Spearman corr
cor_test <- cor.test(
  combined_retrospective_decide_g_casted$multilab,
  combined_retrospective_decide_g_casted$exploratory,
  data = combined_retrospective_decide_g_casted,
  method = "spearman", 
  exact = FALSE
)

cor_test

exploratory_means <- combined_retrospective_decide_g_casted[, .(
  mean_exploratory_g = mean(exploratory, na.rm = TRUE),
  sd_exploratory_g    = sd(exploratory, na.rm = TRUE),
  n_projects          = sum(!is.na(exploratory))
), by = dataset]

print(exploratory_means)

##### Meta-regression ####

decide_iv_wide <- dcast(decide_iv,
                        dataset + unique_lab ~ stage,
                        value.var = "iv_score")
setDT(decide_iv_wide)
# remove s100a1 project (no letter) as it doesn't have confirmatory data
decide_iv_wide <- decide_iv_wide[trimws(unique_lab) != ""]
setnames(decide_iv_wide, "Confirmatory", "multi_lab_iv")
setnames(decide_iv_wide, "Exploratory", "exploratory_iv")
setnames(decide_iv_wide, "unique_lab", "id")

retrospective_iv_wide <- dcast(retrospective_iv,
                               dataset + id ~ stage,
                               value.var = "iv_score")
setDT(retrospective_iv_wide)
setnames(retrospective_iv_wide, "Exploratory", "exploratory_iv")
setnames(retrospective_iv_wide, "Multi_lab", "multi_lab_iv")

combined_decide_retrospective_iv_wide <- rbind(decide_iv_wide, retrospective_iv_wide)

combined_decide_retrospective_iv_effectsize_wide <- merge(
  combined_decide_retrospective_iv_wide,
  combined_retrospective_decide_g_casted,
  by = c("id", "dataset"),
  all = TRUE
)

setnames(combined_decide_retrospective_iv_effectsize_wide, "multilab", "multi_lab_es")
setnames(combined_decide_retrospective_iv_effectsize_wide, "exploratory", "exploratory_es")

# add CI columns from meta dt
meta_combined <- rbindlist(list(
  meta_results_decide[,
                      .(
                        dataset,
                        id = id,
                        ci_lower_confirmatory = ci_lower,
                        ci_upper_confirmatory = ci_upper
                      )],
  meta_results_ext[,
                   .(
                     dataset,
                     id,
                     ci_lower_confirmatory,
                     ci_upper_confirmatory
                   )]
))

combined_decide_retrospective_iv_effectsize_ci_wide <- merge(
  combined_decide_retrospective_iv_effectsize_wide,
  meta_combined,
  by = c("id", "dataset"),
  all = TRUE
)


# rename shorter
es_dt <- combined_decide_retrospective_iv_effectsize_ci_wide

# Create a Standard Error column from the CI
es_dt$multi_lab_se <- (es_dt$ci_upper_confirmatory - es_dt$ci_lower_confirmatory) / 3.92


# Create the meta-analysis object
m_gen <- metagen(TE = multi_lab_es,         # The Multi-lab effect size
                 seTE = multi_lab_se,       # The Standard Error of Multi-lab effect
                 data = es_dt,
                 studlab = id,          # Optional: Name of the projects
                 sm = "SMD",            # Standardized Mean Difference
                 common = FALSE,        # In this case I am comparing across projects, they are not estimating the same thing like the meta-analysis of multi-lab so I am doing a random effects model here not a fixed effects
                 random = TRUE          # Fixed effects are the datasets, random effects the single projects
                 )                    

# Run Meta-regression with the 2 datasets as a Fixed Effect
# y is the multi-lab effect size of a project i
# fixed effects: exploratory es and dataset
# random effects: project
# Dataset here is as a factor that explains why effects differ
# Note: in meta-regression the random effects are inherited from the meta-analysis object (m_gen)
# they won't appear as (1|effect) like in generalized linear models

# Fixed model, without intercept (adding -1) to better model that if ES is 0 in exploratory it's also 0 in confirmatory)
# This allows to only see the slope difference between the 2 datasets
# The effect size will be only a fx of the exploratory predictor (removed dataset predictor)
m_reg <- metareg(m_gen, ~ exploratory_es + exploratory_es: dataset - 1)

summary(m_reg)


# Extract model results
reg_table <- as.data.table(coef(summary(m_reg)), keep.rownames = "term")

# Compute effective DECIDE slope (sum of both coefficients)
decide_slope <- reg_table[1, estimate] + reg_table[2, estimate]
decide_se    <- sqrt(reg_table[1, se]^2 + reg_table[2, se]^2)  # approximate SE
decide_z     <- decide_slope / decide_se
decide_p     <- 2 * pnorm(abs(decide_z), lower.tail = FALSE)
decide_ci_lb <- decide_slope - 1.96 * decide_se
decide_ci_ub <- decide_slope + 1.96 * decide_se

# Bind computed row
derived_row <- data.table(
  term     = "exploratory_es:datasetDECIDE (effective DECIDE slope)",
  estimate = decide_slope,
  se       = decide_se,
  zval     = decide_z,
  pval     = decide_p,
  ci.lb    = decide_ci_lb,
  ci.ub    = decide_ci_ub
)

reg_table <- rbind(reg_table, derived_row)

# Clean up term names
reg_table[, term := c(
  "Exploratory ES — Retrospective slope",
  "Exploratory ES × Confirmatory (slope difference)",
  "Exploratory ES — Confirmatory effective slope"
)]

# Store formatted table
reg_table_formatted <- reg_table[, .(
  Term     = term,
  β        = estimate,
  SE       = se,
  z        = zval,
  p        = pval,
  `95% CI` = paste0("[", round(ci.lb, 2), ", ", round(ci.ub, 2), "]")
)]

# Render and save
reg_table_formatted |>
  kbl(
    digits  = 3,
    align   = c("l", "r", "r", "r", "r", "r"),
    caption = paste0(
      "Meta-regression results. Residual heterogeneity: ",
      "τ² = ", round(m_reg$tau2, 3), ", ",
      "I² = ", round(m_reg$I2, 1), "%, ",
      "QE(df = ", m_reg$k - m_reg$p, ") = ", round(m_reg$QE, 2),
      ", p < 0.001"
    )
  ) |>
  kable_classic(full_width = FALSE, html_font = "Cambria") |>
  row_spec(3, italic = TRUE, color = "grey40") |>
  footnote(
    general = paste0(
      "Test of moderators: QM(df = 2) = ", round(m_reg$QM, 2),
      ", p < 0.001. Model fitted without intercept. ",
      "Effective Confirmatory slope computed as sum of coefficients with approximate SE."
    ),
    general_title = ""
  ) |>
  save_kable(file = file.path(save_dir_external, "meta_regression_table.html"))

write.table(
  reg_table_formatted,
  file      = file.path(save_dir_external, "meta_regression_table.tsv"),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)

#### P-value plots________________________####

combined_retrospective_decide_pval <- rbindlist(list(
  retrospective_all[,
               .(
                 id,
                 hedges_g_exploratory,
                 se_g_exploratory,
                 p_val_exploratory,
                 g_pooled_confirmatory,
                 se_g_confirmatory,
                 pval_confirmatory,
                 dataset
               )],
  meta_wide[,
            .(
              id = project_letter,
              hedges_g_exploratory = exploratory_g,
              se_g_exploratory = exploratory_se,
              p_val_exploratory = exploratory_pval,
              g_pooled_confirmatory = confirmatory_g_pooled,
              se_g_confirmatory = confirmatory_se,
              pval_confirmatory = confirmatory_pval,
              dataset = "DECIDE"
            )]
  ), use.names = TRUE
  )

# reshape to long format and rename for the boxplot
combined_pval_long <- melt(
  combined_retrospective_decide_pval,
  id.vars = c("id", "dataset"),
  measure.vars = c("p_val_exploratory", "pval_confirmatory"),
  variable.name = "stage",
  value.name = "p_value")


setDT(combined_pval_long)
combined_pval_long[, stage := fifelse(
  stage == "p_val_exploratory",
  "exploratory", "confirmatory"
)]

# create group variable
combined_pval_long[, group := interaction(dataset, stage, sep = "_")]

combined_pval_long[, group := factor(
  group,
  levels = c(
    "retrospective_exploratory",
    "retrospective_confirmatory",
    "DECIDE_exploratory",
    "DECIDE_confirmatory"
    )
  )
  ]


stat_table_pval <- make_pw_wilcox_stat_tbl(
  data = combined_pval_long,
  y = "p_value",
  group = "group"
)

stat_table_pval_selected <- stat_table_pval[
  group1 == "retrospective_exploratory" & group2 == "retrospective_confirmatory" |
    group1 == "DECIDE_exploratory" & group2 == "DECIDE_confirmatory" |
    group1 == "retrospective_confirmatory" & group2 == "DECIDE_confirmatory"
  
]
cols_dataset <- c("DECIDE" = "#00AFBB", "retrospective" = "#E7B800")
stat_table_pval_selected <-
  stat_table_pval_selected <-
  stat_table_pval_selected |>
  arrange(y.position, group1, group2) |>
  mutate(y.position = c(0.7, 0.9, 1.7))

combined_pval_long_plot <-
  ggplot(
    combined_pval_long,
    aes(x = group, y = p_value, fill = dataset)
  ) +
  geom_boxplot(
    aes(alpha = stage),
    outlier.shape = NA,
    width = 0.55,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) +
  scale_fill_manual(values = cols_dataset) +
  scale_alpha_manual(
    values = c("exploratory" = 0.35, "confirmatory" = 0.85)
  ) +
  scale_x_discrete(
    labels = c(
      "retrospective_exploratory"  = "Exploratory",
      "retrospective_confirmatory" = "Multi-lab",
      "DECIDE_exploratory"         = "Exploratory",
      "DECIDE_confirmatory"        = "Multi-lab"
    )
  ) +
  labs(
    x = NULL,
    y = "p-values (two sided t-test)",
    fill = NULL
  ) +
  # retro exp vs retro conf
  annotate("segment", x = 1, xend = 2, y = 8, yend = 8) +
  annotate("segment", x = 1, xend = 1, y = 8, yend = 5) +
  annotate("segment", x = 2, xend = 2, y = 8, yend = 5) +
  annotate("text", x = 1.5, y = 8, label = stat_table_pval_selected$p_label[1],
           vjust = -0.5, size = 4) +
  # decide exp vs decide conf
  annotate("segment", x = 3, xend = 4, y = 8, yend = 8) +
  annotate("segment", x = 3, xend = 3, y = 8, yend = 5) +
  annotate("segment", x = 4, xend = 4, y = 8, yend = 5) +
  annotate("text", x = 3.5, y = 8, label = stat_table_pval_selected$p_label[2],
           vjust = -0.5, size = 4) +
  # retro conf vs decide conf
  annotate("segment", x = 2, xend = 4, y = 50, yend = 50) +
  annotate("segment", x = 2, xend = 2, y = 50, yend = 15) +
  annotate("segment", x = 4, xend = 4, y = 50, yend = 15) +
  annotate("text", x = 3, y = 50, label = stat_table_pval_selected$p_label[3],
           vjust = -0.5, size = 4) +
  scale_y_log10(
    breaks = c(0.001, 0.01, 0.05, 0.1, 0.5, 1),
    labels = c("0.001", "0.01", "0.05", "0.1", "0.5", "1"),
    limits = c(NA, 120)
  ) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "grey40"
  ) +
  annotate(
    "text", x = 1.5, y = 0, label = "Retrospective",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  annotate(
    "text", x = 3.5, y = 0, label = "Confirmatory",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  # ggbreak::scale_y_break(
  #   breaks = c(0.00005, 0.0005),
  #   scales = 0.2
  # ) +
  # coord_cartesian(ylim = c(0.00001, 120)) +
  theme_prism() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.y.right  = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.title.y.right = element_blank(),
    plot.margin = margin(30, 40, 50, 10)
  )

combined_pval_long_plot

ggsave(
  filename = file.path(save_dir_external, "p_values_decide_retrospective.png"),
  plot = combined_pval_long_plot,
  width = 10,
  height = 6,
  dpi = 300
)


kruskal.test(p_value ~ group, data = combined_pval_long)

pairwise.wilcox.test(
  combined_pval_long$p_value,
  combined_pval_long$group,
  p.adjust.method = "BH"   # Benjamini–Hochberg correction
)

combined_pval_long[, significant := fifelse(p_value < 0.05, "yes", "no")]

combined_pval_long[
  ,
  rate_yes := mean(significant == "yes"),
  by = group
]
rate_dt <- combined_pval_long[
  ,
  .(rate_yes = mean(p_value < 0.05, na.rm = TRUE)),
  by = .(group, dataset, stage)
]

## Fisher's exact test

tab_retrospective <- combined_pval_long[
  dataset == "retrospective",
  table(stage, significant)
]
tab_retrospective
p_retrospective <- fisher.test(tab_retrospective)$p.value

tab_decide <- combined_pval_long[
  dataset == "DECIDE",
  table(stage, significant)
]
tab_decide
p_decide <- fisher.test(tab_decide)$p.value

tab_retrospective_decide <- combined_pval_long[
  group %in% c("retrospective_confirmatory",
               "DECIDE_confirmatory"),
  table(group, significant)
]

tab_retrospective_decide
p_retrospective_decide <- fisher.test(tab_retrospective_decide)$p.value
p_raw <- c(p_retrospective, p_decide, p_retrospective_decide)
p_holm <- p.adjust(p_raw, method = "holm")

p_holm

stat_fisher_holm <- data.table(
  group1 = c("retrospective_exploratory", "DECIDE_exploratory", "retrospective_confirmatory"),
  group2 = c("retrospective_confirmatory", "DECIDE_confirmatory", "DECIDE_confirmatory"),
  p = p_raw,
  p.adj = p_holm
)
stat_fisher_holm[, p_label := paste0("p = ", signif(p.adj, 3))]
stat_fisher_holm[, y.position := c(1.1, 1, 1.2)]

##### p-value rate plot ####
significant_rate_plot <-
  ggplot(
    rate_dt,
    aes(
      x = group,
      y = rate_yes,
      fill = dataset,
      alpha = stage
    )
  ) +
  geom_col(
    position = position_dodge2(width = 0.7, preserve = "single"),
    width = 0.6,
    colour = "black",
    linewidth = 0.6
  ) +
  # geom_text(
  #   aes(label = scales::percent(rate_yes, accuracy = 1)),
  #   vjust = -0.35,
  #   size = 5
  # ) +
  ggpubr::stat_pvalue_manual(
    stat_fisher_holm,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4,
    vjust = -0.5
  )+
  scale_fill_manual(values = cols_dataset) +
  scale_alpha_manual(
    values = c("exploratory" = 0.35, "confirmatory" = 0.85)
  ) +
  scale_x_discrete(
    labels = c(
      "retrospective_exploratory"  = "Exploratory",
      "retrospective_confirmatory" = "Multi-lab",
      "DECIDE_exploratory"         = "Exploratory",
      "DECIDE_confirmatory"        = "Multi-lab"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1.2),
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = "Proportion significant (p < 0.05)",
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    "text", x = 1.5, y = 0, label = "Retrospective",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  annotate(
    "text", x = 3.5, y = 0, label = "Confirmatory",
    vjust = 5, size = 5, fontface = "bold"
  ) +
  theme_prism() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.text.x = element_text(vjust = 0.5),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 16),
    plot.margin = margin(10, 40, 50, 10)
  )

significant_rate_plot

ggsave(
  filename = file.path(save_dir_external, "significant_rate_decide_ext.png"),
  plot = significant_rate_plot,
  width = 10,
  height = 6,
  dpi = 300
)

##### Z-curve analysis ####
# Estimate average true power for exploratory and confirmatory
# and expected replication rate from the exploratory

# Print out the ODR (% significant)
combined_pval_long[, 
                   .(ODR = mean(p_value < 0.05, na.rm = TRUE),
                     n   = .N),
                   by = .(dataset, stage)]

# Convert p-values to z-values for z-curve
zcurve_both <- combined_pval_long[!is.na(p_value), .(
  z     = qnorm(1 - p_value / 2),  # two-sided p to z
  stage = stage,
  dataset = dataset
)]

# Cap z at 6 (as z-curve does)
zcurve_both[z > 6, z := 6]

## Fit z-curves: pooling the 2 datasets as it need at least z-scores in the fitting range)

# ── Fit z-curve for exploratory (ERR meaningful here) ────────────────────────
# I can't do z-curve for confirmatory as it need at least 10 significant p-values and we have 9
zfit_exploratory <- zcurve(p = combined_pval_long[
  stage == "exploratory" & !is.na(p_value), p_value
])



# Summary notes:
# model: EM is the fitting algorithm: Expectation - Maximizazion (default)
# ERR: Expected Replication Rate: average post-hoc power across
# ODR: Observed Discovery Rate: it's % of significant
# Q: goodness of fit statistic. Negative Q means the z-scores are closer to the
#    significance threshold (z = 1.96).
summary(zfit_exploratory)
plot(zfit_exploratory)

##  Extract the underlying data and rebuild the plot in ggplot2 for the manuscript
z_seq <- seq(0, 6, by = 0.05)
fitted_expl <- get_fitted_density(zfit_exploratory, z_seq)
fitted_expl[, stage := "exploratory"]



# Scale factor: align density to histogram counts
n_expl <- zcurve_both[stage == "exploratory", .N]
n_conf <- zcurve_both[stage == "confirmatory", .N]
binwidth <- 0.5

fitted_expl[stage == "exploratory", density_scaled := density * n_expl * binwidth]

# ERR label for exploratory only
ERR      <- round(zfit_exploratory$coefficients["ERR"], 2)
ERR_ci   <- round(quantile(zfit_exploratory$coefficients_boot$ERR, c(0.025, 0.975)), 2)
err_label <- paste0("ERR = ", ERR, " [", ERR_ci[1], ", ", ERR_ci[2], "]")

stage_colors <- c("exploratory" = "#CC3300", "confirmatory" = "#00AFBB")
stage_labels <- c("exploratory" = "Exploratory", "confirmatory" = "Multi-lab")

# exploratory fitted curve + both histograms
zcurve_plot <- ggplot() +
  # Histograms for both stages
  geom_histogram(
    data = zcurve_both,
    aes(x = z, fill = stage),
    binwidth = binwidth,
    boundary = 0,
    color = "white",
    alpha = 0.6,
    position = "identity"
  ) +
  # Fitted curve for exploratory only
  geom_line(
    data = fitted_expl[, density_scaled := density * n_expl * binwidth],
    aes(x = z, y = density_scaled),
    color = "#CC3300",
    linewidth = 1.2
  ) +
  geom_vline(
    xintercept = 1.96,
    color = "gray40",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  annotate(
    "text", x = 3.8, y = Inf,
    label = err_label,
    hjust = 0, vjust = 1.5,
    size = 4, color = "#CC3300"
  ) +
  scale_fill_manual(values = stage_colors, labels = stage_labels) +
  scale_x_continuous(breaks = 0:6, limits = c(0, 6.2)) +
  labs(x = "Z-score", y = "Count", fill = NULL, title = NULL) +
  theme_prism() +
  theme(
    legend.position    = "bottom",
    panel.grid.major.x = element_blank()
  )

zcurve_plot

saveRDS(zcurve_plot, file.path(panels_dir, "zcurve_plot.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "zcurve_exploratory_vs_confirmatory.png"),
  plot     = zcurve_plot,
  width    = 8,
  height   = 5,
  dpi      = 300
)

# Summary table
zcurve_summary <- data.table(
  Stage         = "Exploratory",
  N_supplied    = zfit_exploratory$N_obs,
  N_significant = zfit_exploratory$N_sig,
  ODR           = round(zfit_exploratory$N_sig / zfit_exploratory$N_obs, 2),
  ERR           = round(zfit_exploratory$coefficients["ERR"], 2),
  ERR_low       = round(quantile(zfit_exploratory$coefficients_boot$ERR, 0.025), 2),
  ERR_high      = round(quantile(zfit_exploratory$coefficients_boot$ERR, 0.975), 2),
  EDR           = round(zfit_exploratory$coefficients["EDR"], 2),
  EDR_low       = round(quantile(zfit_exploratory$coefficients_boot$EDR, 0.025), 2),
  EDR_high      = round(quantile(zfit_exploratory$coefficients_boot$EDR, 0.975), 2),
  Note          = "ERR/EDR estimated via z-curve EM. Confirmatory stage had insufficient significant results (n=9) for z-curve fitting."
)

print(zcurve_summary)

write.xlsx(
  zcurve_summary,
  file.path(save_dir_external, "zcurve_summary.xlsx"),
  rowNames = FALSE
)

# Replication outcome counts (not in the manuscript)____________________________________----------------------------------------------

# Count total successes per project


# Remove SESOI criterion from decide as it's not possible for the retrospective (no SESOI data available)
replication_matrix_decide[, `CI above SESOI`:= NULL]

# Get logical columns
log_cols <- names(replication_matrix_decide)[sapply(replication_matrix_decide, is.logical)]

# count successes
replication_matrix_decide[, success_counts := rowSums(.SD, na.rm = TRUE), .SDcols = log_cols]
 
replication_matrix_retrospective[, success_counts := rowSums(.SD, na.rm = TRUE), .SDcols = log_cols]

combined_success_counts <- rbindlist(list(
  replication_matrix_retrospective[,
                            .(
                              dataset = "retrospective",
                              success_counts,
                              id = id
                            )],
  replication_matrix_decide[,
                            .(
                              dataset = "DECIDE",
                              success_counts,
                              id = project_letter
                            )]
))


combined_success_counts_plot <- ggplot(combined_success_counts, aes(x = dataset, y = success_counts, fill = dataset)) +
  geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) +
  scale_fill_manual(values = c("darkseagreen", "brown4")) +
  labs(
    title = "Replication success counts",
    x = "Dataset",
    y = "Success counts"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

combined_success_counts_plot

ggsave(
  filename = file.path(save_dir_external, "combined_success_counts_plot.png"),
  plot = combined_success_counts_plot,
  width = 10,
  height = 6,
  dpi = 300
)




# Correlation of best simulation replication criteria (Expl in conf) (not in the manuscript)____________________________________ --------

# Remove study 12 as the IV is duplicated with 11
combined_decide_retrospective_iv_effectsize_ci_wide <- combined_decide_retrospective_iv_effectsize_ci_wide[id != 12]
comb_iv_es <- combined_decide_retrospective_iv_effectsize_ci_wide

# add same sign column
comb_iv_es[, direction_agreement := sign(exploratory_es) == sign(multi_lab_es)
           & exploratory_es != 0 & multi_lab_es != 0]

# add expl es in confirmatory CI criteria
comb_iv_es[, expl_in_ci := exploratory_es >= ci_lower_confirmatory &
             exploratory_es <= ci_upper_confirmatory & direction_agreement]

# Create labeled factor for x-axis
comb_iv_es[, expl_in_ci_label := factor(expl_in_ci, levels = c(FALSE, TRUE), labels = c("not met", "met"))]

# Relabel dataset for facetting
comb_iv_es$dataset <- factor(
  comb_iv_es$dataset,
  levels = c("retrospective", "DECIDE"),
  labels = c("Retrospective", "Confirmatory")
)

# Build combined group variable using relabeled dataset
comb_iv_es$group <- interaction(comb_iv_es$dataset, comb_iv_es$expl_in_ci, sep = " - ")
comb_iv_es$group <- factor(
  comb_iv_es$group,
  levels = c(
    "Retrospective - TRUE",
    "Retrospective - FALSE",
    "Confirmatory - TRUE",
    "Confirmatory - FALSE"
  )
)

stat_tbl_expl_in_ci_vs_iv <- make_pw_wilcox_stat_tbl(
  data = comb_iv_es,
  y = "multi_lab_iv",
  group = "group"
)

stat_tbl_expl_in_ci_vs_iv_selected <- stat_tbl_expl_in_ci_vs_iv[
  group1 == "Retrospective - TRUE" & group2 == "Retrospective - FALSE" |
    group1 == "Confirmatory - TRUE" & group2 == "Confirmatory - FALSE"
]

# Split back the groups for facetting
stat_tbl_expl_in_ci_vs_iv_selected[, dataset := fifelse(
  grepl("Confirmatory", group1), "Confirmatory", "Retrospective"
)]
stat_tbl_expl_in_ci_vs_iv_selected[, group1 := gsub(".*- ", "", group1)]
stat_tbl_expl_in_ci_vs_iv_selected[, group2 := gsub(".*- ", "", group2)]

# Relabel group1 and group2 to match plot x-axis labels
stat_tbl_expl_in_ci_vs_iv_selected[, group1 := factor(group1, levels = c("FALSE", "TRUE"), labels = c("not met", "met"))]
stat_tbl_expl_in_ci_vs_iv_selected[, group2 := factor(group2, levels = c("FALSE", "TRUE"), labels = c("not met", "met"))]

# Match dataset factor levels
stat_tbl_expl_in_ci_vs_iv_selected$dataset <- factor(
  stat_tbl_expl_in_ci_vs_iv_selected$dataset,
  levels = c("Retrospective", "Confirmatory")
)
cols_dataset <- c("Confirmatory" = "#00AFBB", "Retrospective" = "#E7B800")

### Plot
iv_expl_ci_plot <-
  ggplot(
    comb_iv_es,
    aes(x = expl_in_ci_label, y = multi_lab_iv, fill = dataset)
  ) +
  geom_boxplot(
    width = 0.5,
    linewidth = 0.5,
    colour = "black",
    outlier.shape = NA
  ) +
  ggbeeswarm::geom_beeswarm(
    size = 2,
    cex = 1.1,
    colour = "black",
    priority = "density"
  ) +
  facet_wrap(~dataset, nrow = 1) +
  ggpubr::stat_pvalue_manual(
    stat_tbl_expl_in_ci_vs_iv_selected,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    facet.by = "dataset",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4
  ) +
  scale_fill_manual(values = cols_dataset) +
  coord_cartesian(ylim = c(2.5, 9.5), clip = "off") +
  labs(
    x = "Exploratory ES in confirmatory CI",
    y = "Minimal internal validity (mIV)"
  ) +
  guides(fill = "none") +
  theme_prism() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 15, face = "bold"),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 15, 10, 10)
  )

iv_expl_ci_plot

ggsave(
  filename = file.path(save_dir_external, "Expl_in_CI_vs_internal_validity_plot.png"),
  plot = iv_expl_ci_plot,
  width = 10,
  height = 6,
  dpi = 300
)
# Pooling datasets
iv_expl_ci_plot <-
  ggplot(
    comb_iv_es,
    aes(x = expl_in_ci, y = multi_lab_iv, fill = dataset)
  ) +
  geom_boxplot(
    position = position_dodge(width = 0.7),
    outlier.shape = NA,
    width = 0.55,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_beeswarm(
    size = 1.8,
    colour = "black",
    priority = "density",
    cex = 1.2
  ) +
  facet_wrap(~ factor(dataset, levels = c("retrospective", "DECIDE")), nrow = 1) +
  scale_x_discrete(labels = c(`TRUE` = "met", `FALSE` = "not met")) +
  scale_fill_manual(values = cols_dataset) +
  guides(fill = "none") +
  ggpubr::stat_pvalue_manual(
    stat_tbl_expl_in_ci_vs_iv_selected,
    label = "p_label",
    xmin = "group1",
    xmax = "group2",
    y.position = "y.position",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 4,
    step.increase = 0
  ) +
  labs(
    x = "Exploratory ES in confirmatory CI",
    y = "minimal internal validity (mIV)"
  ) +
  theme_prism() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14, face = "bold")
  )
iv_expl_ci_plot

#### Define a continuous variable for the criteria ES within confirmatory CI:
# This will be the column c which is the distance between exploratory ES and confirmatory ES
# If the exploratory is outside the confirmatory CI and/or if there is a sign error, a penalty is added scaled to CI halfwidth
# The higher the c the worse



# add distance between effect sizes
comb_iv_es[, es_distance := abs(exploratory_es - multi_lab_es)]

# CI half-width column to scale the penalties
comb_iv_es[, ci_halfwidth := (ci_upper_confirmatory - ci_lower_confirmatory) / 2]


# Penalties (to be scaled to the CI halfwidth)
s_error = 1 # for opposite sign
out_of_ci = 0.5 # for being outside the

# c score: distance between expl es and conf es: the higher the worse

comb_iv_es[, c := fcase(
  # CASE 1: if expl is within the CI and same sign, score is the abs distance from the confirmatory ES
  direction_agreement & expl_in_ci,
  es_distance,
  
  # CASE 2: if same sign but out of the CI, add penalty
  direction_agreement & !expl_in_ci,
  es_distance + out_of_ci * ci_halfwidth,
  
  # CASE 3: edge case: if s error but still in CI (confirmatory not significant) add penalty. the score is slightly worse when the distance between exploratory ES and confirmatory ES spans 0 as compare with a similar distance but same sign.
  !direction_agreement & expl_in_ci,
  es_distance + s_error * ci_halfwidth,
  
  # CASE 4: if s error and out of CI add both penalties
  !direction_agreement & !expl_in_ci,
  es_distance + (out_of_ci + s_error) * ci_halfwidth,
  
  default = NA_real_
)]


# Distribution is not monotonic so I can't do a Spearman test
ggplot(
  comb_iv_es,
  aes(x = exploratory_iv, y = c)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  labs(
    x = "Exploratory internal validity",
    y = "Replication discordance score (c)"
  ) +
  theme_minimal()

ggplot(
  comb_iv_es,
  aes(x = multi_lab_iv, y = c)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  labs(
    x = "Multi-lab internal validity",
    y = "Replication discordance score (c)"
  ) +
  theme_minimal()

# plot es distance with iv score
ggplot(
  comb_iv_es,
  aes(x = exploratory_iv, y = es_distance)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  labs(
    x = "Exploratory internal validity",
    y = "ES distance"
  ) +
  theme_minimal()

ggplot(
  comb_iv_es,
  aes(x = multi_lab_iv, y = es_distance)
) +
  geom_point(alpha = 0.7) +
  geom_smooth(
    method = "loess",
    se = FALSE,
    color = "black"
  ) +
  labs(
    x = "Multi-lab internal validity",
    y = "ES distance"
  ) +
  theme_minimal()

# correlation es distance in DECIDE
# cor.test(
#   comb_iv_es[dataset == "DECIDE", multi_lab_iv],
#   comb_iv_es[dataset == "DECIDE", es_distance],
#   method = "kendall"
# )

# # correlation in DECIDE
# cor.test(
#   comb_iv_es[dataset == "DECIDE", multi_lab_iv],
#   comb_iv_es[dataset == "DECIDE", c],
#   method = "spearman"
# )
# 
# # correlation in retrospective
# cor.test(
#   comb_iv_es[dataset == "retrospective", multi_lab_iv],
#   comb_iv_es[dataset == "retrospective", c],
#   method = "spearman"
# )

# Kendall correlation: ordinal and tie-robust (more suitable here than Spearman)

# exploratory_iv, both datasets
cor.test(
  comb_iv_es$exploratory_iv,
  comb_iv_es$c,
  method = "kendall",
  exact = FALSE
)
# multi-lab iv, both datasets
cor.test(
  comb_iv_es$multi_lab_iv,
  comb_iv_es$c,
  method = "kendall",
  exact = FALSE
)
# multi-lab iv, retrospective
cor.test(
  comb_iv_es[dataset == "Retrospective", multi_lab_iv],
  comb_iv_es[dataset == "Retrospective", c],
  method = "kendall",
  exact = FALSE
)
# multi-lab iv, DECIDE
cor.test(
  comb_iv_es[dataset == "Confirmatory", multi_lab_iv],
  comb_iv_es[dataset == "Confirmatory", c],
  method = "kendall",
  exact = FALSE
)
# exploratory iv, retrospective
cor.test(
  comb_iv_es[dataset == "Retrospective", exploratory_iv],
  comb_iv_es[dataset == "Retrospective", c],
  method = "kendall",
  exact = FALSE
)
# exploratory iv, DECIDE
cor.test(
  comb_iv_es[dataset == "Confirmatory", exploratory_iv],
  comb_iv_es[dataset == "Confirmatory", c],
  method = "kendall",
  exact = FALSE
)

# categoriting c score
comb_iv_es[, c_cat := fcase(
  c <= quantile(c, 0.50, na.rm = TRUE), "low",
  c <= quantile(c, 0.75, na.rm = TRUE), "moderate",
  default = "high"
)]

comb_iv_es[, c_cat := factor(c_cat, levels = c("low", "moderate", "high"), ordered = TRUE)]

comb_iv_es[, expl_iv_ord := factor(exploratory_iv, ordered = TRUE)]
comb_iv_es[, conf_iv_ord := factor(multi_lab_iv, ordered = TRUE)]

# contingency tables
table(comb_iv_es$expl_iv_ord, comb_iv_es$c_cat)
table(comb_iv_es$conf_iv_ord, comb_iv_es$c_cat)

## delta ES vs delta IV
comb_iv_es[, delta_iv := multi_lab_iv - exploratory_iv]
comb_iv_es[, delta_es := multi_lab_es - exploratory_es]

table(sign(comb_iv_es$delta_es))
table(sign(comb_iv_es$delta_iv))

## Conditional sign test
# Since IV sign is completely skewed (15/20 delta increase), Kendall here is uninformative
# Restrict to IV improved studies to see in how many studies es decrease
comb_iv_es_imp <- comb_iv_es[delta_iv > 0]
table(sign(comb_iv_es_imp$delta_es))

# one-sided binomial sign test
binom.test(
  sum(comb_iv_es_imp$delta_es < 0),
  nrow(comb_iv_es_imp),
  p = 0.5,
  alternative = "greater"
)



# Variance Analysis 2______________________####
# Var. analysis for the retrospective dataset

##### CV controls exploratory vs pooled confirmatory #####

# Exploratory: single-lab cv for controls
ctrl_exploratory_ext <- retrospective_castable[group_es == "ctrl" & stage == "Exploratory",
                                               .(mean_expl      = mean,
                                                 sd_expl        = sd,
                                                 n_expl         = sample_size,
                                                 cv_exploratory = sd / abs(mean)),
                                               by = .(id)]
ctrl_exploratory_ext <- id_plotid_map[ctrl_exploratory_ext, on = "id"]

# Confirmatory: pooled across labs for controls
ctrl_confirmatory_ext <- retrospective_castable[group_es == "ctrl" & stage == "Multi_lab",
                                                {
                                                  pooled_mean <- mean(mean)
                                                  pooled_sd   <- sqrt(sum((sample_size - 1) * sd^2 + sample_size * (mean - pooled_mean)^2) / (sum(sample_size) - 1))
                                                  .(mean_pooled = pooled_mean,
                                                    sd_pooled   = pooled_sd,
                                                    n_total     = sum(sample_size))
                                                },
                                                by = .(id, plot_id)]

# Remove NA project #2 as I can't compute pooled means from % disease free data
ctrl_confirmatory_ext <- ctrl_confirmatory_ext[!is.na(mean_pooled)]
ctrl_confirmatory_ext[, cv_confirmatory := sd_pooled / abs(mean_pooled)]

ctrl_confirmatory_ext[, plot_id := as.factor(plot_id)]
merged_ctrl_ext <- ctrl_exploratory_ext[ctrl_confirmatory_ext, on = .(id, plot_id)]
cv_merged_ctrl_ext <- merged_ctrl_ext[, .(id, plot_id, cv_exploratory, cv_confirmatory)]

var_compare_ext_long <- data.table::melt(cv_merged_ctrl_ext,
                                         measure.vars  = c("cv_exploratory", "cv_confirmatory"),
                                         variable.name = "Stage",
                                         value.name    = "cv")
var_compare_ext_long[, Stage := fifelse(Stage == "cv_exploratory", "Exploratory", "Multi_lab")]

cv_ext_pooled_labs_ctrl_plot <- ggplot(
  var_compare_ext_long,
  aes(x = factor(plot_id), y = cv, color = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, length(unique(var_compare_ext_long$plot_id)) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("Exploratory" = scales::alpha("#CC3300", 0.9), "Multi_lab" = "#E7B800"),
    labels = c("Exploratory" = "Exploratory", "Multi_lab" = "Multi-lab")
  ) +
  labs(title = NULL, x = "Project", y = "Variance (CV)", color = "Stage") +
  theme_prism() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

cv_ext_pooled_labs_ctrl_plot


### CV difference expl vs pooled multi-lab ###

cv_diff_pooled_ext_dt <- cv_merged_ctrl_ext[
  , .(cv_diff = cv_confirmatory - cv_exploratory),
  by = .(plot_id)
]

cv_diff_ext_pooled_test <- t.test(cv_diff_pooled_ext_dt$cv_diff, mu = 0)
wilcox.test(cv_diff_pooled_ext_dt$cv_diff, mu = 0)

p_val <- cv_diff_ext_pooled_test$p.value
sig_label_ext <- fcase(
  p_val < 0.001, "p < 0.001",
  p_val < 0.01,  "p < 0.01",
  p_val < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_val, 3))
)

cv_diff_ctrl_pooled_ext_plot <- ggplot(cv_diff_pooled_ext_dt, aes(x = "", y = cv_diff)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#CC3300", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#CC3300", width = 0.15, size = 2) +
  # annotate("text", x = Inf, y = Inf, label = sig_label_ext,
  #          hjust = 1.1, vjust = 1.5, size = 4) +
  labs(x = NULL, y = "CV difference\n(pooled multi-lab - exploratory)", title = NULL) +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

cv_ext_pooled_labs_ctrl_plot | cv_diff_ctrl_pooled_ext_plot


### Per-project CV test ###

cv_test_per_project_ext <- merged_ctrl_ext[,
                                           {
                                             se_cv_expl <- cv_exploratory / sqrt(2 * n_expl)  * sqrt(1 + 2 * cv_exploratory^2)
                                             se_cv_conf <- cv_confirmatory / sqrt(2 * n_total) * sqrt(1 + 2 * cv_confirmatory^2)
                                             se_diff    <- sqrt(se_cv_expl^2 + se_cv_conf^2)
                                             cv_diff    <- cv_confirmatory - cv_exploratory
                                             t_stat     <- cv_diff / se_diff
                                             df         <- (se_cv_expl^2 + se_cv_conf^2)^2 /
                                               (se_cv_expl^4 / (n_expl - 1) + se_cv_conf^4 / (n_total - 1))
                                             p_value    <- pt(-abs(t_stat), df = df)
                                             .(cv_diff, t_stat, df, p_value)
                                           },
                                           by = .(id, plot_id)
]

cv_test_per_project_ext[, sig_label := fcase(
  p_value < 0.001, "p < 0.001",
  p_value < 0.01,  "p < 0.01",
  p_value < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_value, 3))
)]


##### Mean difference exploratory vs pooled confirmatory (not in the manuscript) #####

merged_ctrl_ext_long <- data.table::melt(merged_ctrl_ext,
                                         measure.vars  = c("mean_expl", "mean_pooled"),
                                         variable.name = "Stage",
                                         value.name    = "mean")
merged_ctrl_ext_long[, Stage := fifelse(Stage == "mean_expl", "Exploratory", "Multi_lab")]

mean_ext_pooled_labs_ctrl_plot <- ggplot(
  merged_ctrl_ext_long,
  aes(x = factor(plot_id), y = mean, color = Stage)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  geom_vline(
    xintercept = seq(1.5, length(unique(merged_ctrl_ext_long$plot_id)) - 0.5, by = 1),
    linetype = "solid", color = "gray80", linewidth = 0.5
  ) +
  scale_color_manual(
    values = c("Exploratory" = scales::alpha("#CC3300", 0.9), "Multi_lab" = "#E7B800"),
    labels = c("Exploratory" = "Exploratory", "Multi_lab" = "Multi-lab")
  ) +
  labs(title = "Control", x = "Project", y = "Mean", color = "Stage") +
  theme_prism() +
  theme(axis.text.x = element_text(hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

mean_ext_pooled_labs_ctrl_plot


### Hedges' g per project ###

merged_ctrl_ext <- merged_ctrl_ext[,
                                   {
                                     mean_diff        <- mean_pooled - mean_expl
                                     sd_ratio         <- sd_pooled / sd_expl
                                     sd_pooled_hedges <- sqrt(((n_expl - 1) * sd_expl^2 + (n_total - 1) * sd_pooled^2) /
                                                                (n_expl + n_total - 2))
                                     d       <- mean_diff / sd_pooled_hedges
                                     df      <- n_expl + n_total - 2
                                     J       <- 1 - (3 / (4 * df - 1))
                                     g       <- d * J
                                     se_g    <- sqrt((n_expl + n_total) / (n_expl * n_total) + g^2 / (2 * (n_expl + n_total)))
                                     t_stat  <- g / se_g
                                     p_value <- 2 * pt(-abs(t_stat), df = df)
                                     .(mean_diff, sd_pooled, sd_expl, sd_ratio, g, se_g, t_stat, p_value)
                                   },
                                   by = .(id, plot_id)
]

merged_ctrl_ext[, sig_label := fcase(
  p_value < 0.001, "p < 0.001",
  p_value < 0.01,  "p < 0.01",
  p_value < 0.05,  "p < 0.05",
  default = paste0("p = ", round(p_value, 3))
)]

overall_test_ext <- t.test(merged_ctrl_ext$g, mu = 0)

hedges_g_ext_plot_ctrl <- ggplot(merged_ctrl_ext, aes(x = "", y = g)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
  geom_boxplot(fill = "#CC3300", alpha = 0.6, outlier.shape = NA) +
  geom_jitter(color = "#CC3300", width = 0.15, size = 2) +
  labs(x = NULL, y = "Hedges' g\n(multi-lab - exploratory)", title = "Control") +
  theme_prism() +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

mean_ext_pooled_labs_ctrl_plot | hedges_g_ext_plot_ctrl


### Forest plot CV difference vs Hedges' g ###

cv_forest_ext <- cv_test_per_project_ext[, .(plot_id,
                                             estimate = cv_diff,
                                             se       = cv_diff / t_stat,
                                             metric   = "CV difference")]

g_forest_ext <- merged_ctrl_ext[, .(plot_id,
                                    estimate = g,
                                    se       = se_g,
                                    metric   = "Hedges' g")]

forest_ext_dt <- rbind(cv_forest_ext, g_forest_ext)
forest_ext_dt[, ci_lower := estimate - 1.96 * se]
forest_ext_dt[, ci_upper := estimate + 1.96 * se]
forest_ext_dt[, plot_id := factor(plot_id, levels = sort(unique(plot_id)))]

forest_ext_plot <- ggplot(forest_ext_dt,
                          aes(x = estimate, y = plot_id, color = metric)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper),
                  position = position_dodge(width = 0.5), size = 0.6) +
  scale_color_manual(
    values = c("CV difference" = scales::alpha("#CC3300", 0.9), "Hedges' g" = "#E7B800"),
    name = NULL
  ) +
  labs(x = "Effect size (with 95% CI)", y = "Project",
       title = "CV difference vs Hedges' g per project") +
  theme_prism() +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(color = "gray90"),
        panel.grid.major.x = element_blank())

forest_ext_plot


### Decomposition plot ###

decomp_ext_long <- data.table::melt(
  merged_ctrl_ext[, .(plot_id, mean_diff, sd_ratio, g)],
  id.vars       = "plot_id",
  variable.name = "metric",
  value.name    = "value"
)

decomp_ext_long[, metric := fcase(
  metric == "mean_diff", "Mean difference (numerator)",
  metric == "sd_ratio",  "SD ratio (denominator)",
  metric == "g",         "Hedges' g (combined)"
)]

decomp_ext_long[, metric := factor(metric, levels = c(
  "Mean difference (numerator)",
  "SD ratio (denominator)",
  "Hedges' g (combined)"
))]

decomp_ext_long[, plot_id := factor(plot_id, levels = rev(sort(unique(plot_id))))]
decomp_ext_long[, ref := fifelse(metric == "SD ratio (denominator)", 1, 0)]

decomp_ext_long[, metric := factor(metric, levels = c(
  "Hedges' g (combined)",
  "Mean difference (numerator)",
  "SD ratio (denominator)"
))]

decomp_ext_plot <- ggplot(decomp_ext_long, aes(x = value, y = plot_id, color = metric)) +
  geom_vline(aes(xintercept = ref), linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "Mean difference (numerator)" = scales::alpha("#003366", 0.9),
    "SD ratio (denominator)"      = "#CC3300",
    "Hedges' g (combined)"        = "#00AFBB"
  )) +
  facet_wrap(~ metric, scales = "free_x", ncol = 3) +
  labs(x = NULL, y = "Project", title = "Decomposition of effect size — Controls (Retrospective)") +
  theme_prism() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        panel.grid.major.y = element_line(color = "gray90"))

decomp_ext_plot

# Effect Size deconstruction__________________####

##### Shrinkage indicators - normalized mean difference (not included in the manuscript) ####
# Indicators of Hedges g, mean diff, sd, ctrl mean diff, ctrl sd
# with formula |exploratory| - |confirmatory| / |exploratory| + |confirmatory|
# < 0: confirmatory inflated
# 0: same magnitude
# 0 - 1: confirmatory smaller than exploratory
# 1: total shrinkage
# >1: sign error

### Hedges' g, mean diff, SD of ctrl vs Treated ###

expl_ctrl_treated_ext <- retrospective_castable[
  stage == "Exploratory" & group_es == "ctrl",
  .(id, g_expl = hedges_g, mean_diff_expl = mean_difference, sd_expl = pooled_sd)
]

conf_g_ext <- meta_results_ext[id_plotid_map, on = "id"][,
                                                         .(id, plot_id, g_conf = g_pooled_confirmatory)]

conf_mean_sd_ext <- retrospective_castable[
  stage == "Multi_lab" & group_es == "ctrl",
  {
    pooled_mean_diff <- mean(mean_difference, na.rm = TRUE)
    pooled_sd_val    <- sqrt(sum((sample_size - 1) * pooled_sd^2 +
                                   sample_size * (mean_difference - pooled_mean_diff)^2,
                                 na.rm = TRUE) /
                               (sum(sample_size, na.rm = TRUE) - 1))
    .(mean_diff_conf = pooled_mean_diff, sd_conf = pooled_sd_val)
  },
  by = id
][id_plotid_map, on = "id"]

shrinkage_ctrl_treated_ext <- expl_ctrl_treated_ext[conf_g_ext, on = "id"]
shrinkage_ctrl_treated_ext <- shrinkage_ctrl_treated_ext[conf_mean_sd_ext, on = .(id, plot_id)]

shrinkage_ctrl_treated_ext[, `:=`(
  g_shrinkage    = (abs(g_expl)         - abs(g_conf))         / (abs(g_expl)         + abs(g_conf)),
  mean_shrinkage = (abs(mean_diff_expl) - abs(mean_diff_conf)) / (abs(mean_diff_expl) + abs(mean_diff_conf)),
  sd_shrinkage   = (abs(sd_expl)        - abs(sd_conf))        / (abs(sd_expl)        + abs(sd_conf))
)]

merged_ctrl_ext_shrinkage <- merged_ctrl_ext[
  ctrl_exploratory_ext[, .(id, plot_id, mean_expl, sd_expl)], on = .(id, plot_id)
][
  ctrl_confirmatory_ext[, .(id, plot_id, mean_pooled, sd_pooled)], on = .(id, plot_id)
]
merged_ctrl_ext_shrinkage[, plot_id := as.factor(plot_id)]

merged_ctrl_ext_shrinkage[, `:=`(
  ctrl_mean_shrinkage = (abs(mean_expl) - abs(mean_pooled)) / (abs(mean_expl) + abs(mean_pooled))
)]

# Remove project 2 (pooled means not computable from % disease free data)
shrinkage_ctrl_treated_ext <- shrinkage_ctrl_treated_ext[plot_id != 2]
merged_ctrl_ext_shrinkage  <- merged_ctrl_ext_shrinkage[plot_id != 2]

multi_lab_dt <- id_plotid_map[multi_lab_dt, on = "id"]
multi_lab_dt               <- multi_lab_dt[plot_id != 2]

# Project factor order (shared across all panels)
project_levels_ext <- rev(sort(unique(as.character(shrinkage_ctrl_treated_ext$plot_id))))



# ── Panel 1: Hedges g (Ctrl vs Treated) ──────────────────────────────────────
p1_ext_dt <- shrinkage_ctrl_treated_ext[,
                                        .(plot_id = as.character(plot_id), value = g_shrinkage,
                                          metric = "Effect Site Shrinkage")]
p1_ext_dt[, plot_id := factor(plot_id, levels = project_levels_ext)]

p1_ext <- ggplot(p1_ext_dt, aes(x = value, y = plot_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = "#00AFBB") +
  facet_wrap(~ metric) +
  labs(x = NULL, y = "Retrospective Dataset") +
  shrinkage_theme()

# ── Panel 2: Mean diff (Ctrl vs Treated) ─────────────────────────────────────
p2_ext_dt <- shrinkage_ctrl_treated_ext[,
                                        .(plot_id = as.character(plot_id), value = mean_shrinkage,
                                          metric = "Mean Shrinkage")]
p2_ext_dt[, plot_id := factor(plot_id, levels = project_levels_ext)]

p2_ext <- ggplot(p2_ext_dt, aes(x = value, y = plot_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = scales::alpha("#003366", 0.9)) +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  shrinkage_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# ── Panel 3: SD (Ctrl vs Treated) ────────────────────────────────────────────
p3_ext_dt <- shrinkage_ctrl_treated_ext[,
                                        .(plot_id = as.character(plot_id), value = sd_shrinkage,
                                          metric = "Variance Inflation")]
p3_ext_dt[, plot_id := factor(plot_id, levels = project_levels_ext)]

p3_ext <- ggplot(p3_ext_dt, aes(x = value, y = plot_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = "#CC3300") +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  shrinkage_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# ── Panel 4: Mean diff (Ctrl Expl vs Conf) ───────────────────────────────────
p4_ext_dt <- merged_ctrl_ext_shrinkage[,
                                       .(plot_id = as.character(plot_id), value = ctrl_mean_shrinkage,
                                         metric = "Control Stability")]
p4_ext_dt[, plot_id := factor(plot_id, levels = project_levels_ext)]

p4_ext <- ggplot(p4_ext_dt, aes(x = value, y = plot_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
  geom_point(size = 3, color = "#E7B800") +
  facet_wrap(~ metric) +
  labs(x = NULL, y = NULL) +
  shrinkage_theme() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.y = element_blank())

# ── Panel 5: Distance method (Niangoran 2023) — mean difference as variable Y ──
# Variable Y = mean difference (ctrl - treated) per center
# SS_i = (n_i - 1)*sd_i^2 + n_i*(mean_diff_i - grand_mean_diff)^2
# where n_i = n_ctrl + n_treated per center, sd_i = pooled_sd
# Di = (SS_i / df1_i) / (SS_total / df2) ~ F(df1_i, df2)

alpha_dist <- 0.05

# One row per center: mean_difference and pooled_sd are on the Ctrl row
distance_es_ext_dt <- retrospective_castable[
  stage == "Multi_lab" & group_es == "ctrl"   # ctrl row holds mean_diff and pooled_sd
][
  id_plotid_map, on = "id"                    # join anonymous plot_id
][
  !is.na(plot_id) & plot_id != 2             # drop project 2
][,
  {
    # Grand mean difference: weighted average across centers by total n
    y_grand  <- sum(EU * mean_difference) / sum(EU)
    
    # Total sample size across all centers
    N_total  <- sum(EU)
    
    # Numerator df per center
    df1      <- EU - 1
    
    # Denominator df: total observations minus 1
    df2      <- N_total - 1
    
    # Per-center SS: within-center (from pooled_sd) + between-center component
    SS_i     <- (EU - 1) * pooled_sd^2 + EU * (mean_difference - y_grand)^2
    
    # Total SS
    SS_total <- sum(SS_i)
    
    # Di statistic
    Di       <- (SS_i / df1) / (SS_total / df2)
    
    # Critical F value per center
    f_crit   <- qf(1 - alpha_dist, df1 = df1, df2 = df2)
    
    # Flag atypical centers
    atypical <- Di > f_crit
    
    .(center, EU, mean_difference, pooled_sd, y_grand, Di, f_crit, atypical)
  },
  by = .(id, plot_id)
]

# Align factor levels with other panels
distance_es_ext_dt[, plot_id := factor(
  as.character(plot_id),
  levels = project_levels_ext
)]

# Strip label to match other panels
distance_es_ext_dt[, metric := "Atypical Labs"]
distance_decide_dt[, atypical := factor(atypical, levels = c("FALSE", "TRUE"))]

p5_ext <- ggplot(
  distance_es_ext_dt,
  aes(x = Di, y = plot_id, color = atypical)
) +
  # Di = 1: center exactly at the grand mean difference
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray60") +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("FALSE" = "#7B2D8B", "TRUE" = "pink2"),
    labels = c("FALSE" = "Typical", "TRUE" = "Atypical (p < 0.05)"),
    drop   = FALSE
  ) +
  facet_wrap(~ metric) +
  labs(x = "Distance statistic (Dᵢ)", y = "Retrospective Dataset" , title = NULL) +
  shrinkage_theme() +
  theme(
    # axis.text.y     = element_blank(),
    # axis.ticks.y    = element_blank(),
    # axis.line.y     = element_blank(),
    strip.text = element_blank(),
    legend.position = "bottom",
    legend.title    = element_blank()
  )

p5_ext

saveRDS(p5_ext, file.path(panels_dir, "p5_ext.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "atypical_labs_ext_plot.png"),
  plot     = p5_ext,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# ── Combine all 5 ────────────────────────────────────────────────────────────
shrinkage_decomp_ext_plot <- p1_ext + p2_ext + p3_ext + p4_ext + p5_ext +
  plot_layout(nrow = 1, widths = c(1.4, 1, 1, 1, 1), guides = "collect") &
  theme(legend.position = "none")


shrinkage_decomp_ext_plot

ggsave(
  filename = file.path(save_dir_external, "shrinkage_decomp_ext_plot.png"),
  plot     = shrinkage_decomp_ext_plot,
  width    = 14,
  height   = 6,
  dpi      = 300
)

# Summary table of all shrinkage indicators

# Hedges g, mean diff, SD shrinkage (from shrinkage_ctrl_treated_ext)
shrinkage_summary <- shrinkage_ctrl_treated_ext[,
                                                .(plot_id,
                                                  g_expl         = round(g_expl, 3),
                                                  g_conf         = round(g_conf, 3),
                                                  g_shrinkage    = round(g_shrinkage, 3),
                                                  mean_diff_expl = round(mean_diff_expl, 3),
                                                  mean_diff_conf = round(mean_diff_conf, 3),
                                                  mean_shrinkage = round(mean_shrinkage, 3),
                                                  sd_expl        = round(sd_expl, 3),
                                                  sd_conf        = round(sd_conf, 3),
                                                  sd_shrinkage   = round(sd_shrinkage, 3))
]

# Ctrl mean shrinkage (from merged_ctrl_ext_shrinkage)
ctrl_mean_summary <- merged_ctrl_ext_shrinkage[,
                                               .(plot_id,
                                                 ctrl_mean_expl       = round(mean_expl, 3),
                                                 ctrl_mean_conf       = round(mean_pooled, 3),
                                                 ctrl_mean_shrinkage  = round(ctrl_mean_shrinkage, 3))
]
ctrl_mean_summary[, plot_id := as.character(plot_id)]
shrinkage_summary[,  plot_id := as.character(plot_id)]

# Distance Di: one row per project — use median Di across centers
distance_summary <- distance_es_ext_dt[,
                                       .(n_centers      = .N,
                                         n_atypical     = sum(atypical),
                                         median_Di      = round(median(Di), 3),
                                         any_atypical   = any(atypical)),
                                       by = .(plot_id)
]
distance_summary[, plot_id := as.character(plot_id)]

# Merge all
summary_table_ext <- shrinkage_summary[
  ctrl_mean_summary, on = "plot_id"
][
  distance_summary, on = "plot_id"
]

setorder(summary_table_ext, plot_id)

# Display with gt
summary_table_ext_dt <- summary_table_ext[, .(
  Project                  = plot_id,
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
)] |>
  gt() |>
  tab_header(title = "Shrinkage indicators — Retrospective dataset") |>
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
  # Hedges g: severe shrinkage (> 0.5) or negative
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Hedges g shrinkage",
      rows    = `Hedges g shrinkage` > 0.5 | `Hedges g shrinkage` < 0
    )
  ) |>
  # Mean diff: severe shrinkage (> 0.5) or negative
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Mean diff shrinkage",
      rows    = `Mean diff shrinkage` > 0.5 | `Mean diff shrinkage` < 0
    )
  ) |>
  # SD: negative inflation beyond -0.1 threshold
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "SD shrinkage",
      rows    = `SD shrinkage` < -0.1
    )
  ) |>
  # Ctrl mean: drift beyond ±0.2
  tab_style(
    style     = cell_fill(color = "#f8d7da"),
    locations = cells_body(
      columns = "Ctrl mean shrinkage",
      rows    = `Ctrl mean shrinkage` > 0.2 | `Ctrl mean shrinkage` < -0.2
    )
  ) |>
  # Distance: atypical centers flagged (only for projects with > 2 centers)
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
  summary_table_ext_dt,
  file.path(save_dir_external, "shrinkage_summary_retrospective.xlsx"),
  rowNames = FALSE
)


##### Factor contribution: Shapley values ####
# Retrospective dataset
retrospective_castable <- id_plotid_map[retrospective_castable, on = "id"]

# Exploratory components (single lab)
expl_raw_ext <- retrospective_castable[
  stage == "Exploratory" & group_es %in% c("ctrl", "treated"),
  .(mean = mean, sd = sd, n = sample_size, group_es),
  by = .(id, plot_id)
]

expl_ctrl_ext <- expl_raw_ext[group_es == "ctrl",
                              .(id, plot_id,
                                mu_ctrl_exp = mean, sd_ctrl_exp = sd, n_ctrl_exp = n)]

expl_treated_ext <- expl_raw_ext[group_es == "treated",
                                 .(id, plot_id,
                                   mu_treated_exp = mean, sd_treated_exp = sd, n_treated_exp = n)]

expl_comp_ext <- expl_ctrl_ext[expl_treated_ext, on = .(id, plot_id)]

expl_comp_ext[, sd_exp_pooled := sqrt(
  ((n_ctrl_exp - 1) * sd_ctrl_exp^2 + (n_treated_exp - 1) * sd_treated_exp^2) /
    (n_ctrl_exp + n_treated_exp - 2)
)]

expl_comp_ext[, g_exp := hedges_g_from_parts(
  mu_ctrl_exp, mu_treated_exp, sd_exp_pooled, n_ctrl_exp, n_treated_exp
)]

# Confirmatory components (pooled across labs)
# conf_ctrl_ext <- retrospective_castable[
#   stage == "Multi_lab" & group_es == "ctrl", {
#     pm  <- mean(mean)
#     psd <- sqrt(sum((sample_size - 1) * sd^2 + sample_size * (mean - pm)^2) /
#                   (sum(sample_size) - 1))
#     .(mu_ctrl_conf = pm, sd_ctrl_conf = psd, n_ctrl_conf = sum(sample_size))
#   }, by = .(id, plot_id)]

# Confirmatory: weighted mean across labs (weighted by sample size)
conf_ctrl_ext <- retrospective_castable[
  stage == "Multi_lab" & group_es == "ctrl", {
    w   <- sample_size / sum(sample_size)          # pesi = n_i / N_total
    pm  <- sum(w * mean)                           # weighted mean
    psd <- sqrt(sum((sample_size - 1) * sd^2 + 
                      sample_size * (mean - pm)^2) /
                  (sum(sample_size) - 1))
    .(mu_ctrl_conf = pm, sd_ctrl_conf = psd, n_ctrl_conf = sum(sample_size))
  }, by = .(id, plot_id)]

# conf_treated_ext <- retrospective_castable[
#   stage == "Multi_lab" & group_es == "treated", {
#     pm  <- mean(mean)
#     psd <- sqrt(sum((sample_size - 1) * sd^2 + sample_size * (mean - pm)^2) /
#                   (sum(sample_size) - 1))
#     .(mu_treated_conf = pm, sd_treated_conf = psd, n_treated_conf = sum(sample_size))
#   }, by = .(id, plot_id)]

# Confirmatory treated: weighted mean across labs (weighted by sample size)
conf_treated_ext <- retrospective_castable[
  stage == "Multi_lab" & group_es == "treated", {
    w   <- sample_size / sum(sample_size)
    pm  <- sum(w * mean)
    psd <- sqrt(sum((sample_size - 1) * sd^2 +
                      sample_size * (mean - pm)^2) /
                  (sum(sample_size) - 1))
    .(mu_treated_conf = pm, sd_treated_conf = psd, n_treated_conf = sum(sample_size))
  }, by = .(id, plot_id)]

conf_comp_ext <- conf_ctrl_ext[conf_treated_ext, on = .(id, plot_id)]

conf_comp_ext[, sd_conf_pooled := sqrt(
  ((n_ctrl_conf - 1) * sd_ctrl_conf^2 + (n_treated_conf - 1) * sd_treated_conf^2) /
    (n_ctrl_conf + n_treated_conf - 2)
)]

conf_comp_ext[, g_conf := hedges_g_from_parts(
  mu_ctrl_conf, mu_treated_conf, sd_conf_pooled, n_ctrl_conf, n_treated_conf
)]

# Merge
decomp_ext <- expl_comp_ext[conf_comp_ext, on = .(id, plot_id)]

# Remove project 2 (% disease free, means not computable)
decomp_ext <- decomp_ext[plot_id != 2]

decomp_ext[, total_shrinkage := abs(g_exp) - abs(g_conf)]

# All 7 coalition |g| values
decomp_ext[, `:=`(
  g0    = abs(g_exp),
  
  g_c   = abs(hedges_g_from_parts(mu_ctrl_conf,  mu_treated_exp,  sd_exp_pooled,
                                  n_ctrl_exp,  n_treated_exp)),
  g_t   = abs(hedges_g_from_parts(mu_ctrl_exp,   mu_treated_conf, sd_exp_pooled,
                                  n_ctrl_exp,  n_treated_exp)),
  g_s   = abs(hedges_g_from_parts(mu_ctrl_exp,   mu_treated_exp,  sd_conf_pooled,
                                  n_ctrl_conf, n_treated_conf)),
  
  g_ct  = abs(hedges_g_from_parts(mu_ctrl_conf,  mu_treated_conf, sd_exp_pooled,
                                  n_ctrl_exp,  n_treated_exp)),
  g_cs  = abs(hedges_g_from_parts(mu_ctrl_conf,  mu_treated_exp,  sd_conf_pooled,
                                  n_ctrl_conf, n_treated_conf)),
  g_ts  = abs(hedges_g_from_parts(mu_ctrl_exp,   mu_treated_conf, sd_conf_pooled,
                                  n_ctrl_conf, n_treated_conf)),
  
  g_cts = abs(g_conf)
)]

# Shapley values
decomp_ext[, `:=`(
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

decomp_ext[, `:=`(
  phi_ctrl    = -phi_ctrl,
  phi_treated = -phi_treated,
  phi_sd      = -phi_sd
)]

# check that the sum of phis is equal to the shrinkage difference
decomp_ext[, check := phi_ctrl + phi_treated + phi_sd]
print(decomp_ext[, .(plot_id,
                     total_shrinkage = round(total_shrinkage, 3),
                     check           = round(check, 3),
                     diff            = round(total_shrinkage - check, 8))])

# Plot
project_ord_ext <- rev(sort(unique(as.character(decomp_ext$plot_id))))

ext_long <- melt(
  decomp_ext[, .(plot_id = as.character(plot_id),
                 `Control stability`      = phi_ctrl,
                 `Treatment response` = phi_treated,
                 `Variance inflation` = phi_sd)],
  id.vars = "plot_id", variable.name = "component", value.name = "phi"
)
setDT(ext_long)

set(ext_long, j = "plot_id", value = factor(ext_long$plot_id, levels = project_ord_ext))
set(ext_long, j = "component", value = factor(ext_long$component,
                                              levels = c("Control stability", "Treatment response", "Variance inflation")))

p_shapley_ext <- ggplot(ext_long,
                        aes(x = phi, y = plot_id, fill = component)) +
  geom_col(position = "stack", width = 0.7) +
  geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
  scale_fill_manual(values = component_colors) +
  scale_x_continuous(
    sec.axis = dup_axis(
      breaks = c(-2, 4),
      labels = c("← reducing", "contributing to shrinkage →"),
      name   = NULL
    )
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
    axis.text.x.top      = element_text(size = 14, margin = margin(b = 8))
  )

p_shapley_ext

saveRDS(p_shapley_ext, file.path(panels_dir, "p_shapley_ext.rds"))   # ADD

ggsave(
  filename = file.path(save_dir_external, "shrinkage_shapley_decomp_ext.png"),
  plot     = p_shapley_ext,
  width    = 10,
  height   = 7,
  dpi      = 300
)

# Total shrinkage
p_total_ext <- ggplot(decomp_ext,
                      aes(x = total_shrinkage,
                          y = factor(as.character(plot_id), levels = project_ord_ext))) +
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

p_total_ext + p_shapley_ext + plot_layout(ncol = 2, widths = c(1, 2))

ggsave(
  filename = file.path(save_dir_external, "shrinkage_shapley_decomp_ext.png"),
  plot     = p_total_ext + p_shapley_ext + plot_layout(ncol = 2, widths = c(1, 2)),
  width    = 14,
  height   = 7,
  dpi      = 300
)
