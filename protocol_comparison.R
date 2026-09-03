source(here("packages.R"))
source(here("all_functions.R"))
data_dir        <- here("data") 
save_dir_decide <- here("results")

protocol_comparison <- fread(file.path(data_dir, "protocol_comparison_anonymized.csv"))
validity_levels <- c("IV", "EV", "SV", "TV")
# Create a unique row_id
# protocol_comparison[, row_id := paste(detail, aspect, validity, sep = " | ")]

# Condense information of detail and aspect in a new column

protocol_comparison[, row_id := fcase(
  # Blinding
  detail == "Data analysis" & aspect == "Blinding",      "Blinded Analysis",
  detail == "Outcome assessment" & aspect == "Blinding", "Blinded Assessment",
  detail == "Intervention giving" & aspect == "Blinding","Blinded Intervention",
  detail == "Group allocation" & aspect == "Blinding",   "Blinded Allocation",
  
  # Controls
  detail == "Sham mock naive",      "Sham/Mock/Naive Control",
  detail == "Comparator",           "Comparator Control",
  detail == "Positive control",     "Positive Control",
  detail == "Negative control",     "Negative Control",
  
  # Population
  detail == "Exclusion criteria",   "Exclusion Criteria",
  detail == "Inclusion criteria",   "Inclusion Criteria",
  detail == "Health status",        "Health Status",
  detail == "Sex",                  "Sex",
  detail == "Age",                  "Age",
  
  # Randomization
  detail == "Randomization method",  "Randomization Method",
  detail == "Data analysis scoring", "Randomized Analysis",
  detail == "Intervention giving" & aspect == "Randomization", "Randomized Intervention",
  detail == "Group allocation" & aspect == "Randomization",    "Randomized Allocation",
  
  # Others
  detail == "Multi lab approach",   "Multi-lab Approach",
  detail == "Statistical advice",   "Statistical Advice",
  detail == "Outcome a priori",    "A priori Outcome",
  detail == "Power calculation",    "Power Calculation",
  detail == "Outcome readout",      "Outcome Readout",
  detail == "Disease simulation",   "Disease Simulation",
  detail == "Scheme",               "Intervention Scheme",
  detail == "Strain",               "Model Strain",
  detail == "Species",              "Model Species",
  
  # Fallback
  default = as.character(detail)
)]


# Set row_id levels in the order of validity_levels then aspect
row_levels <- protocol_comparison[
  order(factor(validity, levels = validity_levels), aspect),
  unique(row_id)
]
# protocol_comparison[, row_id := gsub("_", " ", row_id)]
protocol_comparison[, row_id := factor(row_id, levels = row_levels)]
protocol_comparison[, row_id := fct_rev(row_id)] # to avoid that ggplot invert the order

protocol_comparison[, rel_score := score / max_score]

# Use 'detail' as the actual label
protocol_comparison[, validity := factor(validity, levels = validity_levels)]

# Ensure the details are ordered correctly within each validity group
protocol_comparison <- protocol_comparison[order(validity, aspect, row_id)]
protocol_comparison[, row_id := factor(row_id, levels = unique(row_id))]

protocol_comparison$phase <- factor(
  protocol_comparison$phase,
  levels = c("Exploratory Stage", "Confirmatory Stage")
)

# supervised heatmap ------------------------------------------------------

##### raw plot ####
heatmap_raw <- ggplot(
  protocol_comparison, 
  aes(
    x = project_letter, 
    y = row_id,
    fill = rel_score
  )
) +
  geom_tile(color = "white", linewidth = 0.15) +
  # Facet by Phase (columns) AND Validity (rows)
  facet_grid(
    rows = vars(validity), 
    cols = vars(phase), 
    scales = "free", 
    space = "free" # CRITICAL: This keeps tile sizes uniform
  ) +
  scale_fill_gradient2(
    low      = "#8ECFD4",
    mid      = "#2E7D87",
    high     = "#1B3A3E",
    midpoint = 0.5,
    limits   = c(0, 1),
    breaks   = c(0, 1)
  ) +
  labs(
    x = "pCS",
    y = NULL,
    fill = "Score"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 10, face = "bold", vjust = 1.5),
    # Style the side labels (IV, EV, etc.)
    strip.text.y = element_text(size = 15, face = "bold", angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold"),
    panel.grid = element_blank(),
    panel.spacing = unit(0.2, "lines"), # Reduce gap between facets
    legend.position = "right"
  )

heatmap_raw

saveRDS(heatmap_raw, file.path(save_dir_decide, "panels", "heatmap_raw.rds"))

ggsave(
  filename = file.path(save_dir_decide, "protocol_comparison_heatmap_raw.png"),
  plot = heatmap_raw,
  width = 15,   
  height = 9,  
  dpi = 300,
  bg = "white"
)

# Simplified supervised heatmap (only validity classification) -----------------------

protocol_comparison[, validity_aspect := paste(aspect, validity, sep = " | ")]

protocol_comparison_simplified <- protocol_comparison[, 
                    .(
                      rel_score_mean = mean(rel_score, na.rm = TRUE)),
                    by = .(phase, validity_aspect, project_letter, validity, aspect)]


# Set row levels
row_levels <- protocol_comparison_simplified[
  order(
    factor(validity, levels = validity_levels),
    aspect
  ),
  unique(validity_aspect)
]
protocol_comparison_simplified[, validity_aspect := factor(validity_aspect, levels = rev(row_levels))] # rev because ggplot inverts the order in the heatmap


heatmap_raw_simplified <- ggplot(
  protocol_comparison_simplified,
  aes(
    x = project_letter,
    y = validity_aspect,
    fill = rel_score_mean
  )
) +
  geom_tile(color = "white", linewidth = 0.15, width = 0.95, height = 1) +
  facet_grid(cols = vars(phase)) +
  scale_fill_gradient(low = "#D6F4F6", 
                       high = "#008A93",
                       breaks = c(0, 1)) +
  coord_fixed() +
  labs(x = "Confirmatory Studies", y = NULL, fill = "Score") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 15)),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 24, face = "bold"),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22, margin = margin(b = 8))
  )

heatmap_raw_simplified

ggsave(
  filename = file.path(save_dir_decide, "protocol_comparison_heatmap_raw_simplified.png"),
  plot = heatmap_raw_simplified,
  width = 15,   
  height = 9,  
  dpi = 300,
  bg = "white"
)

# unsupervised heatmap ------------------------------------------------------

# wide matrix: rows = row_id, columns = phase + project
protocol_comparison[, row_id_inv := paste(validity, aspect, detail, sep = " | ")]

# enforce ordering: validity -> aspect -> detail
row_levels_inv <- protocol_comparison[
  order(factor(validity, levels = validity_levels), aspect, row_id),
  unique(row_id)
]
protocol_comparison[, row_id_ordered := factor(row_id, levels = row_levels_inv)]

# cast
protocol_comparison_wide <- dcast(
  protocol_comparison,
  row_id_ordered+ validity ~ phase + project_letter,
  value.var = "rel_score"
)

# exploratory
setnames(
  protocol_comparison_wide,
  old = grep("^Exploratory Stage_", names(protocol_comparison_wide), value = TRUE),
  new = sub("^Exploratory Stage_", "", grep("^Exploratory Stage_", names(protocol_comparison_wide), value = TRUE)) |>
    paste0("[Expl] ")
)

# confirmatory
setnames(
  protocol_comparison_wide,
  old = grep("^Confirmatory Stage_", names(protocol_comparison_wide), value = TRUE),
  new = sub("^Confirmatory Stage_", "", grep("^Confirmatory Stage_", names(protocol_comparison_wide), value = TRUE)) |>
    paste0("[Conf] ")
)

old <- names(protocol_comparison_wide)

new <- old
new <- sub("^([A-Z])\\[Expl\\]\\s*$", "[Expl] \\1", new)
new <- sub("^([A-Z])\\[Conf\\]\\s*$", "[Conf] \\1", new)

setnames(protocol_comparison_wide, old = old, new = new)
setDT(protocol_comparison_wide)
# matrix
# mat <- as.matrix(protocol_comparison_wide[, -1]) # remove id col for matrix conversion
# rownames(mat) <- protocol_comparison_wide$row_id_inv # reattach it for labeling

cols_to_keep <- setdiff(names(protocol_comparison_wide), c("row_id_ordered", "validity"))
mat <- as.matrix(protocol_comparison_wide[, ..cols_to_keep])

rownames(mat) <- protocol_comparison_wide$row_id_ordered

# right-side validity annotation
validity_vec <- protocol_comparison_wide$validity
validity_colors <- c(IV = "#4E79A7", EV = "#F28E2B", SV = "#E15759", TV = "#76B7B2")
row_anno <- rowAnnotation(
  Validity = validity_vec,
  col = list(Validity = validity_colors),
  annotation_name_gp = grid::gpar(fontsize = 12, fontface = "bold"),
  annotation_legend_param = list(
    Validity = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                    labels_gp = gpar(fontsize = 11))
  )
)
# color mapping
# col_fun <- colorRamp2(
#   c(0, 0.5, 1),
#   c("#d8f3dc","#b3cd7a","#3a7d44")
# )
col_fun <- colorRamp2(
  c(0, 0.5, 1),
  c("#8ECFD4","#2E7D87","#1B3A3E")
)

# png(
#   filename = file.path(save_dir_decide, "protocol_comparison_heatmap_unsupervised.png"),
#   width = 13, height = 9, units = "in", res = 300, bg = "white"  
# )

unsupervised_heatmap <- Heatmap(
  mat,
  name = "Score",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = grid::gpar(fontsize = 10),
  column_title_side = "bottom",
  column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
  row_names_side = "left",          # labels now on the left
  right_annotation = row_anno       # validity on the right
)
unsupervised_heatmap

p_unsupervised_heatmap <- grid.grabExpr({
  draw(unsupervised_heatmap,
       padding = unit(c(20, 5, 5, 20), "mm"),
       heatmap_legend_side = "left")
  grid.text(
    "Confirmatory Studies",
    x = 0.5,
    y = unit(12, "mm"),
    gp = gpar(fontsize = 16, fontface = "bold")
  )
})

p_unsupervised_heatmap

saveRDS(p_unsupervised_heatmap, file.path(save_dir_decide, "panels", "p_unsupervised_heatmap.rds"))

# dev.off()






# Plot change from exploratory to confirmatory: delta confirmati --------


# Achieve a wide format table with multiple values per Stage
protocol_comparison_wide_delta <- dcast(
  protocol_comparison,
  validity + aspect + detail + row_id + project_letter ~ phase,
  value.var = "rel_score"
)
protocol_comparison_wide_delta <- as.data.table(protocol_comparison_wide_delta)
setnames(protocol_comparison_wide_delta, "Exploratory Stage", "exploratory")
setnames(protocol_comparison_wide_delta, "Confirmatory Stage", "confirmatory")

# eps <- 0.25
# protocol_comparison_wide_ratio[, log2_fc := log2((confirmatory + eps) / (exploratory + eps))]


protocol_comparison_wide_delta[, diff := confirmatory - exploratory]

protocol_comparison_wide_delta[, diff_scaled :=
                (diff - min(diff, na.rm = TRUE)) /
                (max(diff, na.rm = TRUE) - min(diff, na.rm = TRUE))
]


heatmap_delta <- ggplot(
  protocol_comparison_wide_delta,
  aes(
    x = project_letter,
    y = row_id,
    fill = diff
  )
) +
  geom_tile(color = "white", linewidth = 0.15, width = 0.95, height = 1.2) +
  facet_grid(
    scales = "free_y",
    space = "fixed"
  ) +
  scale_fill_gradient2(
    low  = "#4575B4",   # decreased validity
    mid  = "bisque2",   # no change
    high = "#D73027",   # increased validity
    midpoint = 0,       
    limits = c(-1, 1),
    breaks = c(-1, 0, 1)
  ) +   
  labs(x = "Confirmatory Studies", y = NULL, fill = "delta score \n(Confirmatory - \nExploratory)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.title.x = element_text(size = 22),
    axis.text.y = element_text(size = 15),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 24, face = "bold"),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22, margin = margin(b = 8))
  )

heatmap_delta

ggsave(
  filename = file.path(save_dir_decide, "protocol_comparison_heatmap_delta.png"),
  plot = heatmap_delta,
  width = 15,   
  height = 9,  
  dpi = 300,
  bg = "white"
)

# simplified plot change from exploratory to confirmatory: delta confirmati --------

protocol_comparison_wide_delta[, validity_aspect := paste(aspect, validity, sep = " | ")]

protocol_comparison_delta_simplified <- protocol_comparison_wide_delta[, 
                    .(
                      delta_mean = mean(diff, na.rm = TRUE)),
                    by = .(validity, aspect, project_letter)]

protocol_comparison_delta_simplified[
  ,
  validity_aspect := paste(aspect, validity, sep = " | ")
]

# Set row levels
row_levels <- protocol_comparison_delta_simplified[
  order(
    factor(validity, levels = validity_levels),
    aspect
  ),
  unique(validity_aspect)
]
protocol_comparison_delta_simplified[, validity_aspect := factor(validity_aspect, levels = rev(row_levels))] # rev because ggplot inverts the order in the heatmap



heatmap_delta_simplified <- ggplot(
  protocol_comparison_delta_simplified,
  aes(
    x = project_letter,
    y = validity_aspect,
    fill = delta_mean
  )
) +
  geom_tile(color = "white", linewidth = 0.15, width = 0.95, height = 1.2) +
  facet_grid(
    rows = vars(validity_aspect),
    scales = "free_y",
    space = "fixed"
  ) +
  # scale_fill_gradient2(
  #   low  = "#D6F4F6",
  #   mid  = "#00AFBB",
  #   high = "#008A93",
  #   midpoint = 0.5,
  #   limits = c(0, 1),
  #   breaks = c(0, 1)
  # ) +
  scale_fill_gradient2(
    low  = "#4575B4",   # decreased validity
    mid  = "bisque2",   # no change
    high = "#D73027",   # increased validity
    midpoint = 0,
    limits = c(-1, 1),
    breaks = c(-1, 0, 1)
  ) +
  labs(x = "Confirmatory Studies", y = NULL, fill = "delta score \n(Confirmatory - \nExploratory)") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 22, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 22, face = "bold", margin = margin(t = 15)),
    panel.grid = element_blank(),
    strip.text.x = element_text(size = 24, face = "bold"),
    strip.text.y = element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.text = element_text(size = 22),
    legend.title = element_text(size = 22, margin = margin(b = 8))
  )

heatmap_delta_simplified

ggsave(
  filename = file.path(save_dir_decide, "protocol_comparison_heatmap_delta_simplified.png"),
  plot = heatmap_delta_simplified,
  width = 15,   
  height = 9,  
  dpi = 300,
  bg = "white"
)

# Radar plot --------------------------------------------------------------

# sum up all validities

protocol_comparison_radar <- protocol_comparison_simplified[, .(validity_mean = mean(rel_score_mean)), by = .(phase, validity)]

df_radar <- dcast(protocol_comparison_radar, 
                  phase ~ validity,
                  value.vat = "validity_mean")
setnames(df_radar, "phase", "group")  # required name for ggradar
setDT(df_radar)
# df_radar[, group := as.character(group)]
# df_radar[group == "Exploratory Stage", group := "Exploratory"]
# df_radar[group == "Confirmatory Stage", group := "Confirmatory"]

# set order
df_radar[, group := factor(group, levels = c("Exploratory Stage", "Confirmatory Stage"))]

df_radar <- as_tibble(df_radar)

radar_plot <- ggradar(
  df_radar,
  grid.min = 0,
  grid.mid = 0.5,
  grid.max = 1,
  # values.radar = c("0", "60", "130"),
  group.colours = c("#99E0E5", "#00AFBB"),
  group.line.width = 1.2,
  group.point.size = 3,
  axis.label.size = 5.5,
  # values.radar.size = 8,
  legend.text.size = 4.2
) + 
  theme(
  legend.text  = element_text(size = 12, face = "bold"),
  # legend.title = element_text(size = 12, face = "bold")
)

radar_plot

ggsave(
  filename = file.path(save_dir_decide, "protocol_comparison_radar_plot.png"),
  plot = radar_plot,
  width = 8,   
  height = 4,  
  dpi = 300,
  bg = "white"
)

saveRDS(radar_plot, file.path(save_dir_decide, "panels", "radar_plot.rds")) 

