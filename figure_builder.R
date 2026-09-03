
library(here)

source(here("packages.R"))

library(grid)
library(png)
library(patchwork)
library(magick)

#### 0. Paths ####
# Assumes each panel script now ends with saveRDS(object, file.path(panels_dir,
# "<object_name>.rds")) -- see notes at duplicated names below for the
# simulation panels, which need a dataset suffix at save time.

panels_dir  <- file.path(save_dir_decide, "panels")
canvas_dir  <- file.path(save_dir_decide, "canvas_exports")
fig_out_dir <- file.path(save_dir_decide, "figures_final")

# dir.create(panels_dir,  recursive = TRUE, showWarnings = FALSE)
# dir.create(canvas_dir,  recursive = TRUE, showWarnings = FALSE)
# dir.create(fig_out_dir, recursive = TRUE, showWarnings = FALSE)

#### 1. Helpers####

# Load a saved ggplot/grob panel by object name
load_panel <- function(name) {
  path <- file.path(panels_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing panel file: ", path,
         " -- did the source script save this object via saveRDS()?")
  }
  readRDS(path)
}

# Load a Canvas-exported raster image (PNG) as a grob.
# Falls back to a labeled grey placeholder if the file hasn't been
# exported yet, so the script still runs end-to-end.
load_canvas_panel <- function(filename, label = filename) {
  path <- file.path(canvas_dir, filename)
  if (!file.exists(path)) {
    return(
      grid.grabExpr({
        grid.rect(gp = gpar(fill = "grey90", col = "grey60"))
        grid.text(paste0("MISSING CANVAS EXPORT:\n", label),
                  gp = gpar(fontsize = 10, col = "grey30"))
      })
    )
  }
  img <- magick::image_read(path)
  rasterGrob(as.raster(img), interpolate = TRUE)
}

# Save + report a patchwork figure
save_figure <- function(plot, name, width, height, dpi = 300) {
  out_path <- file.path(fig_out_dir, paste0(name, ".pdf"))
  ggsave(out_path, plot = plot, width = width, height = height,
         units = "in", dpi = dpi, device = cairo_pdf)
  message("Saved: ", out_path)
  invisible(out_path)
}

#### 2. Duplicate-name fix for multi-dataset simulation panels####
# `effect_sizes` and `pr_scatter` are regenerated once per dataset run in
# simulation_replication_criteria.R (Bonapersona / Carneiro / Rosso), all
# under the same object name. This script assumes the SOURCE script has been
# updated to save each run under a dataset-suffixed filename, e.g.:
#
#   saveRDS(effect_sizes, file.path(panels_dir,
#           paste0("effect_sizes_", dataset_name, ".rds")))
#   saveRDS(pr_scatter, file.path(panels_dir,
#           paste0("pr_scatter_", dataset_name, ".rds")))
#
# where dataset_name is set at the top of each of the 3 runs
# ("bonapersona", "carneiro", "rosso"). The main-text panels use
# the Bonapersona run as the primary dataset; Carneiro/Rosso are
# supplementary only.

#### 0b. Shared harmonization theme ####
harmonize_theme <- theme(
  axis.title      = element_text(size = 12),
  axis.text       = element_text(size = 11),
  legend.text     = element_text(size = 11),
  legend.title    = element_text(size = 12, face = "bold"),
  legend.key.size = unit(0.35, "cm"),
  plot.title      = element_text(size = 16, face = "bold", hjust = 0),
  plot.title.position = "plot"
)

#### FIGURE 1####
# a) project timeline
timeline <- image_read(here("results/canvas_exports/timeline.png"))
# timeline_trimmed <- image_trim(timeline)                          # remove existing whitespace
# timeline_padded  <- image_border(timeline_trimmed, "#FFFFFF", "200x200")  # add uniform white padding back in
# image_write(timeline_padded, here("results/canvas_exports/timeline_padded.png"))

fig1a <- load_canvas_panel("timeline.png", "Figure 1a - project timeline")
# b) scatter plot of the exploratory vs pooled ES in the 2 datasets
fig1b <- load_panel("final_scatter_plot")                # confirmatory_vs_retrospective.R
# c) pooled g decide dt
fig1c <- load_panel("hedges_pooled_by_project")         # confirmatory_replication_assessment.R
# d) ES deconstruction decide dt
fig1d <- load_panel("p_shapley_signed")                          # confirmatory_replication_assessment.R
# e) pooled g retrospective dt
fig1e <- load_panel("hedges_pooled_by_project_ext")      # confirmatory_vs_retrospective.R
# f) ES deconstruction retrospective dt
fig1f <- load_panel("p_shapley_ext_signed")                      # confirmatory_vs_retrospective.R

fig1a <- wrap_elements(full = fig1a) + theme(plot.margin = margin(1, 1, 1, 1))
fig1b <- fig1b + theme(plot.margin = margin(5, 5, 5, 5))
fig1c <- fig1c + theme(plot.margin = margin(20, 20, 20, 20))
fig1d <- fig1d + theme(plot.margin = margin(20, 20, 20, 20))
fig1e <- fig1e + theme(plot.margin = margin(20, 20, 20, 20))
fig1f <- fig1f + theme(plot.margin = margin(20, 20, 20, 20))

# fig1c <- fig1c + coord_fixed(ratio = 0.4)
# fig1e <- fig1e + coord_fixed(ratio = 0.4)

design <- "
AB
CD
EF
"
figure1 <- wrap_plots(
  A = fig1a,
  B = fig1b,
  C = fig1c,
  D = fig1d,
  E = fig1e,
  F = fig1f,
  design = design
) +
  plot_layout(heights = c(4, 5, 5), widths = c(0.8, 1.2)) +
  plot_annotation(tag_levels = "A",
                  title = "Figure 1",
                  theme = harmonize_theme) &
  harmonize_theme &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figure1, "Figure1", width = 12, height = 14)


#### FIGURE 2####
# a) Di distance method decide
fig2a <- load_panel("p5_distance_decide")     # confirmatory_vs_retrospective.R
# b) Di distance method retrospective
fig2b <- load_panel("p5_ext")                 # confirmatory_vs_retrospective.R
# c) RMSE plot
fig2c <- load_panel("rmse_ext_decide_plot")   # confirmatory_vs_retrospective.R

fig2a <- fig2a + theme(plot.margin = margin(20, 20, 20, 20))
fig2b <- fig2b + theme(plot.margin = margin(20, 20, 20, 20))
fig2c <- fig2c + theme(plot.margin = margin(20, 20, 20, 20))

design <- "
AB
C#
##
"

figure2 <- wrap_plots(
  A = fig2a,
  B = fig2b,
  C = fig2c,
  design = design
) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A",
                  title = "Figure 2",
                  theme = harmonize_theme) &
  harmonize_theme &
  theme(plot.tag = element_text(size = 10, face = "bold")
        )

figure2[[1]] <- figure2[[1]] + theme(legend.title = element_blank())
figure2[[2]] <- figure2[[2]] + theme(legend.title = element_blank())

save_figure(figure2, "Figure2", width = 12, height = 14)

#### FIGURE 3####
# a) z-curve
fig3a <- load_panel("zcurve_plot")                       # confirmatory_vs_retrospective.R
# b) minimal table - replication criteria
fig3b <- load_canvas_panel("replication_mini_table.png",
                           "Figure 3b - replication criteria table")
# c) FPR simulation bar plots
fig3c <- load_panel("fpr_plot_pooled_n_bonapersona_2021")                  # simulation_replication_criteria.R
# d) PR plot
fig3d <- load_panel("pr_scatter_bonapersona_2021")             # simulation_replication_criteria.R (primary dataset)
# e) replication flags on the 2 dataset
# fig3e <- load_panel("combined_heatmap_plot")              
fig3e <- wrap_elements(full = load_panel("combined_heatmap_plot"))     # confirmatory_vs_retrospective.R


fig3a <- wrap_elements(full = 
                         ggplot() + theme_void() +
                         inset_element(fig3a, left = 0.08, bottom = 0.08, right = 0.98, top = 0.92)
)
fig3b <- wrap_elements(full = fig3b) + theme(plot.margin = margin(1, 1, 1, 1))
# fig3c <- wrap_elements(full = 
#                          ggplot() + theme_void() +
#                          inset_element(fig3c, left = 0.08, bottom = 0.08, right = 0.92, top = 0.92)
#                        )
fig3d <- fig3d + theme(plot.margin = margin(20, 20, 20, 20))
fig3e <- fig3e + theme(plot.margin = margin(20, 20, 20, 20))

design <- "
AB
CD
EE
"
figure3 <- wrap_plots(
  A = fig3a,
  B = fig3b,
  C = fig3c,
  D = fig3d,
  E = fig3e,
  design = design
) +
  plot_layout(heights = c(1.5, 1, 1), widths = c(0.5, 1.3)) +
  plot_annotation(tag_levels = "A",
                  title = "Figure 3",
                  theme = harmonize_theme) &
  harmonize_theme &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figure3, "Figure3", width = 12, height = 14)


#### FIGURE 4####
# NOTE: Fig4 uses gridExtra/grid instead of ggarrange+align="hv" 
# because align="hv" breaks cuses axis numbers across the page

# a) Smallest Detectable Effect bar plots # confirmatory_vs_retrospective.R
fig4a <- load_panel("combined_sde_decide_retrospective_plot") + theme(plot.margin = margin(25, 140, 40, 40))
# b) radar plot   # protocol_comparison.R
fig4b <- load_panel("radar_plot") 
# c) heatmap of all validity instances   # protocol_comparison.R
fig4c <- load_panel("heatmap_raw")
# d) mIV bar plots    # confirmatory_vs_retrospective.R
fig4d <- load_panel("combined_iv_decide_retrospective_plot") + theme(plot.margin = margin(20, 20, 20, 30))

figure4 <- ggarrange(
  fig4a, fig4b,
  fig4c, fig4d,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  widths = c(1.4, 1),
  heights = c(1.25, 1.25),
  align = "hv",
  top = text_grob("Figure 4", face = "bold", size = 16, x = 0, hjust = 0)
)

figure4 <- figure4 & harmonize_theme
figure4 <- figure4 & theme(plot.tag = element_text(size = 10, face = "bold"))

figure4 <- arrangeGrob(
  fig4a, fig4b,
  fig4c, fig4d,
  ncol = 2, nrow = 2,
  widths = c(1.4, 1),
  heights = c(1.25, 1.25)
)

pdf(file.path(fig_out_dir, "Figure4.pdf"), width = 13, height = 11)
grid.draw(figure4)
grid.text("Figure 4", x = unit(0.02, "npc"), y = unit(0.98, "npc"),
          just = c("left", "top"),
          gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()


#### FIGURE 5 ####
# Framework
figure5 <- wrap_elements(
  full = load_canvas_panel("decide_framework.png", "Figure 5 - Framework")) +
    plot_annotation(title = "Figure 5", theme = harmonize_theme)


save_figure(figure5, "Figure5", width = 10, height = 8)


#### FIGURE S1####

# a) CV bar plot
figS1a <- load_panel("cv_decide_vs_retrospective")                  # confirmatory_vs_retrospective.R
# b) g single labs decide
figS1b <- load_panel("hedges_effects_by_project")                   # confirmatory_replication_assessment.R
# c) g single labs retrospective
figS1c <- load_panel("hedges_effects_by_project_retrospective")     # confirmatory_vs_retrospective.R

design <- "
AB
C#
"

figureS1 <- wrap_plots(
  A = figS1a,
  B = figS1b,
  C = figS1c,
  design = design
) + 
  plot_annotation(tag_levels = "A",
                  title = "Figure S1",
                  theme = harmonize_theme) &
  harmonize_theme &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figureS1, "FigureS1", width = 12, height = 10)

#### FIGURE S2 ####
# Experimental Units   
# a) EU bar plot     # confirmatory_vs_retrospective.R
figureS2 <- load_panel("comb_decide_exernal_eu_plot") + 
  labs(title = "Figure S2") +
  harmonize_theme
save_figure(figureS2, "FigureS2", width = 13, height = 9)


#### FIGURE S3####

# unsupervised validities plot
# p_unsupervised_heatmap <- load_panel("p_unsupervised_heatmap")  # protocol_comparison.R

figureS3 <- wrap_elements(full = p_unsupervised_heatmap) +
  plot_annotation(title = "Figure S3", theme = harmonize_theme)

save_figure(figureS3, "FigureS3", width = 13, height = 9)

#### FIGURE S4####

# a) simulation graphical abstract
figS4a_raw <- load_canvas_panel("simulation_graphical_abstract.png",
                            "Figure S4a - simulation graphical abstract")
figS4a <- wrap_elements(full = 
                          ggplot() + theme_void() +
                          inset_element(figS4a_raw, left = 0, right = 1, bottom = 0.15, top = 0.85)
                          )

# b) heatmap simulation replication success
# figS4b <- load_panel("combined_heatmap_pooled_bonapersona_2021")

figS4b_raw <- load_panel("combined_heatmap_pooled_bonapersona_2021")
figS4b <- wrap_elements(full = 
                          ggplot() + theme_void() +
                          inset_element(figS4b_raw, left = 0, right = 1, bottom = 0.25, top = 0.75)
)

# c) simulation shrinkage sensitivity
# figS4c <- load_panel("shrinkage_sensitivity_bonapersona_2021")   # simulation_replication_criteria.R

figS4c_raw <- load_panel("shrinkage_sensitivity_bonapersona_2021") 
figS4c <- wrap_elements(full = 
                          ggplot() + theme_void() +
                          inset_element(figS4c_raw, left = 0, right = 1, bottom = 0.15, top = 0.85)
)


design <- "
AB
C#
"

figureS4 <- wrap_plots(
  A = figS4a,
  B = figS4b,
  C = figS4c,
  design = design
  ) +
  plot_layout(widths = c(1, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A",
                  title = "Figure S4",
                  theme = harmonize_theme) &
  harmonize_theme &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figureS4, "FigureS4", width = 15, height = 13)

#### FIGURE S5####
# a) empirical effect size distribution Bonapersona
figS5a <- load_panel("effect_sizes_bonapersona_2021") + theme(plot.margin = margin(10, 10, 10, 30))

# b) empirical effect size distribution Carneiro
figS5b <- load_panel("effect_sizes_carneiro_2018") + theme(plot.margin = margin(10, 10, 10, 30))

# c) empirical effect size distribution Rosso
figS5c <- load_panel("effect_sizes_rosso_2022") + theme(plot.margin = margin(10, 10, 10, 30))

row1 <- ggarrange(
  figS5a, figS5b, figS5c,
  ncol = 3, nrow = 1,
  labels = c("A", "B", "C"),
  align = "hv"
)

# d) PR plot Carneiro   # simulation_replication_criteria.R
figS5d <- load_panel("pr_scatter_carneiro_2018") + theme(plot.margin = margin(10, 10, 10, 20))

# e) PR plot Rosso     # simulation_replication_criteria.R
figS5e <- load_panel("pr_scatter_rosso_2022") + theme(plot.margin = margin(10, 10, 10, 20))

row2 <- ggarrange(
  figS5d, 
  figS5e,
  ncol = 2, nrow = 1,
  labels = c("D", "E"),
  align = "hv"
)

# f) PR plot Bonapersona by exploratory sample size   # simulation_replication_criteria.R
figS5f <- load_panel("p_pr_trajectories_bonapersona_2021") + theme(plot.margin = margin(10, 80, 10, 80))

row3 <- ggarrange(
  figS5f,  
  ncol = 1, nrow = 1,
  labels = "F"
)

figureS5 <- ggarrange(row1, row2, row3, nrow = 3, heights = c(1, 1.5, 1.5))
figureS5 <- figureS5 & harmonize_theme

figureS5 <- annotate_figure(
  figureS5,
  top = text_grob("Figure S5", face = "bold", size = 16, x = 0, hjust = 0)
)

save_figure(figureS5, "FigureS5", width = 15, height = 12)

