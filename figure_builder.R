source(here("packages.R"))

# Figure 1
# a) object: from Canvas (add placeholder); description: project timeline; script: NA
# b) object: final_scatter_plot; description: scatter plot of the exploratory vs pooled ES in the 2 datasets; script: confirmatory_vs_retrospective.R
# c) object: hedges_pooled_by_project; description: pooled g decide dt; script: confirmatory_replication_assessment.R
# d) object: p_shapley; description: ES deconstruction decide dt; script: confirmatory_replication_assessment.R
# e) object: hedges_pooled_by_project_ext description: pooled g retrospective dt; script: confirmatory_vs_retrospective.R
# f) object: p_shapley_ext; description: ES deconstruction retrospective dt; script: confirmatory_vs_retrospective.R

# Figure 2
# a) object: p5_distance_decide; description: Di distance method decide; script: confirmatory_vs_retrospective.R
# b) object: p5_ext; description: Di distance method retrospective; script: confirmatory_vs_retrospective.R
# c) object: rmse_ext_decide_plot; description: RMSE plot; script: confirmatory_vs_retrospective.R


# figure 3
# a) object: zcurve_plot; description: z-curve; script: confirmatory_vs_retrospective.R
# b) object: from Canvas (add placeholder) description: minimal table-replication criteria; script: NA
# c) object: fpr_plot_pooled_n; description: FPR simulation bar plots; script: simulation_replication_criteria.R
# d) object: pr_scatter; description: PR plot; script: simulation_replication_criteria.R
# e) object: combined_heatmap_plot; description: replication flags on the 2 dataset; script: confirmatory_vs_retrospective.R

# Figure 4
# a) object: combined_sde_decide_retrospective_plot; description: Smallest Detectable Effect bar plots; script: confirmatory_vs_retrospective.R
# b) object: radar_plot; description: radar plot; script: protocol_comparison.R
# c) object: heatmap_raw; description: heatmap of all validity istances; script: protocol_comparison.R
# d) object: combined_iv_decide_retrospective_plot; description: mIV bar plots; script: confirmatory_vs_retrospective.R

# Figure 5: from Canvas (add placeholder); description: Framework; script: NA

# Figure S1
# object: NA ; description: unsupervised validities plot; script: protocol_comparison.R
wrap_elements(full = p_unsupervised_heatmap)

# Figure S2
# a) object: comb_decide_exernal_eu_plot; description: EU bar plot; script: confirmatory_vs_retrospective.R
# b) object: cv_decide_vs_retrospective; description: CV bar plot; script: confirmatory_vs_retrospective.R
# c) object: hedges_effects_by_projecthedges_effects_by_project; description: g single labs decide; script: confirmatory_replication_assessment.R
# d) object: hedges_effects_by_project_retrospective; description: g single labs retrospective; script: confirmatory_vs_retrospective.R

# Figure S3
# a) object: from Canvas (add placeholder) description: simulation graphical abstract; script: NA
# b) object: combined_heatmap; description: heatmap simulation replication success; script: simulation_replication_criteria.R
# c) object: sensitivity_plot; description: simulation shrinkage sensitivity; script: simulation_replication_criteria.R

# Figure S4
# a) object: effect_sizes; description: empirical effect size distribution Bonapersona; script: simulation_replication_criteria.R
# b) object: effect_sizes; description: empirical effect size distribution Carneiro; script: simulation_replication_criteria.R
# c) object: effect_sizes; description: empirical effect size distribution Rosso; script: simulation_replication_criteria.R
# d) object: pr_scatter; description: PR plot Carneiro; script: simulation_replication_criteria.R
# e) object: pr_scatter; description: PR plot Rosso; script: simulation_replication_criteria.R
# f) object: p_pr_trajectories; description: PR plot Bonapersona by exploratory sample size; script: simulation_replication_criteria.R


# a) object: description:
# b) object: description:
# c) object: description:
# d) object: description:
# e) object: description:

source(here("packages.R"))
library(grid)
library(png)
library(patchwork)

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
# ("bonapersona", "carneiro", "rosso"). The main-text panels (Fig 3d) use
# the Bonapersona run as the primary dataset; Carneiro/Rosso are
# supplementary only (Fig S4).

#### FIGURE 1####
# a) project timeline
fig1a <- load_canvas_panel("timeline.jpg", "Figure 1a - project timeline")
# b) scatter plot of the exploratory vs pooled ES in the 2 datasets
fig1b <- load_panel("final_scatter_plot")                # confirmatory_vs_retrospective.R
# c) pooled g decide dt
fig1c <- load_panel("hedges_pooled_by_project")           # confirmatory_replication_assessment.R
# d) ES deconstruction decide dt
fig1d <- load_panel("p_shapley")                          # confirmatory_replication_assessment.R
# e) pooled g retrospective dt
fig1e <- load_panel("hedges_pooled_by_project_ext")       # confirmatory_vs_retrospective.R
# f) ES deconstruction retrospective dt
fig1f <- load_panel("p_shapley_ext")                      # confirmatory_vs_retrospective.R

figure1 <- (wrap_elements(full = fig1a) | fig1b | fig1c) /
  (fig1d | fig1e | fig1f) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figure1, "Figure1", width = 16, height = 10)

#### FIGURE 2####
# a) Di distance method decide
fig2a <- load_panel("p5_distance_decide")     # confirmatory_vs_retrospective.R
# b) Di distance method retrospective
fig2b <- load_panel("p5_ext")                 # confirmatory_vs_retrospective.R
# c) RMSE plot
fig2c <- load_panel("rmse_ext_decide_plot")   # confirmatory_vs_retrospective.R

figure2 <- (fig2a | fig2b | fig2c) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figure2, "Figure2", width = 15, height = 6)

#### FIGURE 3####
# a) z-curve
fig3a <- load_panel("zcurve_plot")                       # confirmatory_vs_retrospective.R
# b) minimal table - replication criteria
fig3b <- load_canvas_panel("figure3_replication_criteria_table.png",
                           "Figure 3b - replication criteria table")
# c) FPR simulation bar plots
fig3c <- load_panel("fpr_plot_pooled_n_bonapersona_2021")                  # simulation_replication_criteria.R
# d) PR plot
fig3d <- load_panel("pr_scatter_bonapersona_2021")             # simulation_replication_criteria.R (primary dataset)
# e) replication flags on the 2 dataset
fig3e <- load_panel("combined_heatmap_pooled_bonapersona_2021")              # confirmatory_vs_retrospective.R

figure3 <- (fig3a | wrap_elements(full = fig3b) | fig3c) /
  (fig3d | fig3e) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figure3, "Figure3", width = 15, height = 10)

#### FIGURE 4####
# a) Smallest Detectable Effect bar plots
fig4a <- load_panel("combined_sde_decide_retrospective_plot")  # confirmatory_vs_retrospective.R
# b) radar plot
fig4b <- load_panel("radar_plot")                              # protocol_comparison.R
# c) heatmap of all validity instances
fig4c <- load_panel("heatmap_raw")                              # protocol_comparison.R
# d) mIV bar plots
fig4d <- load_panel("combined_iv_decide_retrospective_plot")    # confirmatory_vs_retrospective.R

figure4 <- (fig4a | fig4b) /
  (fig4c | fig4d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figure4, "Figure4", width = 12, height = 10)

#### FIGURE 5 (single Canvas panel)####
# Framework
figure5 <- wrap_elements(
  full = load_canvas_panel("decide_framework.png", "Figure 5 - Framework")
)

save_figure(figure5, "Figure5", width = 10, height = 8)

#### FIGURE S1###
# unsupervised validities plot
p_unsupervised_heatmap <- load_panel("p_unsupervised_heatmap")  # protocol_comparison.R

figureS1 <- wrap_elements(full = p_unsupervised_heatmap)

save_figure(figureS1, "FigureS1", width = 13, height = 9)

#### FIGURE S2####
# a) EU bar plot
figS2a <- load_panel("comb_decide_exernal_eu_plot")                 # confirmatory_vs_retrospective.R
# b) CV bar plot
figS2b <- load_panel("cv_decide_vs_retrospective")                  # confirmatory_vs_retrospective.R
# c) g single labs decide
figS2c <- load_panel("hedges_effects_by_project")                   # confirmatory_replication_assessment.R
# d) g single labs retrospective
figS2d <- load_panel("hedges_effects_by_project_retrospective")     # confirmatory_vs_retrospective.R

figureS2 <- (figS2a | figS2b) /
  (figS2c | figS2d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figureS2, "FigureS2", width = 12, height = 10)

##### FIGURE S3####
# a) simulation graphical abstract
figS3a <- load_canvas_panel("simulation_graphical_abstract.png",
                            "Figure S3a - simulation graphical abstract")
# b) heatmap simulation replication success
figS3b <- load_panel("combined_heatmap_pooled_bonapersona_2021")   # simulation_replication_criteria.R
# c) simulation shrinkage sensitivity
figS3c <- load_panel("shrinkage_sensitivity_bonapersona_2021")   # simulation_replication_criteria.R

figureS3 <- (wrap_elements(full = figS3a) | figS3b | figS3c) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figureS3, "FigureS3", width = 15, height = 6)

#### FIGURE S4####
# a) empirical effect size distribution Bonapersona
figS4a <- load_panel("effect_sizes_bonapersona_2021")  # simulation_replication_criteria.R
# b) empirical effect size distribution Carneiro
figS4b <- load_panel("effect_sizes_carneiro_2018")     # simulation_replication_criteria.R
# c) empirical effect size distribution Rosso
figS4c <- load_panel("effect_sizes_rosso_2022")        # simulation_replication_criteria.R
# d) PR plot Carneiro
figS4d <- load_panel("pr_scatter_carneiro_2018")       # simulation_replication_criteria.R
# e) PR plot Rosso
figS4e <- load_panel("pr_scatter_rosso_2022")          # simulation_replication_criteria.R
# f) PR plot Bonapersona by exploratory sample size
figS4f <- load_panel("p_pr_trajectories_bonapersona_2021")         # simulation_replication_criteria.R (Bonapersona, by exploratory n)

figureS4 <- (figS4a | figS4b | figS4c) /
  (figS4d | figS4e | figS4f) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 10, face = "bold"))

save_figure(figureS4, "FigureS4", width = 15, height = 10)


# a) object: description:
# b) object: description:
# c) object: description:
# d) object: description:
# e) object: description: