#### Forest plots ####
confirmatory <- all_stats_wide[Stage == "confirmatory"]
# Label combined Project name and center for forest plot
# confirmatory[, label := paste0(Project_Name, "-", Center)]
#Anonymous center
confirmatory[, center_anon := paste0("Lab ", .GRP), by = .(Project_ID, Center)]
confirmatory[, label := paste0("P", Project_ID, " - ", center_anon)]
confirmatory <- confirmatory[order(project_letter)]

# Create row positions for each study
confirmatory[, row_position := .N:1]

# Identify project groups and their row positions
project_groups <- confirmatory[, .(
  start_row = min(row_position),
  end_row = max(row_position)
), by = project_letter]
project_groups <- project_groups[order(-start_row)]


n_studies <- nrow(confirmatory)

# Meta-analysis with fixed-effect model
meta_results <- all_stats_wide[Stage == "confirmatory", {
  res <- run_meta_fixed(hedges_g, se_g, Center)
  extract_meta_fixed(res)
}, by = .(Project_ID, Project_Name, Stage, project_letter)]

# Build list of metagen objects for forest plot — one per project
meta_res_list <- lapply(
  split(all_stats_wide[Stage == "confirmatory"], by = "project_letter"),
  function(dt) {
    res <- run_meta_fixed(dt$hedges_g, dt$se_g, dt$Center)
    update(res, common = TRUE, random = FALSE)
  }
)

# Combine all projects into one metagen object for forest plot
res_combined <- do.call(metabind, meta_res_list)

pdf(here("results", "forest_plots_hedges.pdf"),
    width = 8, height = max(8, n_studies * 0.3))

# Draw initial forest plot
forest(res_combined,
       xlim   = c(-6, 8),
       cex    = 0.8,
       header = TRUE,
       rows   = n_studies:1,
       addfit = FALSE,# removes the overall ES
       ylim   = c(-2, n_studies + 3))

# Add background shading for alternating projects
colors <- c("#F5F5F5", "#FFFFFF")  # Light gray and white
for(i in 1:nrow(project_groups)) {
  rect(
    xleft = par("usr")[1], 
    xright = par("usr")[2],
    ybottom = project_groups$start_row[i] - 0.5,
    ytop = project_groups$end_row[i] + 0.5,
    col = colors[(i %% 2) + 1],
    border = NA
  )
}

# Redraw forest plot on top of shading
par(new = TRUE)
forest(res,
       xlim = c(-6, 8),
       cex = 0.8,
       header = TRUE,
       rows = n_studies:1,
       addfit = FALSE, # removes the overall ES
       ylim = c(-2, n_studies + 3))
dev.off()
