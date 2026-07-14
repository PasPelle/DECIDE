source(here("packages.R"))

# To run last after confimratory_vs_retrospective and simulation
# Check the sample size change from exploratory to confirmatory
multilab_scale_factors <- combined_sde[, .(
  median_n = median(n1_total + n2_total, na.rm = TRUE)
), by = .(dataset, stage)]

multilab_scale_factors <- dcast(multilab_scale_factors, dataset ~ stage, value.var = "median_n")
setDT(multilab_scale_factors)
multilab_scale_factors[, scale_factor := `Multi-lab` / Exploratory]

cat("Empirical multilab scale factors:\n")
print(multilab_scale_factors[, .(dataset, Exploratory, `Multi-lab`, scale_factor)])

# Use median across datasets for simulation
multilab_scale_factor <- median(multilab_scale_factors$scale_factor)
cat("\nUsing median scale factor for simulation:", round(multilab_scale_factor, 2), "\n")


# Comparing simulation with real datasets confirmatory and retrospective

#### Confirm that the CI in the simulation are much larger than the real data

# CI half-width from real confirmatory data
ci_hw_real <- rbind(
  meta_results_decide[, .(
    ci_halfwidth = (ci_upper - ci_lower) / 2,
    dataset = "DECIDE confirmatory",
    k = 2  # mostly 2 labs
  )],
  meta_results_ext[, .(
    ci_halfwidth = (ci_upper_confirmatory - ci_lower_confirmatory) / 2,
    dataset = "Retrospective multilab",
    k = NA
  )]
)

cat("Real data CI half-widths:\n")
print(ci_hw_real[, .(mean = mean(ci_halfwidth), 
                     median = median(ci_halfwidth),
                     min = min(ci_halfwidth),
                     max = max(ci_halfwidth)), 
                 by = dataset])

# CI half-width from simulation (single lab confirmatory)
cat("\nSimulation CI half-widths:\n")
simulation_results[!is.na(confirmatory_ci_lower), .(
  mean   = mean((confirmatory_ci_upper - confirmatory_ci_lower) / 2, na.rm = TRUE),
  median = median((confirmatory_ci_upper - confirmatory_ci_lower) / 2, na.rm = TRUE),
  min    = min((confirmatory_ci_upper - confirmatory_ci_lower) / 2, na.rm = TRUE),
  max    = max((confirmatory_ci_upper - confirmatory_ci_lower) / 2, na.rm = TRUE)
), by = exploratory_n]



# Compare total N between exploratory and confirmatory in real data
# The real data have lower SDE because of larger n at multilab stage
rbind(
  all_stats_wide[Stage == "exploratory", .(
    stage = "exploratory",
    total_n = mean(n1 + n2, na.rm = TRUE),
    median_n = median(n1 + n2, na.rm = TRUE)
  )],
  all_stats_wide[Stage == "confirmatory", .(
    stage = "confirmatory",
    total_n = mean(n1 + n2, na.rm = TRUE),
    median_n = median(n1 + n2, na.rm = TRUE)
  )]
)

# And for retrospective
rbind(
  retrospective_exploratory[, .(
    stage = "exploratory",
    total_n = mean(n1_exploratory + n2_exploratory, na.rm = TRUE),
    median_n = median(n1_exploratory + n2_exploratory, na.rm = TRUE)
  )],
  multi_lab_dt[, .(
    stage = "multilab",
    total_n = mean(n1_confirmatory + n2_confirmatory, na.rm = TRUE),
    median_n = median(n1_confirmatory + n2_confirmatory, na.rm = TRUE)
  )]
)
