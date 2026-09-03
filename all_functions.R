####### Functions  #####

source(here("packages.R"))

# Helper functions ----

# Hedges correction factor J
hedges_J <- function(n1, n2) {
  df <- (as.numeric(n1) + as.numeric(n2) - 2)
  # J(df) = 1 - 3/(4df - 1)  ==  1 - 3/(4(n1+n2) - 9)
  1 - (3 / (4 * df - 1))
}

# convert Cohen's d to Hedges' g 
d_to_g <- function(d, n1, n2) {
  as.numeric(d) * hedges_J(n1, n2)
}

# Small Telescope helper function: compute d33 for a two-sample t-test ----
compute_d33 <- function(n1, n2, alpha = 0.05, power = 0.33) {
  n1 <- as.numeric(n1); n2 <- as.numeric(n2)
  n_total <- pmin(n1, n2) #get the minimum n if sample size is disbalanced
  out <- vapply(n_total, function(n) {
    if (!is.finite(n) || n <= 1) return(NA_real_)
    pwr::pwr.t.test(
      n = n,
      sig.level = alpha,
      power = power,
      type = "two.sample",
      alternative = "two.sided"
    )$d
  }, numeric(1))
  out
}




# default operator: use the default y if x doesn't exist
`%||%` <- function(x, y) if (!is.null(x)) x else y 


# Simulation: compute SDE as Hedges' g for a given n per group (two-sample t-test)
compute_sde_g <- function(n_per_group, alpha = 0.05, power = 0.8) {
  n <- as.numeric(n_per_group)
  
  d <- vapply(n, function(nn) {
    if (!is.finite(nn) || nn <= 1) return(NA_real_)
    pwr::pwr.t.test(
      n = nn,
      sig.level = alpha,
      power = power,
      type = "two.sample",
      alternative = "two.sided"
    )$d
  }, numeric(1))
  
  # Convert d to g
  d_to_g(d, n1 = n, n2 = n)
}

# Main fx: compute replication methods from a harmonized table -----------------------------------------------------------------------

compute_replication_flags <- function(
    dt,
    # define params that will be used for the replication methods
    cfg = list(alpha = 0.05, z_pi = 1.96, power_d33 = 0.33, power_sde = 0.8,
               sceptical_type = "golden", sceptical_alternative = "two.sided"),
    # Column mapping (so each dataset can pass its own names without renaming upstream)
    # cols is the named list that maps keys and actual column names
    cols = list(
      id = "id",
      g_exp = "g_exp", se_exp = "se_exp", ci_lo_exp = "ci_lo_exp", ci_hi_exp = "ci_hi_exp", p_exp = "p_exp",
      n1_exp = "n1_exp", n2_exp = "n2_exp",
      g_conf = "g_conf", se_conf = "se_conf", ci_lo_conf = "ci_lo_conf", ci_hi_conf = "ci_hi_conf", p_conf = "p_conf",
      # tau2_conf = "tau2_conf",        # for PI: it doesn't apply to the simulation (fixed-effect model)
      SDE_pooled = "SDE_pooled",        # for SDE criterion: it doesn't apply to the simulation (monolab)
      SESOI_g = "SESOI_g"               # for SESOI criterion: it doesn't apply to the retrospective dataset
    ),
    compute_thresholds = TRUE,          # if TRUE, compute d33/g33 from n1_exp/n2_exp when not present
    d33_col = "d33",                    
    g33_col = "g33"                     
) {
  stopifnot(is.data.table(dt))
  
  # Convenience getters: getc returns the real column name from the key
  getc <- function(name) cols[[name]]
  
  # helper that knows what “exists” means
  has_col <- function(key) {
    nm <- getc(key)
    is.character(nm) && length(nm) == 1 && nm %in% names(dt)
  }
  
  # Required columns for the shared methods
  req <- c(getc("g_exp"), getc("g_conf"), getc("ci_lo_conf"), getc("ci_hi_conf"), getc("p_exp"), getc("p_conf"))
  # Missing columns
  missing_req <- req[!req %in% names(dt)]
  if (length(missing_req) > 0) {
    stop("Missing required columns in dt: ", paste(missing_req, collapse = ", "))
  }
  
  alpha <- cfg$alpha %||% 0.05  # use alpha otherwise default it to 0.05
  
  # One-sided critical value for superiority tests
  z_sup <- qnorm(1 - alpha)
  
  # One-sided bounds for confirmatory effect
  dt[, `:=`(
    conf_lower_1s = get(getc("g_conf")) - z_sup * get(getc("se_conf")),
    conf_upper_1s = get(getc("g_conf")) + z_sup * get(getc("se_conf"))
  )]
  
  ## Start computing replication methods ## 
  #  === Direction agreement ===
  dt[, direction_agreement := sign(get(getc("g_exp"))) == sign(get(getc("g_conf")))]
  
  #  === Exploratory g in confirmatory CI (+ direction) ===
  dt[, ci_agreement :=
       get(getc("g_exp")) >= get(getc("ci_lo_conf")) &
       get(getc("g_exp")) <= get(getc("ci_hi_conf")) &
       direction_agreement
  ]
  
  # # === CI overlap (+ direction) ===
  # # Only if exploratory CI exists; otherwise set NA
  # dt[, ci_overlap :=
  #        !(get(getc("ci_hi_exp")) < get(getc("ci_lo_conf")) |
  #            get(getc("ci_hi_conf")) < get(getc("ci_lo_exp"))) &
  #        direction_agreement
  #   ]
  
  # === Confirmatory significant + same direction ===
  dt[, ttest_sig :=
       # get(getc("p_exp")) < alpha &
       get(getc("p_conf")) < alpha &
       direction_agreement
  ]
  
  # # === Prediction interval containment (+ direction) ===
  # # PI only meaningful if tau2_conf exists; otherwise set NA
  # if (has_col("tau2_conf")) {
  #   z_pi <- cfg$z_pi %||% 1.96
  #   dt[, `:=`(
  #     PI_lower_confirmatory = get(getc("g_conf")) - z_pi * sqrt(get(getc("se_conf"))^2 + get(getc("tau2_conf"))),
  #     PI_upper_confirmatory = get(getc("g_conf")) + z_pi * sqrt(get(getc("se_conf"))^2 + get(getc("tau2_conf")))
  #   )]
  #   dt[, exploratory_within_confirmatory_PI :=
  #        get(getc("g_exp")) >= PI_lower_confirmatory &
  #        get(getc("g_exp")) <= PI_upper_confirmatory &
  #        direction_agreement
  #   ]
  # } else {
  #   dt[, `:=`(
  #     PI_lower_confirmatory = NA_real_,
  #     PI_upper_confirmatory = NA_real_,
  #     exploratory_within_confirmatory_PI = NA
  #   )]
  # }
  
  # === Sceptical p-value ===
  dt[, `:=`(
      exploratory_z = get(getc("g_exp")) / get(getc("se_exp")),
      confirmatory_z = get(getc("g_conf")) / get(getc("se_conf")),
      variance_ratio = (get(getc("se_exp")) / get(getc("se_conf")))^2
    )]
  
  dt[, sceptical_p := ifelse(
      !is.finite(exploratory_z) | !is.finite(confirmatory_z) | !is.finite(variance_ratio) |
        get(getc("se_exp")) == 0 | get(getc("se_conf")) == 0,
      NA_real_,
      pSceptical(
        zo = exploratory_z,
        zr = confirmatory_z,
        c  = variance_ratio,
        alternative = cfg$sceptical_alternative %||% "two.sided",
        type = cfg$sceptical_type %||% "golden"
      )
    )]
    
  dt[, sceptical_sig := sceptical_p < alpha]
  
  # === Small telescopes ===
  
  # Compute d33 (Cohen's d threshold at 33% power)
  dt[, (d33_col) := compute_d33(
    n1 = get(getc("n1_exp")),
    n2 = get(getc("n2_exp")),
    alpha = alpha,
    power = cfg$power_d33 # %||% 0.33
  )]
  
  # Convert d33 -> g33 using Hedges correction based on actual (possibly unbalanced) n1+n2
  dt[, (g33_col) := d_to_g(
    d  = get(d33_col),
    n1 = get(getc("n1_exp")),
    n2 = get(getc("n2_exp"))
  )]
  
  # # Apply the CI-superiority small-telescope criterion
  # dt[, small_telescope_confirmed := direction_agreement & (
  #   (get(getc("g_exp")) > 0 & conf_lower_1s >  get(g33_col)) |
  #     (get(getc("g_exp")) < 0 & conf_upper_1s < -get(g33_col))
  # )]
  
  # SUCCESS = We cannot reject that the effect is as large as g33.
  # This happens if the one-sided CI of the replication includes or exceeds g33.
  dt[, small_telescope_success := TRUE] # Default to success
  
  dt[(get(getc("g_exp")) > 0), small_telescope_success := 
       !(conf_upper_1s < get(g33_col))]
  
  dt[(get(getc("g_exp")) < 0), small_telescope_success := 
       !(conf_lower_1s > -get(g33_col))]
  
  # Final flag should also respect direction agreement
  dt[, small_telescope_confirmed := small_telescope_success & direction_agreement]
  
  # === Smallest Detectable Effect (SDE) ===
  # point estimate exceeds SESOI (liberal)
  dt[, sde_confirmed := direction_agreement & (
    (get(getc("g_exp")) > 0 & get(getc("g_conf")) >=  get(getc("SDE_pooled"))) |
      (get(getc("g_exp")) < 0 & get(getc("g_conf")) <= -get(getc("SDE_pooled")))
  )]
  
  # CI one-sided superiority (conservative)
  # dt[, sde_confirmed := direction_agreement & (
  #     (get(getc("g_exp")) > 0 & conf_lower_1s > get(getc("SDE_pooled"))) |
  #       (get(getc("g_exp")) < 0 & conf_upper_1s < -get(getc("SDE_pooled")))
  # )]
  
  #  === Smallest Effect Size of Interest (SESOI) ===
  # SESOI: CI beyond SESOI (+ direction) (optional) ----
  
  if (has_col("SESOI_g")) {
    
    # point estimate exceeds SESOI (liberal)
    dt[, ci_above_sesoi := direction_agreement & (
      (get(getc("g_exp")) > 0 & get(getc("g_conf")) >=  get(getc("SESOI_g"))) |
        (get(getc("g_exp")) < 0 & get(getc("g_conf")) <= -get(getc("SESOI_g")))
    )]
    
    # CI one-sided superiority (conservative)
    # dt[, ci_above_sesoi := direction_agreement & (
    #   (get(getc("g_exp")) > 0 & conf_lower_1s >= get(getc("SESOI_g"))) |
    #     (get(getc("g_exp")) < 0 & conf_upper_1s <= -get(getc("SESOI_g")))
    # )]
    
  } else {
    dt[, ci_above_sesoi := NA]
  }
  
  dt
}

# Function to make a pairwise Wilcoxon table, useful for p-value annotations with ggpubr package
make_pw_wilcox_stat_tbl <- function(data, y, group,
                                    p_adjust = "BH",
                                    step_frac = 0.08,
                                    digits = 3) {
  dt <- as.data.table(data)
  
  if (!y %chin% names(dt)) stop("Column `", y, "` not found.")
  if (!group %chin% names(dt)) stop("Column `", group, "` not found.")
  
  yy <- as.numeric(dt[[y]])
  gg <- droplevels(as.factor(dt[[group]]))
  
  ok <- !is.na(yy) & !is.na(gg)
  yy <- yy[ok]
  gg <- droplevels(gg[ok])
  
  lvl <- levels(gg)
  
  pw <- pairwise.wilcox.test(
    x = yy,
    g = gg,
    p.adjust.method = p_adjust
  )
  
  p_long <- as.data.table(as.table(pw$p.value))
  setnames(p_long, c("group2", "group1", "p_adj"))
  p_long <- p_long[!is.na(p_adj)]
  
  # Extract W and n per pair
  pair_stats <- p_long[, {
    x1 <- yy[gg == group1]
    x2 <- yy[gg == group2]
    wt <- wilcox.test(x1, x2)
    .(W = wt$statistic, n1 = length(x1), n2 = length(x2))
  }, by = .(group1, group2)]
  
  p_long <- pair_stats[p_long, on = .(group1, group2)]
  
  # Exact formatted p-values
  p_long[, p_label := paste0(
    "P = ",
    formatC(p_adj, format = "f", digits = digits)
  )]
  
  # Order brackets: shorter spans lower
  p_long[, i1 := match(group1, lvl)]
  p_long[, i2 := match(group2, lvl)]
  p_long[, span := abs(i2 - i1)]
  setorder(p_long, span, i1, i2)
  
  # y positions
  y_rng <- range(yy, na.rm = TRUE)
  y_span <- diff(y_rng)
  ymax <- y_rng[2]
  step <- step_frac * ifelse(y_span == 0, max(1, ymax), y_span)
  
  p_long[, y.position := ymax + step * seq_len(.N)]
  
  p_long[, .(
    group1,
    group2,
    W,
    n1,
    n2,
    p.adj = p_adj,
    p = p_adj,
    p_label,
    y.position
  )]
}

## Meta analysis per project and stage
# Function to avoid scoping issues inside data.table
run_meta_fixed <- function(te, se, studlab) {
  meta::metagen(
    TE = te,                  # Vector of ES
    seTE = se,                # Standard error of the ES
    studlab = studlab,        # Label of the study
    sm = "SMD",               # Standardized Mean Difference (these are continuous outcomes)
    
    # Turn off random effects and enable fixed effects (common)
    # method.tau = "REML",      # Restricted Maximum Likelihood: method to estimate between study variance
    # method.random.ci = "classic",  # Selected classic (z-distribution) instead of Hartung-Knapp (t-distribution) to have the same z scale as the replication assessment
    common = TRUE,           # Fixed-effect model 
    random = FALSE           # Random-effects model
  )
}

# Helper to extract fixed-effects results 
extract_meta_fixed <- function(res) {
  list(
    g_pooled = res$TE.common,
    se       = res$seTE.common,
    ci_lower = res$lower.common,
    ci_upper = res$upper.common,
    pval     = res$pval.common,
    I2       = res$I2,
    tau2     = res$tau2
  )
}

# Shared theme
shrinkage_theme <- function() {
  list(
    theme_prism(),
    theme(
      legend.position    = "none",
      strip.text         = element_text(size = 9),
      panel.grid.major.y = element_line(color = "gray90")
    )
  )
}

# Helper: base plot with shared structure
make_shrinkage_plot <- function(dt, x_var, project_var, title_label, color) {
  ggplot(dt, aes(x = .data[[x_var]], y = .data[[project_var]])) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = 1, linetype = "dotted", color = "gray40") +
    geom_point(size = 3, color = color) +
    facet_wrap(~ metric) +
    labs(x = NULL, y = NULL) +
    shrinkage_theme()
}

# Single helper used throughout — takes pooled SD as input
hedges_g_from_parts <- function(mu_ctrl, mu_treated, sd_pooled, n_ctrl, n_treated) {
  d <- (mu_ctrl - mu_treated) / sd_pooled
  J <- 1 - (3 / (4 * (n_ctrl + n_treated) - 9))
  d * J
}

# Z-curve function to reconstruct fitted densities to build a ggplot
get_fitted_density <- function(zfit, z_seq) {
  weights <- zfit$fit$weights
  mus     <- zfit$fit$mu
  density_vals <- rowSums(
    sapply(seq_along(mus), function(i) {
      weights[i] * dnorm(z_seq, mean = mus[i], sd = 1)
    })
  )
  data.table(z = z_seq, density = density_vals)
}

# Run the Shapley weighting formula on 8 coalition values
compute_shapley <- function(v0, vc, vt, vs, vct, vcs, vts, vcts) {
  phi_ctrl <- (2/6) * (vc - v0) +
    (1/6) * (vct - vt) +
    (1/6) * (vcs - vs) +
    (2/6) * (vcts - vts)
  
  phi_treated <- (2/6) * (vt - v0) +
    (1/6) * (vct - vc) +
    (1/6) * (vts - vs) +
    (2/6) * (vcts - vcs)
  
  phi_sd <- (2/6) * (vs - v0) +
    (1/6) * (vcs - vc) +
    (1/6) * (vts - vt) +
    (2/6) * (vcts - vct)
  
  list(phi_ctrl = -phi_ctrl, phi_treated = -phi_treated, phi_sd = -phi_sd)
}
# Functions for study simulation ------------------------------------------


# Function to simulate a single study and return all statistics
simulate_study_complete <- function(true_effect, n_per_group) {
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

# Functions for study simulation ------------------------------------------


# Function to simulate a single study and return all statistics
simulate_study_complete <- function(true_effect, n_per_group) {
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
  exploratory_results <- simulate_study_complete(true_effect, exploratory_n)
  
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
  
  # Calculate true effect for confirmatory study by applying the shrinkage factor:
  # Note: shrinkage factors won't be modeled here as they can have many sources like 
  # winner's curse, publication bias, flexible experimental methods/analyses.
  true_confirmatory_effect <- true_effect * (1 - shrinkage_factor)
  
  # Run confirmatory study
  confirmatory_results <- simulate_study_complete(true_confirmatory_effect, 
                                                  confirmatory_n)
  
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
      
      # Exploratory results
      exploratory_p = x$exploratory$p_value,
      exploratory_ss = x$exploratory$sample_size,
      exploratory_ci_lower = x$exploratory$ci_lower,
      exploratory_ci_upper = x$exploratory$ci_upper,
      exploratory_true_g = x$true_effect,  
      exploratory_observed_g = x$exploratory$observed_g,
      exploratory_se_g = x$exploratory$se_g,
      
      # Confirmatory results
      confirmatory_p = x$confirmatory$p_value,
      confirmatory_ss = x$confirmatory$sample_size,
      confirmatory_ci_lower = x$confirmatory$ci_lower,
      confirmatory_ci_upper = x$confirmatory$ci_upper,
      confirmatory_true_g = x$true_confirmatory_effect,
      confirmatory_observed_g = x$confirmatory$observed_g,
      confirmatory_se_g = x$confirmatory$se_g
    )
  })
  
  rbindlist(dt_list)
}


  
# Direct fixed-effects pooling — much faster than metagen, this is for the simulation
pool_fixed <- function(g_vals, se_vals) {
  # precision weight
  w        <- 1 / se_vals^2 
  
  # weighted mean and SE
  g_pool   <- sum(w * g_vals) / sum(w)
  se_pool  <- sqrt(1 / sum(w))
  
  z        <- g_pool / se_pool
  p_val    <- 2 * pnorm(-abs(z))
  list(
    observed_g = g_pool,
    se_g       = se_pool,
    p_value    = p_val,
    ci_lower   = g_pool - 1.96 * se_pool,
    ci_upper   = g_pool + 1.96 * se_pool
  )
}




# Full research trajectory with a multilab setting
# Confirmatory stage as k labs 
# sample size at confirmatory stage increases 2.8x like the decide and retrospective real data
run_research_trajectory_multilab <- function(true_effect, exploratory_n, shrinkage_factor,
                                    k_labs = 3, 
                                    multilab_scale_factor = 2.8) {
  
  # Run exploratory study — select only significant results (winner's curse)
  exploratory_results <- simulate_study_complete(true_effect, exploratory_n)
  
  if (exploratory_results$p_value >= alpha) {
    return(list(
      exploratory = exploratory_results,
      confirmatory = list(observed_g = NA, se_g = NA, p_value = NA,
                          ci_lower = NA, ci_upper = NA, sample_size = NA),
      confirmatory_ss = NA,
      true_confirmatory_effect = NA
    ))
  }
  
  observed_effect <- exploratory_results$observed_g
  if (abs(observed_effect) < 0.05) {
    return(list(
      exploratory = exploratory_results,
      confirmatory = list(observed_g = NA, se_g = NA, p_value = NA,
                          ci_lower = NA, ci_upper = NA, sample_size = NA),
      confirmatory_ss = NA,
      true_confirmatory_effect = NA
    ))
  }
  
  # Power-calculated n for single lab based on observed exploratory effect
  confirmatory_n_single <- tryCatch({
    ceiling(pwr.t.test(d = abs(observed_effect), power = 0.8,
                       sig.level = alpha, type = "two.sample")$n)
  }, error = function(e) NA)
  
  if (!is.na(confirmatory_n_single) && confirmatory_n_single > 50) {
    confirmatory_n_single <- 50
  }
  
  if (is.na(confirmatory_n_single) || confirmatory_n_single < 5) {
    return(list(
      exploratory = exploratory_results,
      confirmatory = list(observed_g = NA, se_g = NA, p_value = NA,
                          ci_lower = NA, ci_upper = NA, sample_size = NA),
      confirmatory_ss = NA,
      true_confirmatory_effect = NA
    ))
  }
  
  # Scale total N by empirical multi-lab factor and distribute across k labs
  # confirmatory_n_total <- ceiling(confirmatory_n_single * multilab_scale_factor)
  # n_per_lab <- max(5, floor(confirmatory_n_total / k_labs))
  
  # Alternative: use the same sample size per lab, minimum 5 to avoid unstable g and SE
  n_per_lab <- max(5, floor(confirmatory_n_single))
  
  # Apply shrinkage to true effect
  true_confirmatory_effect <- true_effect * (1 - shrinkage_factor)
  
  # Run k independent lab studies
  lab_results <- lapply(seq_len(k_labs), function(i) {
    simulate_study_complete(true_confirmatory_effect, n_per_lab)
  })
  
  g_vals  <- sapply(lab_results, `[[`, "observed_g")
  se_vals <- sapply(lab_results, `[[`, "se_g")
  
  # Pool with fixed-effects meta-analysis
  # res <- run_meta_fixed(te = g_vals, se = se_vals,
  #                       studlab = paste0("Lab", seq_len(k_labs)))
  # 
  # confirmatory_results <- list(
  #   observed_g  = res$TE.common,
  #   se_g        = res$seTE.common,
  #   p_value     = res$pval.common,
  #   ci_lower    = res$lower.common,
  #   ci_upper    = res$upper.common,
  #   sample_size = n_per_lab * k_labs
  # )
  
  # Pool with fixed-effects meta-analysis, using pool_fixed
  res <- pool_fixed(g_vals, se_vals)
  
  confirmatory_results <- list(
    observed_g  = res$observed_g,
    se_g        = res$se_g,
    p_value     = res$p_value,
    ci_lower    = res$ci_lower,
    ci_upper    = res$ci_upper,
    sample_size = n_per_lab * k_labs
  )
  
  return(list(
    exploratory              = exploratory_results,
    confirmatory             = confirmatory_results,
    confirmatory_ss          = n_per_lab * k_labs,
    true_confirmatory_effect = true_confirmatory_effect
  ))
}

