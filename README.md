## Reproducibility

This project uses [`renv`](https://rstudio.github.io/renv/) to manage package versions.

### Setup

1. Clone the repository
2. Open `DECIDE.Rproj` in RStudio
3. Install exact package versions:
```r
   renv::restore()
```
4. Place data files in the `/data` folder (see Data Availability below)
5. Run the full analysis pipeline:
```r
   source("main.R")
```

### R version
R 4.5.x — see `renv.lock` for exact package versions.



####################
#### SIMULATION ####
####################

GOAL: To stress-test different “replication methods” under realistic preclinical conditions, and see which ones are the most credible.

To recapitulate a preclinical trajectory I have included the following:

o	true effects come from a realistic empirical distribution (I have tested 3 datasets)
o	exploratory studies are small (n = 5, 10, 15, 20 per group),
o	confirmatory studies are powered on the observed exploratory effect,
o	and the true effect in the confirmatory phase is shrunk by some factor (s) δconf=δtrue⋅(1−s

Replication success rate is evaluated across:

o	effect size classes (small/medium/large true effects),
o	shrinkage levels (0, 20, 50, 80, 99%),
o	exploratory sample sizes.
o	False Positive Rate (FPR) when the confirmatory true effect is essentially zero ((s = 0.99)):

Precision–Recall / F1 by treating:

	Low shrinkage (s = 0 or 0.2) as “genuine replications” (ground truth = 1),
	Very high shrinkage (s=0.99) as “failed replications” (ground truth = 0),

Shrinkage sensitivity 
This is an indicator of the negative correlation between shrinkage level (s) and success. 
The objective is to test, as the true effect collapses towards 0, whether the criterion appropriately 
stop detecting successful replications, or it is blind to shrinkage.

The best methods will have: 

o	high success for genuinely stable effects,
o	low FPR when effects have shrunk
o	robustness to small exploratory (n).


