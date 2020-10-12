# ----------------- 0: load libraries -----------------
library(dplyr)
library(qgraph)
library(psychonetrics)
library(bootnet)

# data in "data_cole", "data_cole_covariates", and "data_cole_long"
# "cole" stands for "Cologne Chart of Life Events"

# ----------------- 1: estimate networks -----------------
# estimate network (without CTQ domains as covariates)
cole_nw <-
  estimateNetwork(
    data_cole,
    default = "EBICglasso",
    tuning = 0,
    corMethod = "spearman",
    missing = "pairwise"
  )

# estimate network (including CTQ domains as covariates)
cole_nw_covariates <-
  estimateNetwork(
    data_cole_covariates,
    default = "EBICglasso",
    tuning = 0,
    corMethod = "spearman",
    missing = "pairwise"
  )

# ----------------- 2: longitudinal modeling (panelgvar) -----------------
# note:
# PANSS G6: depression
# PANSS G3: guilt feelings
# burden of life events (details in paper)

# define design matrix (across 7 visits)
design <- matrix(
  c(
    "PANSS_G06_T0",
    "PANSS_G03_T0",
    "Bur_T0",
    
    
    "PANSS_G06_IV3",
    "PANSS_G03_IV3",
    "Bur_IV3",
    
    "PANSS_G06_IV6",
    "PANSS_G03_IV6",
    "Bur_IV6",
    
    "PANSS_G06_T1",
    "PANSS_G03_T1",
    "Bur_T1",
    
    
    "PANSS_G06_IV12",
    "PANSS_G03_IV12",
    "Bur_IV12",
    
    
    "PANSS_G06_IV15",
    "PANSS_G03_IV15",
    "Bur_IV15",
    
    "PANSS_G06_T2",
    "PANSS_G03_T2",
    "Bur_T2"
  ),
  ncol = 7,
  byrow = F
)

# define the panelgvar-model
Model <-
  panelgvar(
    data_cole_long,
    vars = design,
    estimator = "FIML",
    verbose = TRUE,
    storedata = TRUE
  )

# run the model, including modelsearch for optimal (final) model
Model_pruned_stepup <-
  Model %>% runmodel %>% modelsearch(addalpha = 0.05, prunealpha = 0.05)

# extract temporal network
temporal <- Model_pruned_stepup %>% getmatrix("PDC")

# bootstrap this analysis (1000 times)
# depending on number of vars & settings, this may take quite long
set.seed(1)
Bootstraps <- lapply(1:1000, function(x) {
  bootstrapped_model <-
    panelgvar(
      data_cole_long[sample(1:nrow(data_cole_long),
                            nrow(data_cole_long),
                            TRUE), ],
      vars = design,
      estimator = "FIML",
      verbose = TRUE,
      storedata = TRUE
    )  %>%
    runmodel %>%
    modelsearch(addalpha = 0.05, prunealpha = 0.05)
  return(bootstrapped_model)
})


# bootstrapped results: temporal
# check if individual edges are included (y/n)
resBoots_temp <-
  lapply(Bootstraps, function(x)
    ifelse(getmatrix(x, "PDC") > 0, 1, 0))

# how often (%) is edge included, over all bootstrap iterations
apply(simplify2array(resBoots_temp), 1:2, mean) 
