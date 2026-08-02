# -------------------------------------------------------------------------
#' 
#' Code to reproduce analysis and figures presented in:
#' "Reconciliation of titer differences between SARS-CoV-2 neutralizing antibody assays"
#' E. Stadler, F. Mordant, M. Gartner, S. S. Docken, M. P. Davenport, D. S. Khoury, K. Subbarao
#' Journal of Virology, 2026
#' 
#' Code written by E. Stadler
#' 
# -------------------------------------------------------------------------


# analysis set-up ---------------------------------------------------------

### analysis:
standardize <- FALSE # standardize all titers (relative to NIBSC 21/338)
stand_sample <- 28 # sample used for standardization: NIBSC 21/338
stand_IU <- 4705 # standard sample IU/ml

### output control:
print_results <- TRUE # print figures, tables, and analysis results
save_results <- TRUE # save figures, tables, etc.


# setup -------------------------------------------------------------------

source("setup.R")


# processing --------------------------------------------------------------

#' Assays:
#' Cytopathic effect (CPE), Immunospot (IS), and Pseudovirus (PV) assays
#' 
#' CPE data overview:
#' original: CPE assay without removing inoculum
#' rem.low: CPE assay with removing inoculum and lower viral inoculum (same TCID50 as "original")
#' rem.high: CPE assay with removing inoculum & 8-times higher virus inoculum (compared to "rem.low" and "original")

# load, clean and format data for analysis:
source("processing/01_Assays-raw-data.R") # IS, CPE (original), PV assay
source("processing/02_CPE-removing-inoculum.R") # CPE assay with removing inoculum (both lower or higher virus inoculum)


# analysis ----------------------------------------------------------------

# Estimate neutralization titers (ND50s) for each assay:
source("analysis/01_IS-fit.R") 
source("analysis/02_CPE-fit.R") 
source("analysis/03_PV-fit.R") 

# Comparisons of the different assays: (run the above codes for estimating ND50s first)
source("analysis/04_Assay-comparison.R")

