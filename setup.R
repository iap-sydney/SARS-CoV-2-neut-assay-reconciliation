# -------------------------------------------------------------------------
#' 
#' This file is used to set up everything that's needed across the project.
#' It loads libraries and helper functions, sets defaults, and creates an output folder structure.
#' 
# -------------------------------------------------------------------------


# load packages -----------------------------------------------------------

library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(reshape2)
library(stringr)
library(MethComp)
library(censReg)
library(MASS)
library(writexl)
library(lme4)


# load local functions ----------------------------------------------------

source("helper-functions/01_General-helper-functions.R")
source("helper-functions/02_Data-processing.R")
source("helper-functions/03_CPE-helper-functions.R")
source("helper-functions/04_PV-helper-functions.R")
source("helper-functions/05_IS-helper-functions.R")


# set defaults ------------------------------------------------------------

alpha_CI <- 0.95 # confidence intervals (CIs) are 95% CIs

plot_width <- 16
plot_height <- 9


# create output folder & folder structure ---------------------------------

if(!dir.exists("output")){dir.create("output")}

output_subdirs = c('Figures','Tables','Fitting') # add all output sub-folders here

for (ii in 1:length(output_subdirs)) {
  if(!dir.exists(paste0("output/", output_subdirs[ii]))){
    dir.create(paste0("output/", output_subdirs[ii]))
  }
}

# clean up:
rm(output_subdirs,ii)
