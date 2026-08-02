# -------------------------------------------------------------------------
#' 
#' Helper functions to fit to the CPE control data to estimate the number of
#' infectious virions per well
#' 
# -------------------------------------------------------------------------


# negative log-likelihood function for CPE VC fit -------------------------

nllh.CPE.vc <- function(data,par){
  # INPUTS:
  # data              CPE assay VC (virus control) data
  # par               parameter to be estimated: number of virions per well in the control samples
  # OUTPUT: negative log-likelihood for the fit to the CPE virus control data
  
  -sum(log(dbinom(data$n.neg,data$n.rep,exp(-par*data$dilution))))
}
