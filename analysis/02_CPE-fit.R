# -------------------------------------------------------------------------
#' 
#' Visualize data and fit to all CPE assay data
#' 
# -------------------------------------------------------------------------


# analysis set-up ---------------------------------------------------------

# CPE data overview:
# original: CPE assay without removing inoculum
# rem.low: CPE assay with removing the inoculum and lower virus inoculum (same as "original")
# rem.high: CPE assay with removing the inoculum & 8-times higher virus inoculum (compared to "rem.low")

# specify CPE assays to analyse:
analyse.CPEs <- c("original","rem.low","rem.high")


# negative log-likelihood fit and visualization ---------------------------

# fit for all CPE assays to be fitted and visualization of the fit:
for(i in 1:length(analyse.CPEs)){
  # data for analysis and control data:
  if(analyse.CPEs[i]=="original"){
    CPE.ND50s <- CPE.original.ND50s
    CPE.vc <- CPE.original.vc
  }else if(analyse.CPEs[i]=="rem.low"){
    CPE.ND50s <- CPE.rem.low.ND50s
    CPE.vc <- CPE.rem.low.vc
  }else if(analyse.CPEs[i]=="rem.high"){
    CPE.ND50s <- CPE.rem.high.ND50s
    CPE.vc <- CPE.rem.high.vc
  }
  # save table with ND50s per assay run and for all assay runs:
  if(save_results){
    # table with estimated ND50's for each sample:
    write_xlsx(CPE.ND50s,glue::glue("output/Tables/CPE-ND50s(data-{analyse.CPEs[i]})_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
  }
  
  ### fit to the virus control data to estimate the infectious virions per well:
  tempfit.vc <- nlm(function(par){nllh.CPE.vc(data=CPE.vc,par=par)},c(4),iterlim=1e3,hessian = TRUE)
  n.virions.var <- solve(tempfit.vc$hessian) # variance of the estimated number of virions
  
  if(print_results){
    print(glue::glue("Least squares fit to all CPE control data ({analyse.CPEs[i]}):"))
    print(glue::glue("Estimated number of virions per well: {round(tempfit.vc$estimate,3)}"))
    print(glue::glue("Virions per well 95% CI: {round(tempfit.vc$estimate-qnorm((1+alpha_CI)/2)*sqrt(n.virions.var),3)} to {round(tempfit.vc$estimate+qnorm((1+alpha_CI)/2)*sqrt(n.virions.var),3)}"))
    print(glue::glue("Estimated TCID50 per well in the virus control samples: {round(tempfit.vc$estimate/log(2),3)}"))
  }
  if(save_results){
    save(tempfit.vc,file = glue::glue("output/Fitting/CPE-fit-controls(data-{analyse.CPEs[i]})_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.RData"))
  }
  
  # clean up:
  rm(CPE.vc,tempfit.vc,n.virions.var)
  
  
  ### linear mixed-effects model to estimate variability between assay runs:
  if(analyse.CPEs[i]=="original"){
    tmp.CPE.ND50s <- pivot_longer(CPE.ND50s[,names(CPE.ND50s)%in%c("sample","assay_run1","assay_run2")], cols = starts_with("assay"),
                                  names_to = "run", values_to = "ND50", names_pattern = "assay_run(.)")
    CPE.mod <- lmer(log10(ND50) ~ (1 | sample), data=tmp.CPE.ND50s, REML = TRUE) # model without random effect for run
    if(print_results){print(summary(CPE.mod))}
    run_SD <- sigma(CPE.mod) # residual SD only
    run_var_table <- data.frame(fold_diff_95perc = 10^((2)*run_SD), 
                                mean_run_fold_diff = mean(summary(pmax(CPE.ND50s$assay_run1/CPE.ND50s$assay_run2,CPE.ND50s$assay_run2/CPE.ND50s$assay_run1))))
    
    if(print_results){print(run_var_table)}
    if(save_results){
      write_xlsx(run_var_table,glue::glue("output/Tables/CPE-run-variability(data-{analyse.CPEs[i]})_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
    }
    
    # clean up:
    rm(tmp.CPE.ND50s,CPE.mod,run_SD,run_var_table)
  }
}

# clean up:
rm(i,CPE.ND50s)


# clean up ----------------------------------------------------------------

rm(analyse.CPEs)
