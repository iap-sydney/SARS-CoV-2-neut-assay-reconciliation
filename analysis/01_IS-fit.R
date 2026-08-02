# -------------------------------------------------------------------------
#' 
#' Visualize data and fit a curve to the IS assay data
#' 
# -------------------------------------------------------------------------


# specify analysis and outputs --------------------------------------------

### estimating ND50s:
# initial parameter values:
initial.slope.IS <- log(3) # slope = 3, the fitting function requires a log-transform
initial.ND50.IS <- log(1200) # ND50 = 1200, also log-transformed

### figures:
show_fit <- TRUE # visualize data with the fit (TRUE) or without (FALSE)


# fit to each assay run ---------------------------------------------------

# least-squares fit to the data to estimate the ND50 for each sample separately per assay run

### data
fit.data.IS <- IS.data[IS.data$plate<5,] # plate 5 contains only control data (negative samples) and is thus excluded from fitting
IS.ND50s <- CPE.original.ND50s[c(1:28),names(CPE.original.ND50s)%in%c("sample","sample_name")]
IS.ND50s$sample <- as.numeric(IS.ND50s$sample)

### fit for each assay run
for(run in unique(fit.data.IS$assay)){
  tmp.fit.data.IS <- fit.data.IS[fit.data.IS$assay==run,]
  
  # fit as above:
  n.plates <- 4 # total number of plates (4 per assay run)
  n.samples <- max(as.numeric(tmp.fit.data.IS$sample[!tmp.fit.data.IS$sample%in%c("VOC","NVC")]),na.rm = TRUE)
  
  ### fitting:
  # initial parameter values:
  par.init <- c(rep(mean(as.numeric(IS.data$spots[IS.data$sample=="VOC"])),length(unique(tmp.fit.data.IS$plate))*length(unique(tmp.fit.data.IS$assay))), # maximum number of spots per plate
                rep(initial.slope.IS,1), # slope for each of the 30 samples that are used in the fit
                rep(initial.ND50.IS,n.samples)) # ND50 for each sample
  # minimize least-squares:
  tempfit <- nlm(function(par){leastsq(data=tmp.fit.data.IS,par=par)},par.init,iterlim=1e3,hessian = TRUE)
  
  # make a table with estimated ND50's for each sample:
  tmp.IS.ND50s <- data.frame(sample=as.numeric(unique(tmp.fit.data.IS$sample[!tmp.fit.data.IS$sample%in%c("VOC","NVC")])),
                         ND50=exp(tempfit$estimate[n.plates+1+c(1:n.samples)]))
  
  # print & save parameter estimates:
  if(print_results){
    print(glue::glue("Least squares fit to IS data (assay run {run}):"))
    print(glue::glue("Sum of least squares: {tempfit$minimum}"))
    print(glue::glue("Exit flat: {tempfit$code}"))
    print(glue::glue("Iterations: {tempfit$iterations}"))
    print(glue::glue("AIC: {round(AIC_ls(tempfit$minimum,nrow(tmp.fit.data.IS),length(tempfit$estimate)),3)}"))
    print(glue::glue("Number of estimated parameters: {length(tempfit$estimate)}"))
    print(tmp.IS.ND50s)
    print(glue::glue("Slope (95% CI): {round(exp(tempfit$estimate[5]),2)} ({round(exp(tempfit$estimate[5]-qnorm((1+alpha_CI)/2)*sqrt(diag(solve(tempfit$hessian[5,5])))),2)} - {round(exp(tempfit$estimate[5]+qnorm((1+alpha_CI)/2)*sqrt(diag(solve(tempfit$hessian[5,5])))),2)})"))
  }
  if(save_results){
    # least squares fitting result:
    save(tempfit,file = glue::glue("output/Fitting/IS-assay{run}-fit_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.RData"))
  }

  # add to ND50s table:
  IS.ND50s <- left_join(IS.ND50s,tmp.IS.ND50s, by = 'sample')
  names(IS.ND50s)[names(IS.ND50s)=="ND50"] <- paste0("assay_run",run)
  
  ### visualize the fit:
  # make curve data for plot:
  dilution <- 10^seq(log10(min(fit.data.IS$corrected.dilution,na.rm=TRUE)),log10(max(fit.data.IS$corrected.dilution,na.rm = TRUE)),length.out=1e3)
  curve.fit.data.IS <- data.frame(dilution=dilution)
  for(i in 1:n.samples){
    # compute the curve y-values and save them as "tmp":
    tmp <- log_fun(x=dilution,par=tempfit$estimate[c(n.plates+1,n.plates+1+i)]) # relative number of spots
    tmp2 <- tempfit$estimate[ceiling(i/7)]*tmp # fit to assay
    
    # add it to the data:
    if(i==1){
      curve.fit.data.IS <- dplyr::mutate(curve.fit.data.IS,prob=tmp2,sample=rep(i,1e3),assay=rep(run,1e3))
    }else{
      curve.fit.data.IS <- rbind(curve.fit.data.IS,data.frame(dilution=dilution,prob=tmp2,sample=rep(i,1e3),assay=rep(run,1e3)))
    }
  }
  
  # figures for each assay (data and VOC):
  fig.IS.fit <- plot.IS.fit(fit.data = fit.data.IS[fit.data.IS$assay==run,], assay = run, curve.data = curve.fit.data.IS)
  fig.fit.VOC <- plot.IS.fit.VOC(fit.data = fit.data.IS[fit.data.IS$assay==run,], assay = run, fit.output = tempfit)
  
  # combine all figures:
  fig.fit <- fig.IS.fit + fig.fit.VOC + 
    plot_layout(design = "11111112")
  
  if(print_results){print(fig.fit)}
  if(save_results){
    ggsave(glue::glue("output/Figures/IS-assay{run}-fit_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
           width=plot_width*1.5,height=10*1.5,units="cm",plot=fig.fit)
  }
  
  # clean up:
  rm(tmp.fit.data.IS,tmp.IS.ND50s,dilution,curve.fit.data.IS,i,tmp,tmp2,fig.IS.fit,fig.fit.VOC,fig.fit)
}

# add GMT ND50 and CIs:
IS.ND50s <- IS.ND50s %>% 
  rowwise() %>% 
  mutate(ND50 = 10^(mean(log10(c(assay_run2,assay_run3)),na.rm=TRUE)),
         ND50_lower = 10^(mean(log10(c(assay_run2,assay_run3)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run2,assay_run3)))/sqrt(2)),
         ND50_upper = 10^(mean(log10(c(assay_run2,assay_run3)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run2,assay_run3)))/sqrt(2)))

# add standardized titers (IU/ml): 
IS.ND50s <- IS.ND50s %>% 
  mutate(assay_run2_IU = assay_run2/IS.ND50s$assay_run2[IS.ND50s$sample==stand_sample]*stand_IU,
         assay_run3_IU = assay_run3/IS.ND50s$assay_run3[IS.ND50s$sample==stand_sample]*stand_IU) %>% 
  rowwise() %>% 
  mutate(ND50_IU = 10^(mean(log10(c(assay_run2_IU,assay_run3_IU)),na.rm=TRUE)),
         ND50_IU_lower = 10^(mean(log10(c(assay_run2_IU,assay_run3_IU)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run2_IU,assay_run3_IU)))/sqrt(2)),
         ND50_IU_upper = 10^(mean(log10(c(assay_run2_IU,assay_run3_IU)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run2_IU,assay_run3_IU)))/sqrt(2))) %>%
  ungroup() %>% as.data.frame()

# save table with ND50s per assay run and for all assay runs:
if(save_results){
  # table with estimated ND50's for each sample:
  write_xlsx(IS.ND50s,glue::glue("output/Tables/IS-ND50s_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

### linear mixed-effects model to estimate variability between assay runs:
tmp.IS.ND50s <- pivot_longer(IS.ND50s[,names(IS.ND50s)%in%c("sample","assay_run2","assay_run3")], cols = starts_with("assay"),
                             names_to = "run", values_to = "ND50", names_pattern = "assay_run(.)")
IS.mod <- lmer(log10(ND50) ~ (1 | sample), data=tmp.IS.ND50s, REML = TRUE) # model without random effect for run
if(print_results){print(summary(IS.mod))}
run_SD <- sigma(IS.mod) # residual SD only
run_var_table <- data.frame(fold_diff_95perc = 10^((2)*run_SD),
                            mean_run_fold_diff = mean(summary(pmax(IS.ND50s$assay_run2/IS.ND50s$assay_run3,IS.ND50s$assay_run3/IS.ND50s$assay_run2))))

if(print_results){print(run_var_table)}
if(save_results){
  write_xlsx(run_var_table,glue::glue("output/Tables/IS-run-variability_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

# clean up:
rm(run,n.plates,n.samples,par.init,tempfit,tmp.IS.ND50s,IS.mod,
   run_SD,run_var_table)


# cleanup -----------------------------------------------------------------

rm(initial.slope.IS,initial.ND50.IS,show_fit,fit.data.IS)
