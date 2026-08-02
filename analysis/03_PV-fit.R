# -------------------------------------------------------------------------
#' 
#' Visualize data and fit to the Pseudovirus (PV) assay data
#' 
# -------------------------------------------------------------------------


# specify analysis and outputs --------------------------------------------

### estimating ND50s:
# initial parameter values:
initial.slope.PV <- log(3) # slope = 3, log-transformed

### figures:
show_fit <- TRUE # visualize data with the fit (TRUE) or without (FALSE)


# fit each run separately to estimate variability -------------------------

# least squares fit of logistic function to the data

### data:
fit.data <- PV.data[PV.data$sample!="CC",]
PV.ND50s <- CPE.original.ND50s[,names(CPE.original.ND50s)%in%c("sample","sample_name")]
PV.ND50s$sample <- as.numeric(PV.ND50s$sample)

### fitting:
for(run in unique(fit.data$assay)){
  tmp.fit.data <- fit.data[fit.data$assay==run,] # & !fit.data$sample%in%c("29","30")
  n.samples <- length(unique(tmp.fit.data$sample[!is.na(tmp.fit.data$dilution)]))
  
  # initial parameter values:
  par.init <- c(rep(mean(tmp.fit.data$outcome[tmp.fit.data$sample=="VC"]),length(unique(tmp.fit.data$plate))), # maximal fluorescence for each plate in each assay
                initial.slope.PV, # slope
                log(initial.PV.ND50s(tmp.fit.data))) # ND50 initial guess for each sample
  
  # minimize least squares:
  tempfit <- nlm(function(par){leastsq.PV(data=tmp.fit.data,par=par)},par.init,iterlim=1e3,hessian = TRUE)#,print.level = 2)
  
  # make a table with estimated ND50's for each sample:
  tmp.PV.ND50s <- data.frame(sample=as.numeric(unique(tmp.fit.data$sample[!tmp.fit.data$sample%in%c("VC","CC")])),
                             ND50=exp(tail(tempfit$estimate,n.samples)))
  
  # print & save parameter estimates:
  if(print_results){
    print(glue::glue("Least squares fit to PV data (assay run {run}):"))
    print(glue::glue("Sum of least squares: {tempfit$minimum}"))
    print(glue::glue("Exit flag: {tempfit$code}"))
    print(glue::glue("Iterations: {tempfit$iterations}"))
    print(glue::glue("AIC: {round(AIC_ls(tempfit$minimum,nrow(tmp.fit.data),length(tempfit$estimate)),3)}"))
    print(glue::glue("Number of estimated parameters: {length(tempfit$estimate)}"))
    print(tmp.PV.ND50s)
    print(glue::glue("Slope (95% CI): {round(exp(tempfit$estimate[6]),2)} ({round(exp(tempfit$estimate[6]-qnorm((1+alpha_CI)/2)*sqrt(diag(solve(tempfit$hessian[6,6])))),2)} - {round(exp(tempfit$estimate[6]+qnorm((1+alpha_CI)/2)*sqrt(diag(solve(tempfit$hessian[6,6])))),2)})"))
  }
  if(save_results){
    # least squares fitting result:
    save(tempfit,file = glue::glue("output/Fitting/PV-assay{run}-fit_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.RData"))
  }
  
  # add to ND50s table:
  PV.ND50s <- left_join(PV.ND50s,tmp.PV.ND50s, by = 'sample')
  names(PV.ND50s)[names(PV.ND50s)=="ND50"] <- paste0("assay_run",run)
  
  ### visualize the fit:
  # make curve data for plot:
  dilution <- 10^seq(log10(min(tmp.fit.data$dilution,na.rm=TRUE)),log10(max(tmp.fit.data$dilution,na.rm = TRUE)),length.out=1e3)
  curve.tmp.fit.data <- data.frame(dilution=dilution)
  for(i in 1:length(unique(tmp.fit.data$sample[!is.na(tmp.fit.data$dilution)]))){ # for each sample
    # compute the curve y-values and save them as "tmp":
    tmp <- log_fun(x=dilution,par=tempfit$estimate[c(length(unique(tmp.fit.data$plate))*length(unique(tmp.fit.data$assay))+1, # slope
                                                     length(unique(tmp.fit.data$plate))*length(unique(tmp.fit.data$assay))+1+i)]) # ND50 
    tmp2 <- tempfit$estimate[ceiling(i/6)]*tmp # fit to assay
    
    # add it to the data:
    if(i==1){
      curve.tmp.fit.data <- dplyr::mutate(curve.tmp.fit.data,prob=tmp2,sample=rep(i,1e3),assay=rep(run,1e3))
    }else{
      curve.tmp.fit.data <- rbind(curve.tmp.fit.data,data.frame(dilution=dilution,prob=tmp2,sample=rep(i,1e3),assay=rep(run,1e3)))
    }
  }
  
  # visualize data with fit:
  fig.PV <- plot.PV(data = PV.data[PV.data$assay==run,],assay = run,add.fit = TRUE,fit.data = curve.tmp.fit.data)
  fig.PV.contr <- plot.PV.contr(data = PV.data[PV.data$assay==run,],assay = run,add.fit = TRUE,fit.output = tempfit)
  
  # combine figures:
  fig.PV <- fig.PV + fig.PV.contr + 
    plot_layout(design = "1111112")
  
  if(print_results){print(fig.PV)}
  if(save_results){
    ggsave(glue::glue("output/Figures/PV-assay{run}-fit_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
           width=plot_width*1.5,height=10*1.5,units="cm",plot=fig.PV)
  }
  
  # clean up:
  rm(tmp.fit.data,n.samples,par.init,tempfit,tmp.PV.ND50s,dilution,curve.tmp.fit.data,
     i,tmp,tmp2,fig.PV,fig.PV.contr)
  
}

# add GMT ND50 and CIs:
PV.ND50s <- PV.ND50s %>% 
  rowwise() %>% 
  mutate(ND50 = 10^(mean(log10(c(assay_run1,assay_run2)),na.rm=TRUE)),
         ND50_lower = 10^(mean(log10(c(assay_run1,assay_run2)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)),
         ND50_upper = 10^(mean(log10(c(assay_run1,assay_run2)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)))

# add standardized titers (IU/ml): 
PV.ND50s <- PV.ND50s %>% 
  mutate(assay_run1_IU = assay_run1/PV.ND50s$assay_run1[PV.ND50s$sample==stand_sample]*stand_IU,
         assay_run2_IU = assay_run2/PV.ND50s$assay_run2[PV.ND50s$sample==stand_sample]*stand_IU) %>% 
  rowwise() %>% 
  mutate(ND50_IU = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)),na.rm=TRUE)),
         ND50_IU_lower = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2)),
         ND50_IU_upper = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2))) %>%
  ungroup() %>% as.data.frame()

# save table with ND50s per assay run and for all assay runs:
if(save_results){
  # table with estimated ND50's for each sample:
  write_xlsx(PV.ND50s,glue::glue("output/Tables/PV-ND50s_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

### linear mixed-effects model to estimate variability between assay runs:
tmp.PV.ND50s <- pivot_longer(PV.ND50s[c(1:28),names(PV.ND50s)%in%c("sample","assay_run1","assay_run2")], cols = starts_with("assay"),
                             names_to = "run", values_to = "ND50", names_pattern = "assay_run(.)")
PV.mod <- lmer(log10(ND50) ~ (1 | sample), data=tmp.PV.ND50s, REML = TRUE) # model without random effect for run
if(print_results){print(summary(PV.mod))}
run_SD <- sigma(PV.mod) # residual SD only
run_var_table <- data.frame(fold_diff = 10^((2)*run_SD), 
                            mean_run_fold_diff = mean(summary(pmax(PV.ND50s$assay_run1[c(1:28)]/PV.ND50s$assay_run2[c(1:28)],PV.ND50s$assay_run2[c(1:28)]/PV.ND50s$assay_run1[c(1:28)]))))

if(print_results){print(run_var_table)}
if(save_results){
  write_xlsx(run_var_table,glue::glue("output/Tables/PV-run-variability_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

# clean up:
rm(run,tmp.PV.ND50s,PV.mod,run_SD,run_var_table)


# cleanup -----------------------------------------------------------------

rm(initial.slope.PV,show_fit,fit.data)
