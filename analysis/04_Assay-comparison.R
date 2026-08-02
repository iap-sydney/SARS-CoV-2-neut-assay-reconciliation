# -------------------------------------------------------------------------
#' 
#' Comparison of different neutralization assays
#' data: fits to Immunospot (IS), Pseudovirus (PV), and various CPE assays
#' output: figures and tables comparing these assays
#' 
# -------------------------------------------------------------------------


# set-up for analysis -----------------------------------------------------

# data source: use latest fitting data with the specified date
# (run the codes to fit to the data first)

# options for figures:
show_1to1_relationship <- TRUE # show the 1:1 relationship between titres
show_linear_fit <- TRUE # show linear fit to the data: all best fits use standard linear regression
show_prediction <- TRUE # show the predicted relationship between titres

# assays to plot vs IS ND50s:
IS_comparison <- c("PV","CPE_original","CPE_rem_low","CPE_rem_high")

# list of comparisons between CPE assays: (first assay is shown on the x-axis, second on the y-axis)
CPE_comparisons <- list(c("CPE_original","CPE_rem_low"),c("CPE_original","CPE_rem_high"))


# load and combine different ND50 data ------------------------------------

# all ND50s:
CPE.original.ND50s$sample <- as.numeric(CPE.original.ND50s$sample)
CPE.rem.low.ND50s$sample <- as.numeric(CPE.rem.low.ND50s$sample)
CPE.rem.high.ND50s$sample <- as.numeric(CPE.rem.high.ND50s$sample)

if(standardize){
  ND50s <- IS.ND50s[c(1:28)[-stand_sample],names(IS.ND50s)%in%c("sample","sample_name","ND50_IU","ND50_IU_lower","ND50_IU_upper")] %>% 
    rename(IS_ND50 = ND50_IU, IS_ND50_lower = ND50_IU_lower, IS_ND50_upper = ND50_IU_upper) %>% 
    left_join(PV.ND50s[,names(PV.ND50s)%in%c("sample","sample_name","ND50_IU","ND50_IU_lower","ND50_IU_upper")] %>% rename(PV_ND50 = ND50_IU, PV_ND50_lower = ND50_IU_lower, PV_ND50_upper = ND50_IU_upper)) %>% 
    left_join(CPE.original.ND50s[,names(CPE.original.ND50s)%in%c("sample","sample_name","ND50_IU","ND50_IU_lower","ND50_IU_upper")] %>% rename(CPE_original_ND50 = ND50_IU, CPE_original_ND50_lower = ND50_IU_lower, CPE_original_ND50_upper = ND50_IU_upper)) %>%
    left_join(CPE.rem.low.ND50s[,names(CPE.rem.low.ND50s)%in%c("sample","sample_name","ND50_IU","ND50_IU_lower","ND50_IU_upper")] %>% rename(CPE_rem_low_ND50 = ND50_IU, CPE_rem_low_ND50_lower = ND50_IU_lower, CPE_rem_low_ND50_upper = ND50_IU_upper)) %>%
    left_join(CPE.rem.high.ND50s[,names(CPE.rem.high.ND50s)%in%c("sample","sample_name","ND50_IU","ND50_IU_lower","ND50_IU_upper")] %>% rename(CPE_rem_high_ND50 = ND50_IU, CPE_rem_high_ND50_lower = ND50_IU_lower, CPE_rem_high_ND50_upper = ND50_IU_upper))
}else{
  ND50s <- IS.ND50s[,names(IS.ND50s)%in%c("sample","sample_name","ND50","ND50_lower","ND50_upper")] %>% 
    rename(IS_ND50 = ND50, IS_ND50_lower = ND50_lower, IS_ND50_upper = ND50_upper) %>% 
    left_join(PV.ND50s[,names(PV.ND50s)%in%c("sample","sample_name","ND50","ND50_lower","ND50_upper")] %>% rename(PV_ND50 = ND50, PV_ND50_lower = ND50_lower, PV_ND50_upper = ND50_upper)) %>% 
    left_join(CPE.original.ND50s[,names(CPE.original.ND50s)%in%c("sample","sample_name","ND50","ND50_lower","ND50_upper")] %>% rename(CPE_original_ND50 = ND50, CPE_original_ND50_lower = ND50_lower, CPE_original_ND50_upper = ND50_upper)) %>%
    left_join(CPE.rem.low.ND50s[,names(CPE.rem.low.ND50s)%in%c("sample","sample_name","ND50","ND50_lower","ND50_upper")] %>% rename(CPE_rem_low_ND50 = ND50, CPE_rem_low_ND50_lower = ND50_lower, CPE_rem_low_ND50_upper = ND50_upper)) %>%
    left_join(CPE.rem.high.ND50s[,names(CPE.rem.high.ND50s)%in%c("sample","sample_name","ND50","ND50_lower","ND50_upper")] %>% rename(CPE_rem_high_ND50 = ND50, CPE_rem_high_ND50_lower = ND50_lower, CPE_rem_high_ND50_upper = ND50_upper))
}

# IS assay:
IS_2_fit <- load_fit_data(assay_name = "IS-assay2-fit")
IS_3_fit <- load_fit_data(assay_name = "IS-assay3-fit")
IS_fit_slope <- exp(mean(c(IS_2_fit$estimate[5],IS_3_fit$estimate[5])))
IS_fit_slope_var_log <- var(c(IS_2_fit$estimate[5],IS_3_fit$estimate[5]))/2

# CPE assays: original 
CPE_original_virions <- load_fit_data(assay_name = "CPE-fit-controls\\(data-original\\)")
CPE_original_vir_var <- as.numeric(solve(CPE_original_virions$hessian))
CPE_original_virions <- CPE_original_virions$estimate

# CPE assays: wash.low
CPE_rem_low_virions <- load_fit_data(assay_name = "CPE-fit-controls\\(data-rem.low\\)")
CPE_rem_low_vir_var <- as.numeric(solve(CPE_rem_low_virions$hessian))
CPE_rem_low_virions <- CPE_rem_low_virions$estimate

# CPE assays: wash.high 
CPE_rem_high_virions <- load_fit_data(assay_name = "CPE-fit-controls\\(data-rem.high\\)")
CPE_rem_high_vir_var <- as.numeric(solve(CPE_rem_high_virions$hessian))
CPE_rem_high_virions <- CPE_rem_high_virions$estimate

# save table with all ND50s and 95% CIs:
output.ND50s <- ND50s
output.ND50s$PV_ND50 <- paste(round(output.ND50s$PV_ND50,2)," (",round(output.ND50s$PV_ND50_lower,2)," to ",round(output.ND50s$PV_ND50_upper,2),")",sep="")
output.ND50s$IS_ND50 <- paste(round(output.ND50s$IS_ND50,2)," (",round(output.ND50s$IS_ND50_lower,2)," to ",round(output.ND50s$IS_ND50_upper,2),")",sep="")
output.ND50s$CPE_original_ND50 <- paste(round(output.ND50s$CPE_original_ND50,2)," (",round(output.ND50s$CPE_original_ND50_lower,2)," to ",round(output.ND50s$CPE_original_ND50_upper,2),")",sep="")
output.ND50s$CPE_rem_low_ND50 <- paste(round(output.ND50s$CPE_rem_low_ND50,2)," (",round(output.ND50s$CPE_rem_low_ND50_lower,2)," to ",round(output.ND50s$CPE_rem_low_ND50_upper,2),")",sep="")
output.ND50s$CPE_rem_high_ND50 <- paste(round(output.ND50s$CPE_rem_high_ND50,2)," (",round(output.ND50s$CPE_rem_high_ND50_lower,2)," to ",round(output.ND50s$CPE_rem_high_ND50_upper,2),")",sep="")
output.ND50s <- output.ND50s %>% dplyr::select(-contains("lower")) %>% dplyr::select(-contains("upper"))
output.ND50s[output.ND50s=="NA (NA to NA)"] <- NA

if(save_results){
  writexl::write_xlsx(output.ND50s,path = glue::glue("output/Tables/ND50s-all{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

# save a table with estimated number of virions in the different CPE assays:
n.vir.table <- data.frame(data = c("Original CPE (without inoculum removal)","CPE (with inoculum removal, lower virus inoculum)","CPE (with inoculum removal, higher virus inoculum)"),
                          n.virions.per.well = round(c(CPE_original_virions,CPE_rem_low_virions,CPE_rem_high_virions),2),
                          n.virions.CIs = c(paste0(round(CPE_original_virions-qnorm((1+alpha_CI)/2)*sqrt(CPE_original_vir_var),2)," to ",round(CPE_original_virions+qnorm((1+alpha_CI)/2)*sqrt(CPE_original_vir_var),2)),
                                            paste0(round(CPE_rem_low_virions-qnorm((1+alpha_CI)/2)*sqrt(CPE_rem_low_vir_var),2)," to ",round(CPE_rem_low_virions+qnorm((1+alpha_CI)/2)*sqrt(CPE_rem_low_vir_var),2)),
                                            paste0(round(CPE_rem_high_virions-qnorm((1+alpha_CI)/2)*sqrt(CPE_rem_high_vir_var),2)," to ",round(CPE_rem_high_virions+qnorm((1+alpha_CI)/2)*sqrt(CPE_rem_high_vir_var),2))))

if(save_results){
  writexl::write_xlsx(n.vir.table,path = glue::glue("output/Tables/CPE-number-of-virions-per-well_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

# clean up:
rm(IS_2_fit,IS_3_fit,output.ND50s,n.vir.table)


# CPE vs IS vs PV -----------------------------------------

# data for visualization:
summary.ND50s <- ND50s[,names(ND50s)%in%c("sample","sample_name","CPE_original_ND50","IS_ND50","PV_ND50")]
summary.ND50s <- summary.ND50s %>% tidyr::pivot_longer(cols=c("CPE_original_ND50","IS_ND50","PV_ND50"))
# add censoring information: censored if the estimated ND50 is below the lowest dilution in the assay
summary.ND50s <- dplyr::mutate(summary.ND50s, censored=rep(0,nrow(summary.ND50s)),
                               value_log10 = log10(value))
# summary.ND50s$censored[] <- 1 # no samples' ND50s are below the lowest dilution of that assay and sample

# paired t-test comparisons:
my_comparisons <- ggpubr::compare_means(value_log10~name, data=summary.ND50s, method = "t.test", paired = TRUE)
my_comparisons <- my_comparisons %>% mutate(print.p = paste("p =",round(p,4)))
my_comparisons$print.p[my_comparisons$print.p=="p = 0"] <- "p < 0.0001"

# specify data to visualize:
show_names <- c("CPE_original_ND50","IS_ND50","PV_ND50")
show_labels <- c("CPE","IS","PV")
if(length(show_names)==3){
  if(standardize){
    y_pos <- c(4.1,4.1*10^0.02,4.1*10^0.04)
  }else{
    y_pos <- c(4.2,4.2*1.06,4.2*1.06^2)
  }
}else{
  if(standardize){
    y_pos <- 2
  }else{
    y_pos <- 4.2
  }
}

# figure:
fig <- ggplot(summary.ND50s[summary.ND50s$name%in%show_names,],
              aes(x=as.factor(name),y=value,group=sample,color=factor(sample),shape=as.factor(censored))) + 
  # add limit of detection: # the lowest dilution varies between assays and samples, so don't add a limit of dilution to the figure
  # geom_hline(yintercept = 20, linetype=1, color="gray") +
  # geom_text(aes(x=2,y=25,label="limit of detection"), hjust=0, color="gray") +
  # add data:
  geom_point() +
  geom_line() +
  ggpubr::geom_signif(data = my_comparisons[my_comparisons$group1%in%show_names & my_comparisons$group2%in%show_names,],
                      inherit.aes = FALSE, manual = TRUE, vjust = -0.5, 
                      aes(xmin=group1,xmax=group2,annotations = print.p,y_position=y_pos)) + 
  # legend, axis, theme, etc.:
  scale_shape_manual(values=c(19,1)) +
  {if(standardize) scale_y_log10(breaks=c(2.5e2,5e2,7.5e2,1e3,2.5e3,5e3,7.5e3,1e4), minor_breaks = c(5e2,1e3,5e3,1e4))} + 
  {if(!standardize) scale_y_log10(breaks=2^(c(3:14)), minor_breaks = 2^(c(3:14)))} +
  scale_x_discrete(breaks=show_names,labels=show_labels, expand = c(1e-1,1e-1)) +
  # coord_cartesian(expand = 0.1) +
  labs(x=NULL,y=glue::glue("ND50{ifelse(standardize,' (IU/ml)','')}"),title=glue::glue("Comparison of ND50s across different assays")) +
  theme_bw() + 
  theme(legend.position="none")

if(print_results){print(fig)}
if(save_results){
  ggsave(glue::glue("output/Figures/ND50s-{paste(show_labels,collapse='-')}{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
         width=plot_width,height=plot_width,units="cm",plot=fig)
}

# clean up:
rm(summary.ND50s,my_comparisons,show_names,show_labels,y_pos,fig)


# correlations between ND50s ---------------------------------------------- 

tmp_data <- ND50s[,!(names(ND50s)%in%c("sample","sample_name") | grepl("_lower",names(ND50s)) | grepl("_upper",names(ND50s)))]
tmp_cor <- Hmisc::rcorr(as.matrix(tmp_data), type = "spearman")
cor_table <- as.data.frame(tmp_cor$r)
cor_pval_table <- as.data.frame(tmp_cor$P)

if(print_results){print(cor_table)}
if(save_results){
  writexl::write_xlsx(cor_table,path = glue::glue("output/Tables/Spearman-correlations{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
  writexl::write_xlsx(cor_pval_table,path = glue::glue("output/Tables/Spearman-cor-pvals{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

# clean up:
rm(tmp_data,tmp_cor,cor_table,cor_pval_table)


# fold differences between titers with different assays ------------------- 

### geometric mean of fold-differences between assays: 
tmp_data <- ND50s[,!(names(ND50s)%in%c("sample","sample_name") | grepl("_lower",names(ND50s)) | grepl("_upper",names(ND50s)))]
gmt_fold_diff <- gmt_fold_diff_CIs <- matrix(NA,ncol(tmp_data),ncol(tmp_data))
for(i in 1:ncol(tmp_data)){
  for(j in 1:ncol(tmp_data)){
    gmt_fold_diff[i,j] <- exp(mean(log(unlist(tmp_data[,i])/unlist(tmp_data[,j])), na.rm = T))
    gmt_fold_diff_CIs[i,j] <- paste(round(exp(mean(log(unlist(tmp_data[,i])/unlist(tmp_data[,j])),na.rm = T)-qnorm((1+alpha_CI)/2)*sd(log(unlist(tmp_data[,i])/unlist(tmp_data[,j])),na.rm = T)/nrow(tmp_data)),2),
                                    round(exp(mean(log(unlist(tmp_data[,i])/unlist(tmp_data[,j])),na.rm = T)+qnorm((1+alpha_CI)/2)*sd(log(unlist(tmp_data[,i])/unlist(tmp_data[,j])),na.rm = T)/nrow(tmp_data)),2),sep = " - ")
  }
}
gmt_fold_diff <- as.data.frame(gmt_fold_diff)
gmt_fold_diff_CIs <- as.data.frame(gmt_fold_diff_CIs)
names(gmt_fold_diff) <- rownames(gmt_fold_diff) <- names(gmt_fold_diff_CIs) <- rownames(gmt_fold_diff_CIs) <- names(tmp_data)

if(print_results){
  print(gmt_fold_diff)
  print(gmt_fold_diff_CIs)
}
if(save_results){
  writexl::write_xlsx(gmt_fold_diff,path = glue::glue("output/Tables/Geomean-titer-fold-diff{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
  writexl::write_xlsx(gmt_fold_diff_CIs,path = glue::glue("output/Tables/Geomean-titer-fold-diff-CIs{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.xlsx"))
}

# clean up:
rm(tmp_data,gmt_fold_diff,gmt_fold_diff_CIs,i,j)


# Comparison of different assays' ND50s with IS ND50s ------------- 

fig_list <- vector(mode = "list", length=length(IS_comparison))
for(assay in IS_comparison){
  # make data for prediction of ND50s:
  if(show_prediction & assay%in%c("CPE_original","CPE_rem_low","CPE_rem_high") & !standardize){
    if(assay=="CPE_original"){
      prediction <- predict_ND50s(n_virions = CPE_original_virions, n_virions_var = CPE_original_vir_var)
    }else if(assay=="CPE_rem_low"){
      prediction <- predict_ND50s(n_virions = CPE_rem_low_virions, n_virions_var = CPE_rem_low_vir_var)
    }else if(assay=="CPE_rem_high"){
      prediction <- predict_ND50s(n_virions = CPE_rem_high_virions, n_virions_var = CPE_rem_high_vir_var)
    }
  }else{
    prediction <- NULL
  }
  
  fig_list[[which(IS_comparison==assay)]] <- plot_ND50_comp(assay_y = assay, prediction = prediction)
  
  if(print_results){print(fig_list[[which(IS_comparison==assay)]])}
  if(save_results){
    ggsave(glue::glue("output/Figures/{assay}-vs-IS{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
           width=plot_width,height=plot_width,units="cm",plot=fig_list[[which(IS_comparison==assay)]])
  }
}

# # combine figures to make a figure for the manuscript:
# IS_comparison <- c("CPE_original","CPE_rem_low","CPE_rem_high")
# # run code above to generate panels
# for(i in 1:length(fig_list)){
#   fig_list[[i]] <- fig_list[[i]] + ggtitle(NULL) + # remove titles
#     coord_cartesian(xlim=c(128,16384),ylim=c(8,8192)) + # same axes scales for all panels
#     theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + # changes x-axis tick labels
#     labs(x="IS ND50")
# }
# fig_combined <- fig_list[[1]] + fig_list[[3]] + fig_list[[2]] + plot_layout(ncol = 2) +
#   plot_annotation(title = "Comparison of CPE and IS Assays", tag_levels = "A")
# if(print_results){print(fig_combined)}
# if(save_results){
#   ggsave(glue::glue("output/Figures/Fig-3_CPE-assays-vs_IS_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
#          width=plot_width*1.25,height=plot_width*1.25,units="cm",plot=fig_combined)
# }
# rm(IS_comparison,i,fig_combined)
# 
# # clean up:
# rm(fig_list,assay,prediction)


# comparisons between different CPE assays --------------------------------

for(i in c(1:length(CPE_comparisons))){
  # no predictions for these plots
  
  fig <- plot_ND50_comp(assay_x = CPE_comparisons[[i]][1], assay_y = CPE_comparisons[[i]][2])
  
  if(print_results){print(fig)}
  if(save_results){
    ggsave(glue::glue("output/Figures/{CPE_comparisons[[i]][2]}-vs-{CPE_comparisons[[i]][1]}{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
           width=plot_width,height=plot_width,units="cm",plot=fig)
  }
}

# clean up:
rm(i,fig)


# combine different panels to one figure ----------------------------------

# run the corresponding code above to generate panels for the figure figures to be combined

### Figure 1:
# # panel A: CPE vs IS vs PV (run code above)
# fig1 <- fig + ggtitle(NULL) # save figure and remove title
# # panel B: PV vs IS
# fig2 <- plot_ND50_comp(assay_x = "IS", assay_y = "PV") + ggtitle(NULL)
# # panel C: CPE (original) vs IS without prediction
# fig3 <- plot_ND50_comp(assay_x = "IS", assay_y = "CPE_original") + ggtitle(NULL)
# # combine panels:
# fig_combined <- fig1 + fig2 + fig3 + plot_layout(ncol = 2) +
#   plot_annotation(title = "ND50 Comparison Across Assays", tag_levels = "A")
# if(print_results){print(fig_combined)}
# if(save_results){
#   ggsave(glue::glue("output/Figures/Fig-1_ND50-comparison{ifelse(standardize,'-standardized','')}_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
#          width=plot_width*1.25,height=plot_width*1.25,units="cm",plot=fig_combined)
# }
# # clean up:
# rm(fig1,fig2,fig3,fig_combined)


# cleanup -----------------------------------------------------------------

rm(show_1to1_relationship,show_linear_fit,show_prediction,IS_comparison,CPE_comparisons,
   ND50s,IS_fit_slope,IS_fit_slope_var_log,CPE_original_virions,CPE_original_vir_var,
   CPE_rem_low_virions,CPE_rem_low_vir_var,CPE_rem_high_virions,CPE_rem_high_vir_var)
