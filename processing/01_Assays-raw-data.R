# -------------------------------------------------------------------------
#' 
#' Loads and processes the raw data for the CPE assay (original), Immunospot 
#' (IS) assay, and Pseudovirus (PV) assay
#' 
# -------------------------------------------------------------------------


# parameters for data processing ------------------------------------------

CPE.correction <- 2 # add an extra 2-fold dilution to the CPE assay dilutions
# Adding the virus dilution to the serum leads to an additional 2-fold dilution to be included here.
# This is already included in the dilutions reported in the IS and PV assays and only needs to be accounted for in the CPE assay.


# load all raw data -------------------------------------------------------

# summary data:
summary.data <- read_excel("raw-data/FM-24-003 & 004 - Summary - Neutralisation assay comparison.xlsx",sheet="Summary")

# CPE assays raw back titration data (original CPE assay):
raw.CPE.1.contr <- read_excel("raw-data/FM-24-004 - CPE MN (WHO IE Discover only; Vic01) results.xlsx", sheet = "MN results 1 (2x)")
raw.CPE.2.contr <- read_excel("raw-data/FM-24-004 - CPE MN (WHO IE Discover only; Vic01) results.xlsx", sheet = "MN results 2 (2x)")

# IS assays raw data:
raw.IS.2.data <- read_excel("raw-data/FM-24-003 & 004 - Summary - Neutralisation assay comparison.xlsx",sheet="ImmunoSpot assay 2 raw data")
raw.IS.3.data <- read_excel("raw-data/FM-24-003 & 004 - Summary - Neutralisation assay comparison.xlsx",sheet="ImmunoSpot assay 3 raw data")

# PV assay raw data: 
raw.PV.1.data <- read_excel("raw-data/FM-24-003 & 004 - Summary - Neutralisation assay comparison.xlsx",sheet="Pseudoneut assay 1 raw data")
raw.PV.2.data <- read_excel("raw-data/FM-24-003 & 004 - Summary - Neutralisation assay comparison.xlsx",sheet="Pseudoneut assay 2 raw data")


# process IS assay data -------------------------------------------

# re-format IS data:
IS.2.data <- data_processing(raw.data = raw.IS.2.data,experiment="IS")
IS.3.data <- data_processing(raw.data = raw.IS.3.data,experiment="IS")

### combine IS assays and add variables:
IS.data <- rbind(dplyr::mutate(IS.2.data,assay=rep(2,nrow(IS.2.data))),
                 dplyr::mutate(IS.3.data,assay=rep(3,nrow(IS.3.data))))

# add starting dilution to the IS assay data:
starting.dilution.data <- summary.data[c(3:32),c(1,8,11)]
names(starting.dilution.data) <- c("sample","assay2","assay3")

corrected.dilution <- sapply(c(1:nrow(IS.data)),function(x){ifelse(any(starting.dilution.data$sample==as.numeric(IS.data$sample[x])),
                                                                   IS.data$dilution[x]*as.numeric(starting.dilution.data[starting.dilution.data$sample==as.numeric(IS.data$sample[x]),ifelse(IS.data$assay[x]==2,2,3)])/20,
                                                                   NA)})
IS.data <- dplyr::mutate(IS.data, corrected.dilution)

# add grouping variable to IS data:
IS.data <- dplyr::mutate(IS.data, assay.plate.sample.rep = paste(IS.data$assay,IS.data$plate,IS.data$sample,IS.data$rep,sep="."))


# process CPE assay data --------------------------------------------------

# CPE ND50s calculated using Reed-Muench method:
CPE.original.ND50s <- summary.data[c(3:32),c(1,2,5,6)]
names(CPE.original.ND50s) <- c("sample","sample_name","assay_run1","assay_run2")
CPE.original.ND50s <- CPE.original.ND50s %>% mutate(across(c(assay_run1,assay_run2), as.numeric))

# replace ND50s below the starting dilution with NA:
CPE.original.ND50s$assay_run1[CPE.original.ND50s$assay_run1 <= 20] <- NA
CPE.original.ND50s$assay_run2[CPE.original.ND50s$assay_run2 <= 20] <- NA

# add CIs:
CPE.original.ND50s <- CPE.original.ND50s %>% 
   rowwise() %>% 
   mutate(ND50 = 10^(mean(log10(c(assay_run1,assay_run2)),na.rm=TRUE)),
          ND50_lower = 10^(mean(log10(c(assay_run1,assay_run2)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)),
          ND50_upper = 10^(mean(log10(c(assay_run1,assay_run2)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2))) %>% 
   ungroup()

# CPE dilution correction:
CPE.original.ND50s[,names(CPE.original.ND50s)%in%c("assay_run1","assay_run2","ND50","ND50_lower","ND50_upper")] <- 
   CPE.original.ND50s[,names(CPE.original.ND50s)%in%c("assay_run1","assay_run2","ND50","ND50_lower","ND50_upper")]*CPE.correction

# add standardized titers (IU/ml): 
CPE.original.ND50s <- CPE.original.ND50s %>% 
   mutate(assay_run1_IU = assay_run1/CPE.original.ND50s$assay_run1[CPE.original.ND50s$sample==stand_sample]*stand_IU,
          assay_run2_IU = assay_run2/CPE.original.ND50s$assay_run2[CPE.original.ND50s$sample==stand_sample]*stand_IU) %>% 
   rowwise() %>% 
   mutate(ND50_IU = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)),na.rm=TRUE)),
          ND50_IU_lower = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2)),
          ND50_IU_upper = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2))) %>%
   ungroup() %>% as.data.frame()


# process PV data ---------------------------------------------------------

### PV assay 1 data:
# make data frame:
PV.1.data <- data.frame(plate=rep(c(1:5),rep(96,5)),
                        row=rep(c("A","B","C","D","E","F","G","H"),12*5),
                        column=rep(rep(c(1:12),rep(8,12)),5),
                        sample=rep(c(1:30),rep(16,30)), # add controls below
                        rep=rep(rep(c(1,2),c(8,8)),6*5),
                        dilution=rep(c(20,40,80,160,320,640,1280,NA),12*5),
                        outcome=rep(NA,5*96))

# add controls:
PV.1.data$sample[PV.1.data$row=="H" & PV.1.data$column%in%c(1:6)] <- "CC"
PV.1.data$sample[PV.1.data$row=="H" & PV.1.data$column%in%c(7:12)] <- "VC"
PV.1.data$sample <- factor(PV.1.data$sample,levels = c(1:30,"CC","VC"))
PV.1.data$rep[is.na(PV.1.data$dilution)] <- NA

# adjust for starting dilution: starting dilution is 1:200 for all samples but sample 30 with starting dilution 1:20
PV.1.data$dilution[PV.1.data$sample%in%c(1:29)] <- PV.1.data$dilution[PV.1.data$sample%in%c(1:29)]*10 

# fill the outcome variable:
for(i in 1:nrow(PV.1.data)){
   PV.1.data$outcome[i] <- raw.PV.1.data[which(raw.PV.1.data[c(2:nrow(raw.PV.1.data)),5]==PV.1.data$column[i])[PV.1.data$plate[i]]+1, # select correct column and plate
                                         which(raw.PV.1.data[13,]==PV.1.data$row[i])] # select correct row
}

### PV assay 2 data:
# make data frame:
PV.2.data <- data.frame(plate=rep(c(1:5),rep(96,5)),
                        row=rep(c("A","B","C","D","E","F","G","H"),12*5),
                        column=rep(rep(c(1:12),rep(8,12)),5),
                        sample=rep(c(1:30),rep(16,30)), # add controls below
                        rep=rep(rep(c(1,2),c(8,8)),6*5),
                        dilution=rep(c(20,40,80,160,320,640,1280,NA),12*5),
                        outcome=rep(NA,5*96))

# add controls:
PV.2.data$sample[PV.2.data$row=="H" & PV.2.data$column%in%c(1:6)] <- "CC"
PV.2.data$sample[PV.2.data$row=="H" & PV.2.data$column%in%c(7:12)] <- "VC"
PV.2.data$sample <- factor(PV.2.data$sample,levels = c(1:30,"CC","VC"))
PV.2.data$rep[is.na(PV.2.data$dilution)] <- NA

# adjust for starting dilution: starting dilution is 1:200 for all samples but samples 8, 11, 14, 29, 30 (starting dilution 1:20) and samples 9, 15, 18, 23, 26 (starting dilution 1:1000)
PV.2.data$dilution[PV.2.data$sample%in%c(1:7,10,12,13,16,17,19:22,24,25,27,28)] <- PV.2.data$dilution[PV.2.data$sample%in%c(1:7,10,12,13,16,17,19:22,24,25,27,28)]*10 
PV.2.data$dilution[PV.2.data$sample%in%c(9,15,18,23,26)] <- PV.2.data$dilution[PV.2.data$sample%in%c(9,15,18,23,26)]*50 

# fill the outcome variable:
for(i in 1:nrow(PV.2.data)){
   PV.2.data$outcome[i] <- raw.PV.2.data[which(raw.PV.2.data[c(2:nrow(raw.PV.2.data)),5]==PV.2.data$column[i])[PV.2.data$plate[i]]+1, # select correct column and plate
                                         which(raw.PV.2.data[13,]==PV.2.data$row[i])] # select correct row
}

# combine data for assays 1 & 2:
PV.data <- rbind(dplyr::mutate(PV.1.data,assay=rep(1,nrow(PV.1.data))),
                 dplyr::mutate(PV.2.data,assay=rep(2,nrow(PV.2.data))))

# make outcome numeric:
PV.data$outcome <- as.numeric(PV.data$outcome)


# process CPE assay back titration data -----------------------------------

# make control data: 
CPE.1.vc <- data.frame(dilution = c(8*10^c(-1:-6),1), # Back titration: 10-fold dilutions of virus stock, CPE assay virus controls: 8-fold dilution of virus stock # dilutions relative to virus controls
                       n.pos=c(str_count(as.character(raw.CPE.1.contr[c(39:40),4:9]), pattern ="\\+"), # back titration data
                               str_count(as.character(raw.CPE.1.contr[c(7:36),15]), pattern ="\\+") + # virus control
                                  sum(str_count(as.character(raw.CPE.1.contr[c(35:36),c(5:14)]), pattern ="\\+"))), # negative control samples (29 and 30)
                       n.neg=c(str_count(as.character(raw.CPE.1.contr[c(39:40),4:9]), pattern ="-"), # back titration data
                               str_count(as.character(raw.CPE.1.contr[c(7:36),15]), pattern ="-") + # virus control
                                  sum(str_count(as.character(raw.CPE.1.contr[c(35:36),c(5:14)]), pattern ="-")))) # negative control samples (29 and 30)
CPE.1.vc <- dplyr::mutate(CPE.1.vc,n.rep=n.pos+n.neg)

CPE.2.vc <- data.frame(dilution = c(8*10^c(-1:-6),1), # Back titration: 10-fold dilutions of virus stock, CPE assay virus controls: 8-fold dilution of virus stock # dilutions relative to virus controls
                       n.pos=c(str_count(as.character(raw.CPE.2.contr[c(39:40),4:9]), pattern ="\\+"), # back titration data
                               str_count(as.character(raw.CPE.2.contr[c(7:36),15]), pattern ="\\+") + # virus control
                                  sum(str_count(as.character(raw.CPE.2.contr[c(35:36),c(5:14)]), pattern ="\\+"))), # negative control samples (29 and 30)
                       n.neg=c(str_count(as.character(raw.CPE.2.contr[c(39:40),4:9]), pattern ="-"), # back titration data
                               str_count(as.character(raw.CPE.2.contr[c(7:36),15]), pattern ="-") + # virus control
                                  sum(str_count(as.character(raw.CPE.2.contr[c(35:36),c(5:14)]), pattern ="-")))) # negative control samples (29 and 30)
CPE.2.vc <- dplyr::mutate(CPE.2.vc,n.rep=n.pos+n.neg)

# combine both into one data frame:
CPE.original.vc <- rbind(dplyr::mutate(CPE.1.vc,assay=rep(1,nrow(CPE.1.vc))),
                         dplyr::mutate(CPE.2.vc,assay=rep(2,nrow(CPE.2.vc))))


# sample names ------------------------------------------------------------

sample_names <- CPE.original.ND50s[,names(CPE.original.ND50s)%in%c("sample","sample_name")]


# control data for virus inoculum -----------------------------------------

# load data:
control.data <- read_excel("raw-data/Vic01 and Wuhan-1 PV titration results used for assays.xlsx")
IS.control.data <- control.data[control.data$virus=="Vic01",names(control.data)!="RLU (pseudovirus)"]
PV.control.data <- control.data[control.data$virus=="Wuhan-1 pseudovirus",names(control.data)!="Spot counts (Immunospot)"]

### IS control data
IS.VOC.spots <- data.frame(min = min(as.numeric(unlist(IS.data$spots[IS.data$sample=="VOC"]))), max = max(as.numeric(unlist(IS.data$spots[IS.data$sample=="VOC"]))))
IS.dilution <- 300 # dilution used in the final assay runs
# linear fit to the data:
lm_mod <- lm(log10(`Spot counts (Immunospot)`) ~ log10(dilution), data=IS.control.data[IS.control.data$dilution>=243 & IS.control.data$dilution<=6561,]) # range selected based on adjusted R^2
lm_mod_confb <- predict(lm_mod, data.frame(dilution=10^seq(log10(243),log10(6561),length.out=1e3)), se.fit = TRUE, interval = "confidence", level = alpha_CI)
lm_mod_region <- data.frame(dilution=10^seq(log10(243),log10(6561),length.out=1e3), estim = 10^lm_mod_confb$fit[,1], lower=10^lm_mod_confb$fit[,2], upper=10^lm_mod_confb$fit[,3])

# figure:
fig.IS <- ggplot(IS.control.data[IS.control.data$`Spot counts (Immunospot)`>0,],aes(x=dilution,y=`Spot counts (Immunospot)`)) +
   # show linear fit:
   geom_ribbon(data=lm_mod_region,inherit.aes=FALSE,aes(x=dilution,ymin=lower, ymax=upper),fill="black", alpha = 0.1) +
   geom_abline(slope = lm_mod$coef[2],intercept = lm_mod$coef[1], color="black") +
   # show range of data on the VOC plates:
   annotate("rect",xmin=1, xmax=1e6, ymin=IS.VOC.spots$min, ymax=IS.VOC.spots$max, fill='orangered', alpha = 0.3) +
   annotate("text",x=3, y=exp(mean(log(c(IS.VOC.spots$min,IS.VOC.spots$max)))), label = "Range of VOC spots", color = "orangered", hjust = 0) +
   # show dilution in the final assay run:
   geom_vline(xintercept = IS.dilution, color="orangered", alpha=0.8) +
   annotate("text",x=IS.dilution, y=350, label = "Dilution used in the IS assay", color = "orangered", hjust = 0) +
   # show data: 
   geom_point() + 
   # linear fit annotation:
   annotate("text",x=3,y=50,label=paste("Best linear fit:"), hjust=0) +
   annotate("text",x=3,y=50*0.8,label=paste("adjusted R-squared: ",round(summary(lm_mod)$adj.r.squared,3),sep=""), hjust=0) +
   annotate("text",x=3,y=50*0.8^2,label=paste("Intercept: ",round(lm_mod$coef[1],3)," (95% CI: ",round(confint(lm_mod)[1,1],3)," to ",round(confint(lm_mod)[1,2],3),")",sep=""), hjust=0) +
   annotate("text",x=3,y=50*0.8^3,label=paste("Slope: ",round(lm_mod$coef[2],3)," (95% CI: ",round(confint(lm_mod)[2,1],3)," to ",round(confint(lm_mod)[2,2],3),")",sep=""), hjust=0) +
   # legend, axis, theme, etc.:
   scale_x_log10(breaks=unique(IS.control.data$dilution)) + 
   scale_y_log10() + 
   coord_cartesian(xlim=c(3,177147),ylim=c(1,500)) +
   labs(x="Dilution of virus inoculum", y="Spot count", 
        title="IS assay: spot count vs dilution of the virus inoculum") +
   theme_bw() + 
   theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

if(print_results){print(fig.IS)}
if(save_results){
   ggsave(glue::glue("output/Figures/IS-spots-vs-dilution_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
          width=plot_width,height=plot_width,units="cm",plot=fig.IS)
}

### PV control data
PV.control.RLU <- data.frame(min = min(PV.data$outcome[PV.data$sample=="VC"]), max = max(PV.data$outcome[PV.data$sample=="VC"]))
PV.dilution <- 12.9 # dilution in the final assay runs
# linear fit to the data:
lm_mod <- lm(log10(`RLU (pseudovirus)`) ~ log10(dilution), data=PV.control.data[PV.control.data$dilution<=1e4,]) 
lm_mod_confb <- predict(lm_mod, data.frame(dilution=10^seq(1,4,length.out=1e3)), se.fit = TRUE, interval = "confidence", level = alpha_CI)
lm_mod_region <- data.frame(dilution=10^seq(1,4,length.out=1e3), estim = 10^lm_mod_confb$fit[,1], lower=10^lm_mod_confb$fit[,2], upper=10^lm_mod_confb$fit[,3])

# figure:
fig.PV <- ggplot(PV.control.data,aes(x=dilution,y=`RLU (pseudovirus)`)) +
   # show linear fit:
   geom_ribbon(data=lm_mod_region,inherit.aes=FALSE,aes(x=dilution,ymin=lower, ymax=upper),fill="black", alpha = 0.1) +
   geom_abline(slope = lm_mod$coef[2],intercept = lm_mod$coef[1], color="black") +
   # show range of data on the VOC plates:
   annotate("rect",xmin=1e-1, xmax=1e14, ymin=PV.control.RLU$min, ymax=PV.control.RLU$max, fill='orangered', alpha = 0.3) +
   annotate("text",x=1e4, y=exp(mean(log(c(PV.control.RLU$min,PV.control.RLU$max)))), label = "Range of RLU for controls", color = "orangered", hjust = 0) +
   # show dilution in the final assay run:
   geom_vline(xintercept = PV.dilution, color="orangered", alpha=0.8) +
   annotate("text",x=PV.dilution, y=350, label = "Dilution used in the PV assay", color = "orangered", hjust = 0) +
   # show data: 
   geom_point() + 
   # linear fit annotation:
   annotate("text",x=1e4,y=10^4.9,label=paste("Best linear fit:"), hjust=0) + 
   annotate("text",x=1e4,y=10^4.7,label=paste("adjusted R-squared: ",round(summary(lm_mod)$adj.r.squared,3),sep=""), hjust=0) + 
   annotate("text",x=1e4,y=10^4.5,label=paste("Intercept: ",round(lm_mod$coef[1],3)," (95% CI: ",round(confint(lm_mod)[1,1],3)," to ",round(confint(lm_mod)[1,2],3),")",sep=""), hjust=0) + 
   annotate("text",x=1e4,y=10^4.3,label=paste("Slope: ",round(lm_mod$coef[2],3)," (95% CI: ",round(confint(lm_mod)[2,1],3)," to ",round(confint(lm_mod)[2,2],3),")",sep=""), hjust=0) + 
   # legend, axis, theme, etc.:
   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))) + 
   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                 labels = trans_format("log10", math_format(10^.x))) +
   coord_cartesian(xlim=c(1e1,1e12),ylim=c(1e2,1.1e6)) +
   labs(x="Dilution of virus inoculum", y="RLU (Relative Luminescence Units)", 
        title="PV assay: RLU vs dilution of the virus inoculum") +
   theme_bw() + 
   theme(legend.position="none",axis.text.x = element_text(angle = 0, vjust = 0, hjust = 0.5))

if(print_results){print(fig.PV)}
if(save_results){
   ggsave(glue::glue("output/Figures/PV-RLU-vs-dilution_{format(Sys.time(),'%Y-%m-%d_%H%M%S')}.pdf"),
          width=plot_width,height=plot_width,units="cm",plot=fig.PV)
}


### clean up
rm(control.data,IS.control.data,PV.control.data,IS.VOC.spots,IS.dilution,lm_mod,
   lm_mod_confb,lm_mod_region,fig.IS,PV.control.RLU,PV.dilution,fig.PV)


# cleanup -----------------------------------------------------------------

rm(summary.data,i,raw.IS.2.data,raw.IS.3.data,starting.dilution.data,corrected.dilution,
   IS.2.data,IS.3.data,CPE.correction,raw.PV.1.data,raw.PV.2.data,PV.1.data,PV.2.data,
   raw.CPE.1.contr,raw.CPE.2.contr,CPE.1.vc,CPE.2.vc)
