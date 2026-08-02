# -------------------------------------------------------------------------
#' 
#' Loads and processes the raw data of the CPE assay with removing the inoculum
#' and low or high virus inoculum
#' 
# -------------------------------------------------------------------------


# parameters for data processing ------------------------------------------

CPE.correction <- 2 # add an extra 2-fold dilution to the CPE assay dilutions
# Adding the virus dilution to the serum leads to an additional 2-fold dilution to be included here.
# This is already included in the dilutions reported in the IS and PV assays and only needs to be accounted for in the CPE assay.


# -------------------------------------------------------------------------
#' 
#' LOWER VIRAL INOCULUM
#' 
# -------------------------------------------------------------------------


# load raw data -----------------------------------------------------------

# CPE data from the assay with removing the inoculum and lower virus inoculum:
raw.CPE.rem.low.1.data <- read_excel("raw-data/FM-24-030 - CPE MN (WHO IE Discover only; Vic01; wash after adsorption)_2024-07-23.xlsx", sheet = "MN results 1")
raw.CPE.rem.low.2.data <- read_excel("raw-data/FM-24-030 - CPE MN (WHO IE Discover only; Vic01; wash after adsorption)_2024-07-23.xlsx", sheet = "MN results 2")


# process CPE assay data --------------------------------------------------

# CPE ND50s calculated using Reed-Muench method:
CPE.rem.low.ND50s <- left_join(raw.CPE.rem.low.1.data[c(7:36),c(2,3,33)],raw.CPE.rem.low.2.data[c(7:36),c(2,3,33)],by=join_by(`...2`,`...3`))
names(CPE.rem.low.ND50s) <- c("sample","sample_name","assay_run1","assay_run2")
CPE.rem.low.ND50s <- CPE.rem.low.ND50s %>% mutate(across(c(assay_run1,assay_run2), as.numeric))

# replace ND50s below the starting dilution with NA:
CPE.rem.low.ND50s$assay_run1[CPE.rem.low.ND50s$assay_run1 <= 20] <- NA
CPE.rem.low.ND50s$assay_run2[CPE.rem.low.ND50s$assay_run2 <= 20] <- NA

# add CIs:
CPE.rem.low.ND50s <- CPE.rem.low.ND50s %>% 
  rowwise() %>% 
  mutate(ND50 = 10^(mean(log10(c(assay_run1,assay_run2)),na.rm=TRUE)),
         ND50_lower = 10^(mean(log10(c(assay_run1,assay_run2)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)),
         ND50_upper = 10^(mean(log10(c(assay_run1,assay_run2)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)))

# CPE dilution correction:
CPE.rem.low.ND50s[,names(CPE.rem.low.ND50s)%in%c("assay_run1","assay_run2","ND50","ND50_lower","ND50_upper")] <- 
  CPE.rem.low.ND50s[,names(CPE.rem.low.ND50s)%in%c("assay_run1","assay_run2","ND50","ND50_lower","ND50_upper")]*CPE.correction

# add standardized titers (IU/ml): 
CPE.rem.low.ND50s <- CPE.rem.low.ND50s %>% 
  mutate(assay_run1_IU = assay_run1/CPE.rem.low.ND50s$assay_run1[CPE.rem.low.ND50s$sample==stand_sample]*stand_IU,
         assay_run2_IU = assay_run2/CPE.rem.low.ND50s$assay_run2[CPE.rem.low.ND50s$sample==stand_sample]*stand_IU) %>% 
  rowwise() %>% 
  mutate(ND50_IU = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)),na.rm=TRUE)),
         ND50_IU_lower = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2)),
         ND50_IU_upper = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2))) %>%
  ungroup() %>% as.data.frame()


# virus control samples and back titration --------------------------------

CPE.rem.low.1.vc <- data.frame(TCID50s=c((50/250)*(0.025/0.25)^(c(0:5))*0.125*10^3.3, # TCID50 per well in "Back titration David's method", virus control, and samples 29 & 30
                                        ((0.02/0.200)^c(0:5))*0.02*10^3.3), # TCID50 per well in "Back titre + wash"
                              n.pos=c(str_count(as.character(raw.CPE.rem.low.1.data[c(7:36),15]), pattern ="\\+") + # virus control
                                        sum(str_count(as.character(raw.CPE.rem.low.1.data[c(35,36),c(5:14)]), pattern ="\\+")) + # samples 29 & 30
                                        str_count(as.character(raw.CPE.rem.low.1.data[42,4]), pattern ="\\+"), # "Back titration David's method", no dilution
                                      str_count(as.character(raw.CPE.rem.low.1.data[42,c(5:9)]), pattern ="\\+"), # "Back titration David's method", dilutions
                                      str_count(as.character(raw.CPE.rem.low.1.data[41,c(4:9)]), pattern ="\\+")), 
                              n.neg=c(str_count(as.character(raw.CPE.rem.low.1.data[c(7:36),15]), pattern ="-") + # virus control
                                        sum(str_count(as.character(raw.CPE.rem.low.1.data[c(35,36),c(5:14)]), pattern ="-")) + # samples 29 & 30
                                        str_count(as.character(raw.CPE.rem.low.1.data[42,4]), pattern ="-"), # "Back titration David's method", no dilution
                                      str_count(as.character(raw.CPE.rem.low.1.data[42,c(5:9)]), pattern ="-"), # "Back titration David's method", dilutions
                                      str_count(as.character(raw.CPE.rem.low.1.data[41,c(4:9)]), pattern ="-")))
CPE.rem.low.1.vc <- dplyr::mutate(CPE.rem.low.1.vc,n.rep=n.pos+n.neg)
CPE.rem.low.1.vc <- dplyr::mutate(CPE.rem.low.1.vc,dilution=TCID50s/CPE.rem.low.1.vc$TCID50s[1])

CPE.rem.low.2.vc <- data.frame(TCID50s=c((50/250)*(0.025/0.25)^(c(0:5))*0.125*10^3.3, # TCID50 per well in "Back titration David's method", virus control, and samples 29 & 30
                                        ((0.02/0.200)^c(0:5))*0.02*10^3.3), # TCID50 per well in "Back titre + wash"
                              n.pos=c(str_count(as.character(raw.CPE.rem.low.2.data[c(7:36),15]), pattern ="\\+") + # virus control
                                        sum(str_count(as.character(raw.CPE.rem.low.2.data[c(35,36),c(5:14)]), pattern ="\\+")) + # samples 29 & 30
                                        str_count(as.character(raw.CPE.rem.low.2.data[c(43,44),4]), pattern ="\\+"), # "Back titration David's method", no dilution
                                      str_count(as.character(raw.CPE.rem.low.2.data[c(43,44),c(5:9)]), pattern ="\\+"), # "Back titration David's method", dilutions
                                      str_count(as.character(raw.CPE.rem.low.2.data[c(41,42),c(4:9)]), pattern ="\\+")), 
                              n.neg=c(str_count(as.character(raw.CPE.rem.low.2.data[c(7:36),15]), pattern ="-") + # virus control
                                        sum(str_count(as.character(raw.CPE.rem.low.2.data[c(35,36),c(5:14)]), pattern ="-")) + # samples 29 & 30
                                        str_count(as.character(raw.CPE.rem.low.2.data[c(43,44),4]), pattern ="-"), # "Back titration David's method", no dilution
                                      str_count(as.character(raw.CPE.rem.low.2.data[c(43,44),c(5:9)]), pattern ="-"), # "Back titration David's method", dilutions
                                      str_count(as.character(raw.CPE.rem.low.2.data[c(41,42),c(4:9)]), pattern ="-")))
CPE.rem.low.2.vc <- dplyr::mutate(CPE.rem.low.2.vc,n.rep=n.pos+n.neg)
CPE.rem.low.2.vc <- dplyr::mutate(CPE.rem.low.2.vc,dilution=TCID50s/CPE.rem.low.2.vc$TCID50s[1])

# combine virus control to one dataframe:
CPE.rem.low.vc <- rbind(dplyr::mutate(CPE.rem.low.1.vc,assay=rep(1,nrow(CPE.rem.low.1.vc))),
                       dplyr::mutate(CPE.rem.low.2.vc,assay=rep(2,nrow(CPE.rem.low.2.vc))))


# cleanup -----------------------------------------------------------------

rm(raw.CPE.rem.low.1.data,raw.CPE.rem.low.2.data,CPE.rem.low.1.vc,CPE.rem.low.2.vc)


# -------------------------------------------------------------------------
#' 
#' HIGHER VIRAL INOCULUM
#' 
# -------------------------------------------------------------------------


# load raw data -----------------------------------------------------------

# CPE data from the assay with removing the inoculum and higher virus inoculum:
raw.CPE.rem.high.1.data <- read_excel("raw-data/FM-24-031 - CPE MN (WHO IE Discover only; Vic01; wash after adsorption)_2024-07-23.xlsx", sheet = "MN results 1")
raw.CPE.rem.high.2.data <- read_excel("raw-data/FM-24-031 - CPE MN (WHO IE Discover only; Vic01; wash after adsorption)_2024-07-23.xlsx", sheet = "MN results 2")


# process CPE assay data --------------------------------------------------

# CPE ND50s calculated using Reed-Muench method:
CPE.rem.high.ND50s <- left_join(raw.CPE.rem.high.1.data[c(7:36),c(2,3,33)],raw.CPE.rem.high.2.data[c(7:36),c(2,3,33)],by=join_by(`...2`,`...3`))
names(CPE.rem.high.ND50s) <- c("sample","sample_name","assay_run1","assay_run2")
CPE.rem.high.ND50s <- CPE.rem.high.ND50s %>% mutate(across(c(assay_run1,assay_run2), as.numeric))

# replace ND50s below the starting dilution with NA:
CPE.rem.high.ND50s$assay_run1[CPE.rem.high.ND50s$assay_run1 <= 20] <- NA
CPE.rem.high.ND50s$assay_run2[CPE.rem.high.ND50s$assay_run2 <= 20] <- NA

# add CIs:
CPE.rem.high.ND50s <- CPE.rem.high.ND50s %>% 
  rowwise() %>% 
  mutate(ND50 = 10^(mean(log10(c(assay_run1,assay_run2)),na.rm=TRUE)),
         ND50_lower = 10^(mean(log10(c(assay_run1,assay_run2)),na.rm=TRUE)-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)),
         ND50_upper = 10^(mean(log10(c(assay_run1,assay_run2)),na.rm=TRUE)+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1,assay_run2)))/sqrt(2)))

# CPE dilution correction:
CPE.rem.high.ND50s[,names(CPE.rem.high.ND50s)%in%c("assay_run1","assay_run2","ND50","ND50_lower","ND50_upper")] <- 
  CPE.rem.high.ND50s[,names(CPE.rem.high.ND50s)%in%c("assay_run1","assay_run2","ND50","ND50_lower","ND50_upper")]*CPE.correction

# add standardized titers (IU/ml): 
CPE.rem.high.ND50s <- CPE.rem.high.ND50s %>% 
  mutate(assay_run1_IU = assay_run1/CPE.rem.high.ND50s$assay_run1[CPE.rem.high.ND50s$sample==stand_sample]*stand_IU,
         assay_run2_IU = assay_run2/CPE.rem.high.ND50s$assay_run2[CPE.rem.high.ND50s$sample==stand_sample]*stand_IU) %>% 
  rowwise() %>% 
  mutate(ND50_IU = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)),na.rm=TRUE)),
         ND50_IU_lower = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))-qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2)),
         ND50_IU_upper = 10^(mean(log10(c(assay_run1_IU,assay_run2_IU)))+qnorm((1+alpha_CI)/2)*sd(log10(c(assay_run1_IU,assay_run2_IU)))/sqrt(2))) %>%
  ungroup() %>% as.data.frame()


# virus control samples and back titration --------------------------------

# TCID50 is 8-times as high as in the previous CPE update
CPE.rem.high.1.vc <- data.frame(TCID50s=c((50/250)*(0.025/0.25)^(c(0:5))*8*0.125*10^3.3, # TCID50 per well in "Back titration David's method", virus control, and samples 29 & 30
                                         ((0.02/0.200)^c(0:5))*8*0.02*10^3.3), # TCID50 per well in "Back titre + wash"
                               n.pos=c(str_count(as.character(raw.CPE.rem.high.1.data[c(7:36),15]), pattern ="\\+") + # virus control
                                         sum(str_count(as.character(raw.CPE.rem.high.1.data[c(35,36),c(5:14)]), pattern ="\\+")) + # samples 29 & 30
                                         str_count(as.character(raw.CPE.rem.high.1.data[c(43,44),4]), pattern ="\\+"), # "Back titration David's method", no dilution
                                       str_count(as.character(raw.CPE.rem.high.1.data[c(43,44),c(5:9)]), pattern ="\\+"), # "Back titration David's method", dilutions
                                       str_count(as.character(raw.CPE.rem.high.1.data[c(41,42),c(4:9)]), pattern ="\\+")), # "Back titre + wash", dilutions
                               n.neg=c(str_count(as.character(raw.CPE.rem.high.1.data[c(7:36),15]), pattern ="-") + # virus control
                                         sum(str_count(as.character(raw.CPE.rem.high.1.data[c(35,36),c(5:14)]), pattern ="-")) + # samples 29 & 30
                                         str_count(as.character(raw.CPE.rem.high.1.data[c(43,44),4]), pattern ="-"), # "Back titration David's method", no dilution
                                       str_count(as.character(raw.CPE.rem.high.1.data[c(43,44),c(5:9)]), pattern ="-"), # "Back titration David's method", dilutions
                                       str_count(as.character(raw.CPE.rem.high.1.data[c(41,42),c(4:9)]), pattern ="-"))) # "Back titre + wash", dilutions
CPE.rem.high.1.vc <- dplyr::mutate(CPE.rem.high.1.vc,n.rep=n.pos+n.neg)
CPE.rem.high.1.vc <- dplyr::mutate(CPE.rem.high.1.vc,dilution=TCID50s/CPE.rem.high.1.vc$TCID50s[1])

CPE.rem.high.2.vc <- data.frame(TCID50s=c((50/250)*(0.025/0.25)^(c(0:5))*8*0.125*10^3.3, # TCID50 per well in "Back titration David's method", virus control, and samples 29 & 30
                                         ((0.02/0.200)^c(0:5))*8*0.02*10^3.3), # TCID50 per well in "Back titre + wash"
                               n.pos=c(str_count(as.character(raw.CPE.rem.high.2.data[c(7:36),15]), pattern ="\\+") + # virus control
                                         sum(str_count(as.character(raw.CPE.rem.high.2.data[c(35,36),c(5:14)]), pattern ="\\+")) + # samples 29 & 30
                                         str_count(as.character(raw.CPE.rem.high.2.data[c(43,44),4]), pattern ="\\+"), # "Back titration David's method", no dilution
                                       str_count(as.character(raw.CPE.rem.high.2.data[c(43,44),c(5:9)]), pattern ="\\+"), # "Back titration David's method", dilutions
                                       str_count(as.character(raw.CPE.rem.high.2.data[c(41,42),c(4:9)]), pattern ="\\+")), # "Back titre + wash", dilutions
                               n.neg=c(str_count(as.character(raw.CPE.rem.high.2.data[c(7:36),15]), pattern ="-") + # virus control
                                         sum(str_count(as.character(raw.CPE.rem.high.2.data[c(35,36),c(5:14)]), pattern ="-")) + # samples 29 & 30
                                         str_count(as.character(raw.CPE.rem.high.2.data[c(43,44),4]), pattern ="-"), # "Back titration David's method", no dilution
                                       str_count(as.character(raw.CPE.rem.high.2.data[c(43,44),c(5:9)]), pattern ="-"), # "Back titration David's method", dilutions
                                       str_count(as.character(raw.CPE.rem.high.2.data[c(41,42),c(4:9)]), pattern ="-"))) # "Back titre + wash", dilutions
CPE.rem.high.2.vc <- dplyr::mutate(CPE.rem.high.2.vc,n.rep=n.pos+n.neg)
CPE.rem.high.2.vc <- dplyr::mutate(CPE.rem.high.2.vc,dilution=TCID50s/CPE.rem.high.2.vc$TCID50s[1])

# combine virus control to one dataframe:
CPE.rem.high.vc <- rbind(dplyr::mutate(CPE.rem.high.1.vc,assay=rep(1,nrow(CPE.rem.high.1.vc))),
                        dplyr::mutate(CPE.rem.high.2.vc,assay=rep(2,nrow(CPE.rem.high.2.vc))))


# cleanup -----------------------------------------------------------------

rm(CPE.correction,raw.CPE.rem.high.1.data,raw.CPE.rem.high.2.data,CPE.rem.high.1.vc,
   CPE.rem.high.2.vc)

