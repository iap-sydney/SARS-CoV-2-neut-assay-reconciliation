# -------------------------------------------------------------------------
#' 
#' Helper functions for fitting a logistic curve to the Pseudovirus (PV) data
#' and to visualize data and fit
#' 
# -------------------------------------------------------------------------


# initial ND50 guess ------------------------------------------------------

initial.PV.ND50s <- function(fit.data){
  # uses the raw data for fitting to calculate a crude ND50 estimate (for fit to individual assay run)
  # INPUTS:
  # fit.data  raw data for fitting
  # OUTPUT: 
  # ND50s     estimate of initial ND50 parameter values
  
  samples <- unique(fit.data$sample[!fit.data$sample%in%c("CC","VC")])
  ND50s <- rep(NA,length(samples))
  
  for(i in 1:length(samples)){
    tmp.plate <- unique(fit.data$plate[fit.data$sample==samples[i]])
    plate.VC <- mean(fit.data$outcome[fit.data$plate==tmp.plate & fit.data$sample=="VC"])
    
    if(all(fit.data$outcome[fit.data$sample==samples[i]]>=plate.VC)){
      # if all outcomes are above the mean plate VC, then use the lowest dilution as ND50 guess:
      ND50s[i] <- min(fit.data$dilution[fit.data$sample==samples[i]])
    }else{
      # choose the dilution for which the outcome is closest to half the plate.VC:
      ND50s[i] <- fit.data$dilution[which.min(abs(plate.VC/2-fit.data$outcome[fit.data$sample==samples[i]])) + which.min(fit.data$sample==samples[i])-1]
    }
  }
  
  ND50s
}


# least squares fitting function ------------------------------------------

leastsq.PV <- function(data,par){
  # INPUTS:
  # data      dataframe for the PV data
  # par       parameters to be fitted; the order of the parameters is: m for each plate & assay, slope for each/all sample/s, ND50 for each sample
    # m       maximal fluorescence for each plate
    # slope   slope parameter for the logistic function, log-transformed, one slope for all samples
    # ND50    parameter for the logistic function, log-transformed, one for each sample
  # OUTPUT: sum of squared errors for least-squares fit to the data

  n.plates <- length(unique(data$plate))*length(unique(data$assay))
  n.slopes <- 1

  ### control experiments:
  if(length(unique(data$assay))>1){
    tmp.control <- sum((as.numeric(data$outcome[data$sample=="VC"])-par[as.numeric(data$plate[data$sample=="VC"])+max(as.numeric(data$plate))*(as.numeric(data$assay)-1)])^2)
  }else{
    tmp.control <- sum((as.numeric(data$outcome[data$sample=="VC"])-par[as.numeric(data$plate[data$sample=="VC"])])^2)
  }

  ### other samples:

  # only consider samples with corrected dilution information: (i.e., exclude controls)
  data.tmp <- data[!is.na(data$dilution),]

  if(length(unique(data$assay))>1){
    tmp.samples <- sum((as.numeric(data.tmp$outcome)-
                          sapply(c(1:nrow(data.tmp)),function(x){par[as.numeric(data.tmp$plate[x])+max(as.numeric(data$plate))*(as.numeric(data.tmp$assay[x])-1)]* # maximal fluorescence on the same plate
                              log_fun(data.tmp$dilution[x],c(par[n.plates+1],par[n.plates+n.slopes+as.numeric(data.tmp$sample[x])]))}))^2)
  }else{
    tmp.samples <- sum((as.numeric(data.tmp$outcome)-
                          sapply(c(1:nrow(data.tmp)),function(x){par[as.numeric(data.tmp$plate[x])]* # maximal fluorescence on the same plate
                              log_fun(data.tmp$dilution[x],c(par[n.plates+1],par[n.plates+n.slopes+as.numeric(data.tmp$sample[x])]))}))^2)
  }

  ### complete sum of squares:
  tmp.control + tmp.samples
}


# visualization of PV data ----------------------------------------

plot.PV <- function(data,assay,add.fit=FALSE,fit.data=NULL){
  # INPUTS:
  # data        data to be visualized
  # assay       assay number (to restrict the data)
  # add.fit     TRUE or FALSE for adding fit to the data, if TRUE then fit.data needs to be provided (default: FALSE)
  # fit.data    data for the fitted curves (default: NULL)
  # OUTPUT:
  # fig         ggplot figure of data

  fig <- ggplot(data[!is.na(data$dilution) & data$assay==assay,],
                aes(x=as.numeric(dilution), y=as.numeric(outcome), group=rep, color=factor(as.numeric(as.character(sample))))) +
    geom_point() +
    geom_line() +
    # legend, axis, theme, etc.:
    coord_cartesian(ylim=c(0,4.6e5),xlim=c(20,64000)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    labs(x="Dilution", y="RLU (Relative Luminescence Units)", title=paste("PV assay run ",assay,sep="")) +
    theme_bw() +
    theme(legend.position="none") +
    # facet by dose:
    facet_wrap(~ as.numeric(as.character(sample)), ncol=6) +
    # add fit curves to each facet:
    {if(add.fit & !is.null(fit.data)) geom_line(data=fit.data[fit.data$assay==assay,],inherit.aes = FALSE,
                                                aes(x=dilution,y=prob,group=sample),linewidth=1)}

  fig
}

plot.PV.contr <- function(data,assay,add.fit=FALSE,fit.output=NULL){
  # INPUTS:
  # data        data containing control data
  # assay       assay number (to restrict the data)
  # add.fit     TRUE or FALSE for adding fit to the data, if TRUE then fit.output needs to be provided (default: FALSE)
  # fit.output  output from least squares fit to PV data (default: NULL)
  # OUTPUT:
  # fig         ggplot figure of control data

  fig <- ggplot(data[is.na(data$dilution) & data$assay==assay & data$sample!="CC",],
                aes(x=rep(c(1:6),5), y=as.numeric(outcome), group=sample, color=sample)) +
    geom_point() +
    # legend, axis, theme, etc.:
    coord_cartesian(ylim=c(0,4.6e5)) +
    scale_x_continuous(breaks = c(1:6)) +
    labs(x="Number", y="RLU (Relative Luminescence Units)", title="Controls") +
    theme_bw() +
    theme(legend.position="none") +
    # facet by dose:
    facet_wrap(~ plate, ncol=1, labeller = label_both) +
    {if(add.fit & !is.null(fit.output) & length(unique(data$assay))>1) geom_hline(data=data.frame(m=fit.output$estimate[c(1:5)+(assay-1)*5],plate=c(1:5)),
                                                   aes(yintercept=m,group=plate),linewidth=1)} +
    {if(add.fit & !is.null(fit.output) & length(unique(data$assay))==1) geom_hline(data=data.frame(m=fit.output$estimate[c(1:5)],plate=c(1:5)),
                                                   aes(yintercept=m,group=plate),linewidth=1)}

  fig
}
