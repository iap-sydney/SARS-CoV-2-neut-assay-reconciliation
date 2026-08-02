# -------------------------------------------------------------------------
#' 
#' Helper functions to fit a logistic curve to the IS assay data
#' and to visualize the data with fit
#' 
# -------------------------------------------------------------------------


# Least squares fit function ----------------------------------------------

leastsq <- function(data,par,spots.transf="linear"){
  # INPUTS:
  # data          data frame for the IS data
  # par           parameters to be fitted; the order of the parameters is: m for each plate, slope for each/all sample/s, ND50 for each sample
  # m             maximal number of spots on each plate
  # slope         slope parameter for the logistic function, log-transformed, one slope for all samples
  # ND50          parameter for the logistic function, log-transformed, one for each sample
  # spots.transf  fit to the number of spots ("linear") or log(number of spots) ("log")
  # OUTPUT: sum of squared errors for least-squares fit to the data
  
  n.plates <- length(unique(data$plate))*length(unique(data$assay)) # number of plates for each assay times the number of assays
  n.slopes <- 1
  min.assay.num <- min(as.numeric(data$assay))
  
  ### control experiments:
  if(spots.transf=="linear"){
    tmp.control <- sum((as.numeric(data$spots[data$sample=="VOC"])-par[as.numeric(data$plate[data$sample=="VOC"])+(as.numeric(data$assay[data$sample=="VOC"])-min.assay.num)*n.plates/length(unique(data$assay))])^2)
  }else if(spots.transf=="log"){
    tmp.control <- sum((log(as.numeric(data$spots[data$sample=="VOC"]))-par[as.numeric(data$plate[data$sample=="VOC"])+(as.numeric(data$assay[data$sample=="VOC"])-min.assay.num)*n.plates/length(unique(data$assay))])^2)
  }
  
  ### other samples:
  
  # only consider samples with corrected dilution information: (i.e., exclude samples 31 to 35 and controls)
  data.tmp <- data[!is.na(data$corrected.dilution),]
  
  if(spots.transf=="linear"){
    tmp.samples <- sum((as.numeric(data.tmp$spots)-
                          sapply(c(1:nrow(data.tmp)),function(x){par[as.numeric(data.tmp$plate[x])+(as.numeric(data.tmp$assay[x])-min.assay.num)*n.plates/length(unique(data$assay))]* # maximal number of spots on the same plate
                              log_fun(data.tmp$corrected.dilution[x],c(par[n.plates+1],par[n.plates+n.slopes+as.numeric(data.tmp$sample[x])]))}))^2)
    
  }else if(spots.transf=="log"){
    tmp.samples <- sum((log(as.numeric(data.tmp$spots))-
                          sapply(c(1:nrow(data.tmp)),function(x){par[as.numeric(data.tmp$plate[x])+(as.numeric(data.tmp$assay[x])-min.assay.num)*n.plates/length(unique(data$assay))]* # maximal number of spots on the same plate
                              log_fun(data.tmp$corrected.dilution[x],c(par[n.plates+1],par[n.plates+n.slopes+as.numeric(data.tmp$sample[x])]))}))^2)
  }
  
  ### complete sum of squares:
  tmp.control + tmp.samples
}


# visualization of data with fit ------------------------------------------

plot.IS.fit <- function(fit.data,assay,curve.data,num.col=7){
  # Visualization of the IS data with the fit to the data
  # INPUTS:
  # fit.data    data used for fitting
  # assay       assay number (to restrict the data)
  # curve.data  fitted curve data (dataframe with columns: dilution, prob, sample, assay)
  # num.col     number of columns for the samples
  # OUTPUT:
  # fig         ggplot figure of data used for fit and fit to the data

  all.dilutions <- sort(unique(fit.data$corrected.dilution[!is.na(fit.data$corrected.dilution)]))

  if("sample_name"%in%names(fit.data)){
    # helper function for labelling the panels in the facet plot with the "sample_name":
    label_facet <- function(original_names, new_names){
      names(new_names) <- original_names
      return(new_names)
    }

    fig <- ggplot(fit.data[!is.na(fit.data$corrected.dilution) & fit.data$assay==assay,],
                  aes(x=as.numeric(corrected.dilution),y=as.numeric(spots),color=sample,group=assay.plate.sample.rep)) +
      # add data:
      geom_point() +
      geom_line() +
      # legend, axis, theme, etc.:
      {if(length(all.dilutions)<=10) scale_x_log10(breaks = all.dilutions)} +
      {if(length(all.dilutions)>10) scale_x_log10()}  +
      coord_cartesian(ylim=c(0,max(as.numeric(fit.data$spots))+8), xlim=c(min(all.dilutions),max(all.dilutions))) +
      labs(x="Dilution", y="Number of spots", title=paste("Fit to IS assay ",assay-1,sep="")) +
      theme_bw() +
      theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      # facet by dose:
      facet_wrap(~ sample, ncol=num.col, labeller = labeller(sample = label_facet(unique(fit.data$sample[!is.na(fit.data$corrected.dilution) & fit.data$assay==assay]),
                                                                                  unique(fit.data$sample_name[!is.na(fit.data$corrected.dilution) & fit.data$assay==assay])))) +
      # add fit curves to each facet:
      geom_line(data=curve.data[curve.data$assay==assay,],inherit.aes = FALSE,aes(x=dilution,y=prob,group=sample),linewidth=1)
  }else{
    fig <- ggplot(fit.data[!is.na(fit.data$corrected.dilution) & fit.data$assay==assay,],
                  aes(x=as.numeric(corrected.dilution),y=as.numeric(spots),color=sample,group=assay.plate.sample.rep)) +
      # add data:
      geom_point() +
      geom_line() +
      # legend, axis, theme, etc.:
      {if(length(all.dilutions)<=10) scale_x_log10(breaks = all.dilutions)} +
      {if(length(all.dilutions)>10) scale_x_log10()}  +
      coord_cartesian(ylim=c(0,max(as.numeric(fit.data$spots))+8), xlim=c(min(all.dilutions),max(all.dilutions))) +
      labs(x="Dilution", y="Spot count", title=paste("Fit to IS assay run ",assay-1,sep="")) +
      theme_bw() +
      theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      # facet by dose:
      facet_wrap(~ sample, ncol=num.col) +
      # add fit curves to each facet:
      geom_line(data=curve.data[curve.data$assay==assay,],inherit.aes = FALSE,aes(x=dilution,y=prob,group=sample),linewidth=1)
  }

  fig
}

plot.IS.fit.VOC <- function(fit.data,assay,fit.output){
  # Visualization of the IS VOC (virus only control) data with the fit to the data
  # INPUTS:
  # fit.data    data used for fitting
  # assay       assay number (to restrict the data)
  # fit.output  output from least squares fit to IS data
  # OUTPUT:
  # fig         ggplot figure of VOC data with fit

  n.plates <- length(unique(fit.data$plate[fit.data$assay==assay]))
  min.assay.num <- min(as.numeric(fit.data$assay))
  VOC.per.plate <- sum(fit.data$sample[fit.data$assay==assay]=="VOC")/n.plates # assuming there is the same number of VOCs per plate

  fig <- ggplot(fit.data[fit.data$sample=="VOC" & fit.data$assay==assay,],
                aes(x=rep(c(1:VOC.per.plate),n.plates),y=as.numeric(spots))) +
    # add data:
    geom_point(color="gray50") +
    # legend, axis, theme, etc.:
    coord_cartesian(ylim=c(0,max(as.numeric(fit.data$spots))+8)) +
    scale_x_continuous(breaks = 2*c(1:ceiling(VOC.per.plate/2))) +
    labs(x="Repeats", y="Spot count", title=paste("Controls",sep="")) +
    theme_bw() +
    theme(legend.position="none") +
    # facet by plate:
    facet_wrap(~ plate, ncol=1, labeller = label_both) +
    # add fit curves to each facet:
    geom_hline(data=data.frame(m=fit.output$estimate[c(1:n.plates)+(assay-min.assay.num)*length(unique(fit.data$plate[fit.data$assay==assay]))],plate=c(1:n.plates)),
               aes(yintercept=m,group=plate),linewidth=1)

  fig
}
