# -------------------------------------------------------------------------
#' 
#' General helper functions
#' 
# -------------------------------------------------------------------------


# Logistic function -------------------------------------------------------

# probability of neutralizing a virion at titre/dilution x (or the probability
# of not neutralizing a virion at dilution x (concentration is original concentration/x))
# logistic function with maximum 1 and 2 parameters (slope, IC50):
# parameter 1: slope, log-transformed
# parameter 2: ND50, log-transformed
# the titre/dilution x is not log-transformed!
log_fun <- function(x,par){1/(1+exp(-exp(par)[1]*(log10(x)-log10(exp(par[2])))))}
log_fun <- Vectorize(log_fun, vectorize.args = "x")


# AIC for least-squares fit -----------------------------------------------

# based on Banks & Joyner (2017) "AIC under the Framework of Least Squares Estimation"
AIC_ls <- function(leastsquares,n.data,n.par){n.data*log(leastsquares/n.data) + 2*(n.par+1)}


# prediction of ND50s -----------------------------------------------------

# predict the CPE ND50s based on the IS data fit and the number of virions
predict_ND50s <- function(n_virions,n_virions_var=NULL,IS_slope=IS_fit_slope,IS_slope_var=IS_fit_slope_var_log,
                          xmin_log10=0.5,xmax_log10=4.5,alph=alpha_CI,n_boostr=1e5){
  # INPUTS:
  # n_virions             number of virions per well in the CPE assay
  # n_virions_var         variance of the number of virions per well, if there is no variance: n_virions_var = NULL
  # IS_slope              estimated slope of the logistic function fit to the IS data
  # IS_slope_var          variance of the slope for the IS data, on a log-scale
  # xmin_log10            minimal x-position (IS ND50) on a log10-scale
  # xmax_log10            maximal x-position (IS ND50) on a log10-scale
  # alph                  confidence level
  # n_boostr              number of bootstraps
  # OUTPUT:
  # prediction            data frame with IS_ND50 (x-position), CPE_ND50 (prediction), CPE_ND50_lower, and CPE_ND50_upper

  # set up data frame:
  x_log10 <- seq(xmin_log10,xmax_log10,length.out=1e3)
  prediction <- data.frame(IS_ND50 = 10^x_log10)

  # add prediction for CPE ND50s:
  prediction <- mutate(prediction, CPE_ND50 = 10^(x_log10+log(2^(1/n_virions)-1)/IS_slope))

  # Confidence interval on the prediction:
  # bootstrap number of virions and slopes:
  if(is.null(n_virions_var)){
    n.vir.tmp <- rep(n_virions,n_boostr) # no variance in n_virions
  }else{
    n.vir.tmp <- rnorm(n_boostr,mean=n_virions,sd=sqrt(n_virions_var)) # sample the number of virions per well
  }
  slope.tmp <- exp(rnorm(n_boostr,mean=log(IS_slope),sd=sqrt(IS_slope_var))) # IS slope variance is on a log-scale
  # make prediction and take quantiles:
  tmp_intercept_log10 <- quantile(log(2^(1/n.vir.tmp)-1)/slope.tmp,c((1-alph)/2,(1+alph)/2),na.rm = TRUE)
  prediction <- mutate(prediction, CPE_ND50_lower = 10^(x_log10+min(tmp_intercept_log10)),
                       CPE_ND50_upper = 10^(x_log10+max(tmp_intercept_log10)))

  # output:
  prediction
}

# predict the CPE ND50s based on the IS data fit and the number of virions
# predict_ND50s_by_virions <- function(IS_ND50,IS_ND50_var,IS_slope=IS_slope,IS_slope_var=IS_slope_var,
#                                      min_n_virions=1,max_n_virions=4e4,alph=alpha_CI,n_boostr=2e5){
#   # INPUTS:
#   # IS_ND50       IS ND50
#   # IS_ND50_var   variance of the IS ND50, on a log-scale
#   # IS_slope      estimated slope of the logistic function fit to the IS data
#   # IS_slope_var  variance of the slope for the IS data, on a log-scale
#   # min_n_virions         minimal number of virions
#   # max_n_virions         maximal number of virions
#   # alph                  confidence level
#   # n_boostr              number of bootstraps
#   # OUTPUT:
#   # prediction            dataframe with n_virions (x-position), CPE_ND50 (prediction), CPE_ND50_lower, and CPE_ND50_upper
#   
#   # set up data frame:
#   n_virions_log2 <- seq(log2(min_n_virions),log2(max_n_virions),length.out=1e3)
#   prediction <- data.frame(n_virions = 2^n_virions_log2)
#   
#   # add prediction for CPE ND50s:
#   prediction <- mutate(prediction, CPE_ND50 = 10^(log10(IS_ND50)+log(2^(1/prediction$n_virions)-1)/IS_slope))
#   
#   # Confidence interval on the prediction:
#   # bootstrap IS ND50 and slopes:
#   ND50.tmp <- exp(rnorm(n_boostr,mean=log(IS_ND50),sd=sqrt(IS_ND50_var))) # IS ND50 variance is on a log-scale
#   slope.tmp <- exp(rnorm(n_boostr,mean=log(IS_slope),sd=sqrt(IS_slope_var))) # IS slope variance is on a log-scale
#   # make prediction and take quantiles:
#   tmp_lower <- sapply(c(1:nrow(prediction)),function(x){quantile(10^(log10(ND50.tmp)+log(2^(1/prediction$n_virions[x])-1)/slope.tmp),c((1-alph)/2),na.rm = TRUE)})
#   tmp_upper <- sapply(c(1:nrow(prediction)),function(x){quantile(10^(log10(ND50.tmp)+log(2^(1/prediction$n_virions[x])-1)/slope.tmp),c((1+alph)/2),na.rm = TRUE)})
#   prediction <- mutate(prediction, CPE_ND50_lower = tmp_lower, CPE_ND50_upper = tmp_upper)
#   
#   # output:
#   prediction
# }


# visualization: plot vs IS ND50s ---------------------------------

plot_ND50_comp <- function(assay_x = "IS", assay_y, prediction = NULL, show_1to1 = show_1to1_relationship,
                           show_linear = show_linear_fit, data = ND50s, alph = alpha_CI, stand = standardize){
  # INPUTS:
  # assay_x     name of the assay to be shown on the x-axis, default: "IS"
  # assay_y     name of the assay to be shown on the y-axis: "Pseudoneut","CPE_original","CPE_rem_low","CPE_rem_high","CPE_replace","CPE_all_wash"
  # prediction  data for the prediction of the ND50 based on the IS ND50 (generated with "prediction_ND50s"), if prediction = NULL it is not added to the plot
  # show_1to1   show the 1:1 relationship in the plot
  # show_linear show linear fit to the data
  # data        data for plotting (with both the assay and the IS ND50s)
  # alph        confidence level
  # stand       standardized titers (relative to relative to NIBSC 21/338 titer (sample 28))
  # OUTPUT:
  # fig         figure comparing the assay ND50s to the IS ND50s

  ### rename data:
  # assay_x:
  names(data)[names(data)==paste(assay_x,"_ND50",sep="")] <- "x_ND50"
  names(data)[names(data)==paste(assay_x,"_ND50_lower",sep="")] <- "x_lower"
  names(data)[names(data)==paste(assay_x,"_ND50_upper",sep="")] <- "x_upper"
  assay_x_name <- assay_x
  assay_x_details <- ""
  if(grepl("CPE",assay_x)){
    assay_x_details <- gsub("CPE_","",assay_x_name)
    assay_x_name <- "CPE"
    if(assay_x_details=="original"){
      assay_x_details <- " (original)"
    }else if(assay_x_details=="rem_low"){
      assay_x_details <- " (washing, low inoculum)"
    }else if(assay_x_details=="rem_high"){
      assay_x_details <- " (washing, high inoculum)"
    }else if(assay_x_details=="replace"){
      assay_x_details <- " (replacing media, high inoculum)"
    }else if(assay_x_details=="all_wash"){
      assay_x_details <- " (all assay_ys with washing/replacing media)"
    }
  }else if(assay_x=="PV"){
    assay_x_name <- "PV"
  }
  # assay_y:
  names(data)[names(data)==paste(assay_y,"_ND50",sep="")] <- "y_ND50"
  names(data)[names(data)==paste(assay_y,"_ND50_lower",sep="")] <- "y_lower"
  names(data)[names(data)==paste(assay_y,"_ND50_upper",sep="")] <- "y_upper"
  assay_y_name <- assay_y
  assay_y_details <- ""
  if(grepl("CPE",assay_y)){
    assay_y_details <- gsub("CPE_","",assay_y_name)
    assay_y_name <- "CPE"
    if(assay_y_details=="original"){
      assay_y_details <- " (without inoc. removal)"
    }else if(assay_y_details=="rem_low"){
      assay_y_details <- " (inoc. removal, low virus inoc.)"
    }else if(assay_y_details=="rem_high"){
      assay_y_details <- " (inoc. removal, high virus inoc.)"
    }
  }else if(assay_y=="PV"){
    assay_y_name <- "PV"
  }

  ### compute the linear fit if it is to be shown in the figure
  if(show_linear){
    # linear model with estimated slope:
    lm_mod <- lm(log10(y_ND50) ~ log10(x_ND50), data=data) # run linear model
    if(stand){
      x_log10 <- seq(log10(0.01*stand_IU),log10(6*stand_IU),length.out=1e3)
    }else{
      x_log10 <- seq(0.5,4.5,length.out=1e3)
    }
    lm_mod_confb <- predict(lm_mod, data.frame(x_ND50=10^x_log10), se.fit = TRUE, interval = "confidence", level = alph)
    lm_mod_region <- data.frame(x_ND50=10^x_log10, estim = 10^(lm_mod_confb$fit[,1]), lower=10^lm_mod_confb$fit[,2], upper=10^lm_mod_confb$fit[,3])
  }

  ### make the figure
  show_annotations <- TRUE
  annot_pos <- ifelse(stand,0.2*stand_IU,10^2.75)
  fig <- ggplot(data,aes(x=x_ND50,y=y_ND50)) +
    # add 1:1 line:
    {if(show_1to1_relationship) geom_abline(slope = 1, intercept = 0, color="forestgreen",linetype=2, linewidth=1.5)} +
    {if(show_1to1_relationship & !stand) annotate("text",x=10^2.5,y=10^2.75,label="1:1 relationship", hjust=0, color="forestgreen")} +
    {if(show_1to1_relationship & stand) annotate("text",x=0.5*stand_IU,y=0.5*stand_IU,label="1:1 relationship", hjust=0, color="forestgreen")} +
    # add linear fit:
    {if(show_linear) geom_ribbon(data=lm_mod_region,inherit.aes=FALSE,aes(x=x_ND50,ymin=lower, ymax=upper),fill="black", alpha = 0.1)} +
    {if(show_linear) geom_abline(slope = lm_mod$coef[2],intercept = lm_mod$coef[1], color="black")} +
    {if(show_linear & show_annotations) annotate("text",x=annot_pos,y=annot_pos,label=paste("Best fit (log10-scale):"), hjust=0)} +
    {if(show_linear & show_annotations) annotate("text",x=annot_pos,y=annot_pos/1.15,label=paste("Intercept: ",round(lm_mod$coef[1],3)," (95% CI: ",round(confint(lm_mod)[1,1],3)," to ",round(confint(lm_mod)[1,2],3),")",sep=""), hjust=0)} +
    {if(show_linear & show_annotations) annotate("text",x=annot_pos,y=annot_pos/1.15/1.15,label=paste("Slope: ",round(lm_mod$coef[2],3)," (95% CI: ",round(confint(lm_mod)[2,1],3)," to ",round(confint(lm_mod)[2,2],3),")",sep=""), hjust=0)} +
    # add prediction:
    {if(!is.null(prediction)) geom_ribbon(data=prediction, inherit.aes = FALSE, aes(x = IS_ND50, ymin=CPE_ND50_lower, ymax=CPE_ND50_upper),fill="orangered", alpha = 0.2)} +
    {if(!is.null(prediction)) geom_line(data=prediction, inherit.aes = FALSE, aes(x = IS_ND50, y=CPE_ND50), color="orangered",linetype=3,linewidth=1)} +
    {if(!is.null(prediction) & show_annotations) annotate("text",x=4e3,y=10^2.5,label=paste("Prediction (log10-scale):"), hjust=0, color="orangered")} +
    {if(!is.null(prediction) & show_annotations) annotate("text",x=4e3,y=10^2.5/1.15,label=paste("Intercept: ",unique(round(log10(prediction$CPE_ND50)-log10(prediction$IS_ND50),3)),
                                                                              " (95% CI: ",unique(round(log10(prediction$CPE_ND50_lower)-log10(prediction$IS_ND50),3)),
                                                                              " to ",unique(round(log10(prediction$CPE_ND50_upper)-log10(prediction$IS_ND50),3)),")",sep=""), hjust=0, color="orangered")} +
    {if(!is.null(prediction) & show_annotations) annotate("text",x=4e3,y=10^2.5/1.15/1.15,label="Slope: 1", hjust=0, color="orangered")} +
    # add data with errorbars:
    geom_point() +
    geom_errorbarh(aes(xmin=x_lower,xmax=x_upper),width=0.03) +
    geom_errorbar(aes(ymin=y_lower,ymax=y_upper),width=0.03) +
    # legend, axis, theme, etc.:
    {if(!stand) scale_x_log10(breaks=2^(c(3:14)))} +
    {if(!stand) scale_y_log10(breaks=2^(c(3:14)))} +
    {if(stand) scale_x_log10(breaks=40*2^c(0:9))} +
    {if(stand) scale_y_log10(breaks=40*2^c(0:9))} +
    coord_cartesian(xlim=c(2^floor(min(log2(data$x_lower),na.rm = T)),2^ceiling(max(log2(data$x_upper),na.rm = T))),
                    ylim=c(2^floor(min(log2(data$y_lower),na.rm = T)),2^ceiling(max(log2(data$y_upper),na.rm = T)))) +
    labs(x=glue::glue("{assay_x_name}{ifelse(assay_x_name=='CPE',assay_x_details,'')} ND50{ifelse(stand,' (IU/ml)','')}"),
         y=glue::glue("{assay_y_name}{ifelse(assay_y_name=='CPE',assay_y_details,'')} ND50{ifelse(stand,' (IU/ml)','')}"),
         title=glue::glue("ND50 comparison: {assay_y_name}{assay_y_details} vs {assay_x_name}{assay_x_details}")) +
    theme_bw() +
    theme(legend.position="none",axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  fig
}
