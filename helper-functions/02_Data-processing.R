# -------------------------------------------------------------------------
#' 
#' Functions for data processing, i.e. formatting the data in a way that 
#' facilitates analysis.
#' 
# -------------------------------------------------------------------------

data_processing <- function(raw.data,experiment,dilutions=c(20,40,80,160,320,640),dilution_correction=1,samples='Neut ID'){
  # Processes IS or CPE data
  # INPUTS:
  # raw.data     data from a spreadsheet to be formatted
  # experiment   specify the experiment: "IS" or "CPE"
  # dilutions    specify dilutions used in this experiment
  # dilution.correction   corrects the convolution from the spreadsheet (default=1, i.e. no correction) 
  # samples      use 'Neut ID' (default) or 'Sample' to populate the sample column in the processed data
  # OUTPUTS:
  # data         formatted data
  
  ### create a data frame for the data:
  if(experiment=="IS"){
    # number of plates:
    n.plates <- sum(!is.na(unique(as.numeric(as.matrix(raw.data[,which(raw.data=="Plate #",arr.ind = T)[2]])))))
    # data frame:
    data <- data.frame(plate=rep(c(1:n.plates),rep(96,n.plates)), # plates are numbered from 1
                       row=rep(rep(c("A","B","C","D","E","F","G","H"),rep(6,8)),n.plates*2),
                       column=rep(c(rep(c(1:6),8),rep(c(7:12),8)),n.plates),
                       sample=factor(c(sapply(c(1:n.plates),function(x){c((x-1)*7+rep(c(1:7),rep(12,7)),rep(c(rep("VOC",5),"NVC"),2))})),levels = c(1:35,"VOC","NVC")),
                       rep=rep(rep(c(1,2),c(6,6)),8*n.plates),
                       dilution=rep(dilutions,16*n.plates),
                       spots=rep(NA,n.plates*96))
  }else if(experiment=="CPE"){
    if(!any(names(raw.data)=="Sample")){
      # identify row to use for column names and set column names
      names(raw.data) <- as.character(as.matrix(raw.data[which(raw.data=="Sample",arr.ind = T)[1],]))
      # subset data to rows with sample information:
      raw.data <- raw.data[!is.na(raw.data$`Neut ID`) & raw.data$`Neut ID`!="Neut ID",]
    }
    if(samples=='Neut ID'){
      all.samples <- unique(raw.data$`Neut ID`[!is.na(raw.data$`Neut ID`)])
    }else if(samples=='Sample'){
      all.samples <- unique(raw.data$Sample[!is.na(raw.data$Sample)])
    }
    n.samples <- length(all.samples)
    CPE.dilutions <- as.numeric(names(raw.data)[c((1+which(names(raw.data)=="Starting dilution")):(which(names(raw.data)=="Virus control")-1))])
    n.plates <- length(unique(raw.data$`Plate #`)[!is.na(raw.data$`Plate #`)])
    data <- data.frame(sample=rep(all.samples,rep(length(CPE.dilutions),n.samples)),
                       dilution=rep(CPE.dilutions,n.samples),
                       plate=rep(c(1:n.plates),rep(2*length(CPE.dilutions),n.plates)),
                       n.pos=rep(NA,length(CPE.dilutions)*n.samples),
                       n.neg=rep(NA,length(CPE.dilutions)*n.samples))
  }
  
  ### complete data frame
  if(experiment=="IS"){
    index.row <- which(raw.data=="Row",arr.ind = T)[2] # column in raw.data with 'row' (of plate) information
    index.col <- which(raw.data=="Row",arr.ind = T)[1] # row in raw.data with the 'column' (of plate) information
    # fill number of spots:
    for(i in 1:nrow(data)){
      data$spots[i] <- raw.data[which(raw.data[,index.row]==data$row[i])[data$plate[i]], # select correct row and plate
                                which(as.numeric(raw.data[index.col,])==data$column[i])] # select correct column
    }
    # set dilution for control experiments to NA:
    data$dilution[data$sample%in%c("VOC","NVC")] <- NA
  }else if(experiment=="CPE"){
    # add number of positive and negative wells for each sample, dilution, and assay:
    for(i in c(1:nrow(data))){
      if(samples=='Neut ID'){
        data$n.pos[i] <- str_count(as.character(raw.data[which(raw.data$`Neut ID`==data$sample[i]),names(raw.data)==data$dilution[i] & !is.na(names(raw.data))]), pattern ="\\+")
        data$n.neg[i] <- str_count(as.character(raw.data[raw.data$`Neut ID`==data$sample[i],names(raw.data)==data$dilution[i] & !is.na(names(raw.data))]), pattern ='-')
      }else if(samples=="Sample"){
        data$n.pos[i] <- str_count(as.character(raw.data[which(raw.data$Sample==data$sample[i]),names(raw.data)==data$dilution[i] & !is.na(names(raw.data))]), pattern ="\\+")
        data$n.neg[i] <- str_count(as.character(raw.data[raw.data$Sample==data$sample[i],names(raw.data)==data$dilution[i] & !is.na(names(raw.data))]), pattern ='-')
      }
    }
    # add total number of replicates:
    data <- data %>% mutate(n.rep=n.pos+n.neg)
  }
  
  ### dilution correction:
  data$dilution <- data$dilution*dilution_correction
  
  ### output:
  data
}


# loading fit data --------------------------------------------------------

load_fit_data <- function(assay_name, file_path = 'output/Fitting'){
  # INPUTS:
  # file_path   specifies the file location
  # assay_name  name of the assay to extract fit data: "PV-fit", "IS-fit", etc.
  #             Note that if the assay name contains brackets "(",")", they need to be written as "\\(", "\\)"
  # OUTPUT:
  # tempfit     extracted tempfit file

  if(length(list.files(path = file_path,pattern = assay_name, full.names = FALSE))>0){
    latest <- max(gsub("[^0-9]", "",str_extract(sub(paste0(".*",assay_name), "", list.files(path = file_path,pattern = assay_name, full.names = FALSE)), "\\d.*\\d")))
    latest <- unlist(strsplit(latest, split = ""))
    latest <- paste0(paste(latest[c(1:4)],collapse = ''),'-',paste(latest[c(5:6)],collapse = ''),'-',paste(latest[c(7:8)],collapse = ''),'_',paste(latest[c(9:14)],collapse = ''),collapse ='')
    tmp <- load(glue::glue("output/Fitting/{list.files(path = file_path,pattern = paste0(assay_name, '.*', latest), full.names = FALSE)}"))
  }else{
    print("Error: Missing fit data")
  }

  # return fit data:
  return(get(tmp))
}

