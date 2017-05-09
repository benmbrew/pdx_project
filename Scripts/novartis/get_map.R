##########################
# this is the 2nd step in the pipeline - get mapping between rna_seq and microarray for pdxe


##########
# initialize libraries
##########
library(dplyr)
library(lumi)
library(preprocessCore)
library(caret)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdx_project/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# load cleaned full novartis data (lung) and transform with gene_model
##########
rna_full <-readRDS(paste0(data_folder, '/nov_dat_scaled.rda'))

##########
# read in novartis data - list of rna_seq, microarray, and id map between the two
##########
dat_list <- readRDS(paste0(data_folder, '/novartis', '/pdx_rna_seq_microarray.rda'))

# rna_seq 
rna_dat <- dat_list[[1]]

# microarray 
micro_dat <- dat_list[[2]]

# id map 
id_map <- dat_list[[3]]

##########
# transpose matrices and explore
##########
# rna_dat
rna_dat <- as.data.frame(t(rna_dat), stringsAsFactors = F)

# microarray
micro_dat <- as.data.frame(t(micro_dat), stringsAsFactors = F)

##########
# scale data
##########
rna_scaled <- as.data.frame(scale(rna_dat), stringsAsFactors = F)

micro_scaled <- as.data.frame(scale(micro_dat), stringsAsFactors = F)

##########
# get model ids from rna to micro
##########

getId <- function(data, data_type)
{
  if (data_type == 'rna') {
    
    id_map$id <- id_map$biobase.id.RNASeq
    
  } 
  
  if (data_type == 'microarray') {
    
    id_map$id <- id_map$biobase.id.microArray
    
  } 
    data$id <- rownames(data)
    data <- left_join(id_map, data, by = 'id')
    return(data)
}

# scaled
micro_scaled <- getId(micro_scaled, data_type = 'microarray')
rna_scaled <- getId(rna_scaled, data_type = 'rna')

# ##########
# # for each model (id), plot all gene values against each other
# ##########
# ids <- unique(id_map$biobase.id.RNASeq)
# 
# pdf("/home/benbrew/Desktop/models_plot")
# 
# # loop through each id and 
# for(i in ids) {
#  
#   # subset to individual model for scaled data
#   rna_scale <- rna_scaled[rna_scaled$id == i, 5:ncol(rna_scaled)]
#   micro_scale <- micro_scaled[micro_scaled$biobase.id.RNASeq == i, 5:ncol(micro_scaled)]
#   
#   smoothScatter(rna_scale, micro_scale, main = paste0('model', ' ', i, 'scale'),
#                 xlab = 'rna', ylab = 'micro', xlim= c(-5, 10), ylim = c(-5,10))
#   
#   print(i)
#   
# }
# 
# dev.off()

# remove unnedded data
rm(micro_dat, rna_dat, dat_list)

##########
# subset to lung
##########
getLung <- function(data) {
  data <- data[grepl('NSCLC', data$tumor.type),]
  
  return(data)
}

# apply function
micro_scaled <- getLung(micro_scaled)
rna_scaled <- getLung(rna_scaled)

##########
# get common genes between the 3 datasets
##########

# get common features
rna_feat <- colnames(rna_scaled)[5:ncol(rna_scaled)]
rna_feat_full <- colnames(rna_full)[5:ncol(rna_full)]

# get intesection
intersected_feats <- intersect(rna_feat, rna_feat_full)

# subset 3 datasets by intersected features
micro_scaled <- micro_scaled[, c('tumor.type', 'biobase.id.RNASeq', 'biobase.id.microArray', 'id', intersected_feats)]
rna_scaled <- rna_scaled[, c('tumor.type', 'biobase.id.RNASeq', 'biobase.id.microArray', 'id', intersected_feats)]
rna_full <- rna_full[, c('model.id', 'biobase.id', 'tumor.type', 'patient.id', intersected_feats)]

  
##########
# estimate a linear model for each gene between rna and microarray
##########

gene_model <- list()
gene_control_result <- list()

# changed control into original, change rna into micro
for(i in 5:ncol(rna_scaled)) {
  
  rna <- as.data.frame(rna_scaled[, i])
  micro <- as.data.frame(micro_scaled[, i])
  model_data <- data.frame(rna = rna, micro = micro)
  names(model_data) <- c('rna', 'micro')
  gene_model[[i]] <- lm(micro ~ rna, data = model_data)
  
  print(i)
  
}


# now loop throgh rna_full

# # transpose results
# temp <- do.call(rbind, gene_control_result)
# transform_controls <- t(temp)
# 
# # add cg sites
# colnames(transform_controls) <- colnames(beta_overlap[3:ncol(beta_overlap)])
# 
# # add clinical variables
# transform_controls <- as.data.frame(cbind(id = beta_overlap$id, 
#                                           age_sample_collection = beta_overlap$age_sample_collection, 
#                                           transform_controls))



