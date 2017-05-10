##########################
# this is the 3rd step in the pipeline - get novartis, lung batch. read in novartis id_map, lung


##########
# initialize libraries
##########
library(dplyr)
library(lumi)
library(preprocessCore)
library(caret)
library(sva)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdx_project/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in novartis
##########
nov <- readRDS(paste0(data_folder, 'nov_micro_transformed.rda'))

nov <- as.data.frame(nov, stringsAsFactors= F)

##########
# read in lung data
##########
data_full <- readRDS(paste0(data_folder, '/data_full.rds'))
data_dasl <- readRDS(paste0(data_folder, '/data_dasl.rds'))

##########
# correct for batch between lung and novartis
##########

# # get common features 
# data <- data_dasl
# novartis_data <- nov
commonFeat <- function(data, novartis_data) 
{
  # get common features
  features <- colnames(data)[4:ncol(data)]
  novartis_features <- colnames(novartis_data)[3:ncol(novartis_data)]
  intersect_feat <- intersect(features, novartis_features)
  
  # fix lumi, subset
  data$passage =NULL
  data$project_title = NULL
  names(data)[1] <- 'id'
  data <- data[, c('id', intersect_feat)]
  
  # fix novartis, subset
  novartis_data$id <- novartis_data$patient.id
  novartis_data <- novartis_data[, c('id', intersect_feat)]
  
  full_data <- rbind(data, novartis_data)
  
  return(full_data)
  
  
}

# for normal
data_full <- commonFeat(data_full, nov)

data_dasl <- commonFeat(data_dasl, nov)

# get batch variable
data_full$batch <- ifelse(grepl('X', data_full$id), 'nov', 'lung')
data_dasl$batch <- ifelse(grepl('X', data_dasl$id), 'nov', 'lung')

# get id and batch at front 
features <- colnames(data_full)[2:(ncol(data_full) - 1)]
data_full <- data_full[, c('id', 'batch', features)]
data_dasl <- data_dasl[, c('id', 'batch', features)]



# pca 

# pca_data <- data_full
# column_name = 'batch'
# gene_start = 3
# name = 't'
getPCA <- function(pca_data, column_name, gene_start, name) 
{
  # get features sites
  genes <- colnames(pca_data)[gene_start:ncol(pca_data)]
  
  # if (full){
  #   pca_data <- pca_data[!is.na(pca_data$model_id),]
  # }
  
  # put column name with genes
  pca_data <- pca_data[ ,c(column_name, genes)]
  
  # make data numeric 
  pca <- prcomp(pca_data[,2:ncol(pca_data)])
  
  # plot data
 
  col_vec <- c('red', 'green')
  
  
  colors <- col_vec[as.factor(pca_data[, column_name])]
  
  
  plot <- plot(pca$x[,1], 
               pca$x[,2],
               xlab = 'pca 1',
               ylab = 'pca 2',
               cex = 1,
               main = name,
               pch = 16,
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}

getPCA(data_full, 'batch', 3, 'Novartis and Mings Data Together')
getPCA(data_dasl, 'batch', 3, 'dasl data')

##########
# explore difference in ming and novartis
##########



##########
# scale data
##########
data_dasl[, 3:ncol(data_dasl)] <- scale(data_dasl[, 3:ncol(data_dasl)] )
data_full[, 3:ncol(data_full)] <- scale(data_full[, 3:ncol(data_full)] )



##########
# batch 
##########
getCombat <- function(data)
{
  
  data$id <- make.names(data$id, unique=TRUE)
  # put model model_ids in rownames and remove columns
  rownames(data) <- data$id
  mat_data <- data[, 3:ncol(data)]
  # get features 
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # # full
  # if(full) {
  #   data$batch <- ifelse(data$)
  # }
  
  # get intercept
  modcombat <- model.matrix(~1, data = data)
  batch <- as.factor(data$batch)
  combat <- ComBat(dat = mat_data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$model_id <- rownames(final_dat)
  final_dat$batch <- batch
  final_dat <- final_dat[, c('model_id', 'batch', features)]
  final_dat$batch <- as.character(final_dat$batch)
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
}

data_dasl <- getCombat(data_dasl)
data_full <- getCombat(data_full)


getPCA(data_full, 'batch', 3, 'full data')
getPCA(data_dasl, 'batch', 3, 'dasl data')

#########
# save both data sets
#########
saveRDS(data_full, paste0(data_folder, '/data_full.rds'))
saveRDS(data_dasl, paste0(data_folder, '/data_dasl.rds'))

