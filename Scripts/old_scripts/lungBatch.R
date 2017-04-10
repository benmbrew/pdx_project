##########################
# this script will read in lung pdx data and address batch effects

##########
# initiate libraries
##########
library(dplyr)
library(lumi)
library(sva)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdxSNF/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# source functions
##########
source(paste0(code_folder, 'functions.R'))

##########
# read in normalized data
##########

# save data
quan <- readRDS(paste0(data_folder, 'quan.lumi.rda'))
quan_full <- readRDS(paste0(data_folder, 'quan.lumi_full.rda'))

rsn <- readRDS(paste0(data_folder, 'rsn.lumi.rda'))
rsn_full <- readRDS(paste0(data_folder, 'rsn.lumi_full.rda'))

ssn <- readRDS(paste0(data_folder, 'ssn.lumi.rda'))
ssn_full <- readRDS(paste0(data_folder, 'ssn.lumi_full.rda'))

loess <- readRDS(paste0(data_folder, 'loess.lumi.rda'))
loess_full <- readRDS(paste0(data_folder, 'loess.lumi_full.rda'))

rank <- readRDS(paste0(data_folder, 'rank.lumi.rda'))
rank_full <- readRDS(paste0(data_folder, 'rank.lumi_full.rda'))

##########
# create new factor for batches
##########
getBatch <- function(data)
{
  batch_factor <- c("a", "b", "c", "d", "e", 
                    "f", "g", "h", "i", "j")
  
  data$project_title <- batch_factor[as.factor(data$project_title)]
  
  
  return(data)
}

quan <- getBatch(quan)
quan_full <- getBatch(quan_full)

rsn <- getBatch(rsn)
rsn_full <- getBatch(rsn_full)

ssn <- getBatch(ssn)
ssn_full <- getBatch(ssn_full)

loess <- getBatch(loess)
loess_full <- getBatch(loess_full)

rank <- getBatch(rank)
rank_full <- getBatch(rank_full)

##########
# plot pca for batches
##########
# 15.130873
# pca_data <- quan_new
# column_name <- 'batch'
# function needs to take a clinical column, remove others, and plot pcas

getPCA <- function(pca_data, column_name, gene_start, name, full) 
{
  # get features sites
  genes <- colnames(pca_data)[gene_start:ncol(pca_data)]
  
  if (full){
    pca_data <- pca_data[!is.na(pca_data$model_id),]
  }
  
  # put column name with genes
  pca_data <- pca_data[ ,c(column_name, genes)]
  
  # run pca
  data_length <- ncol(pca_data)
  
  # make data numeric 
  pca_data[, 2:ncol(pca_data)] <- apply(pca_data[, 2:ncol(pca_data)], 2, as.numeric)
  pca <- prcomp(pca_data[,2:data_length])
  
  # plot data
  #fill in factors with colors 
  if (full) {
    pca_data[, column_name] <- ifelse(pca_data[, column_name] == 'a', 'a', 'b')
    col_vec <- c('red', 'green')
  } else {
    col_vec <- c('red', 'green', 'blue', 'brown', 'black', 
                 'orange', 'yellow', 'grey', 'lightblue', 'lightgreen')
  }
  
  colors <- col_vec[as.factor(pca_data[, column_name])]
  min_x <- min(pca$x[,1])
  max_x <- max(pca$x[,1])
  min_y <- min(pca$x[,2])
  max_y <- max(pca$x[,2])
  
  plot <- plot(pca$x[,1], 
               pca$x[,2],
               xlab = 'pca 1',
               ylab = 'pca 2',
               cex = 1,
               main = name,
               pch = 16,
               xlim = c(min_x, max_x),
               ylim = c(min_y, max_y),
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}

##########
# examine pca plots of all data
##########

# quan
getPCA(quan, column_name = 'project_title', 6, name = 'quan',  full = F)
getPCA(quan_full, column_name = 'project_title',6, name = 'quan full',  full = T)

# rsn
getPCA(rsn, column_name = 'project_title',6, name = 'rsn',  full = F)
getPCA(rsn_full, column_name = 'project_title',6, name = 'rsn full',  full = T)

# ssn
getPCA(ssn, column_name = 'project_title', 6,name = 'ssn',  full = F)
getPCA(ssn_full, column_name = 'project_title', 6,name = 'ssn full',  full = T)

# loess
getPCA(loess, column_name = 'project_title', 6,name = 'loess',  full = F)
getPCA(loess_full, column_name = 'project_title', 6,name = 'loess full',  full = T)

# rank
getPCA(rank, column_name = 'project_title', 6,name = 'rank',  full = F)
getPCA(rank_full, column_name = 'project_title', 6,name = 'rank full',  full = T)


##########
# run combat to correct for batch effects in between a and everything else
##########
quan_new <- getCombat(quan_full)
rsn_new <- getCombat(rsn_full)
ssn_new <- getCombat(ssn_full)
loess_new <- getCombat(loess_full)
rank_new <- getCombat(rank_full)

# # check data
# max_list <- list()
# for (i in 3:ncol(quan_new)) {
#   max_list[[i]] <- quan_new$model_id[quan_new[,i] == max(quan_new[, i])]
#   print(i)
# }
# 
# temp <- do.call(rbind, max_list)
# temp <- as.data.frame(temp)
# temp$V1 <- as.character(temp$V1)
# temp$V2 <- as.character(temp$V2)

##########
# run pca again to see batch correction
##########
getPCA(quan_new, column_name = 'batch', 3, 'quan corrected', full = T)
getPCA(rsn_new, column_name = 'batch', 3, 'rsn corrected', full = T)
getPCA(ssn_new, column_name = 'batch', 3, 'ssn corrected', full = T)
getPCA(loess_new, column_name = 'batch', 3, 'loess corrected', full = T)
getPCA(rank_new, column_name = 'batch', 3, 'rank corrected', full = T)

##########
# remove outlier - model_id = 110_dup
##########
removeOutlier <- function(data) 
{
  data <- data[data$model_id != '110_dup',]
  return(data)
}

quan_new <- removeOutlier(quan_new)
ssn_new <- removeOutlier(ssn_new)
rsn_new <- removeOutlier(rsn_new)
loess_new <- removeOutlier(loess_new)
rank_new <- removeOutlier(rank_new)

###########
# apply PCA again
###########

getPCA(quan_new, column_name = 'batch', 3, 'quan corrected', full = T)
getPCA(rsn_new, column_name = 'batch', 3, 'rsn corrected', full = T)
getPCA(ssn_new, column_name = 'batch', 3, 'ssn corrected', full = T)
getPCA(loess_new, column_name = 'batch', 3, 'loess corrected', full = T)
getPCA(rank_new, column_name = 'batch', 3, 'rank corrected', full = T)

##########
# save data
##########
saveRDS(quan_new, paste0(data_folder, 'quan_new.rda'))
saveRDS(rsn_new, paste0(data_folder, 'rsn_new.rda'))
saveRDS(ssn_new, paste0(data_folder, 'ssn_new.rda'))
saveRDS(loess_new, paste0(data_folder, 'loess_new.rda'))
saveRDS(rank_new, paste0(data_folder, 'rank_new.rda'))


# 
# # get matrix for arvind
# quan <- readRDS(paste0(data_folder, 'quan_new.rda'))
# 
# # remove batch column
# quan$batch <- NULL
# 
# 
# # put model id in rownmaes and remove column
# rownames(quan) <- quan$model_id
# 
# # remove model id
# quan$model_id <- NULL
# 
# 
# # transpose data
# quan <- as.data.frame(t(quan))
# 
# write.csv(quan, '/home/benbrew/Desktop/lung_processed.csv')
# 
# 
