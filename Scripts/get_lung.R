##########################
# this is the 2nd step in the pipeline - read in raw data - pca, outlier, batch correction.


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
project_folder <- paste0(home_folder, 'pdxSNF/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in id map
##########
id_map <- read.csv(paste0(data_folder, '/pdx_microarray_annotation.csv'), stringsAsFactors = F)
# recall that a.* is a different technology and we gotta do batch correction when combined with other data.

# remove white spaces from array_label
id_map$array_label <- trimws(id_map$array_label, which = "both")


##########
# read in Lung data 
##########

datList <- list()
directory <- dir(micro_array_data)

for (folder in directory){
  # read into list
  datList[[folder]] <- read.csv(paste0(micro_array_data , folder, '/sample_probe_nonnorm_FinalReport.csv'), stringsAsFactors = F)
  print(folder)
}

##########
# loop through data, clean, log2, normalize
##########
for (dat in 1:length(datList)) {
  # get data from list
  data <- datList[[dat]]
  
  data$batch <- dat
  
  # remove first 7 rows
  data <- data[-(1:7),]
  
  # colnames first row, delete first row
  colnames(data) <- data[1,]
  data <- data[-1,]
  
  # recode Target 
  colnames(data)[1] <- 'gene'
  
  # keep data that has colnames 'AVG'
  data <- data[ , grepl('AVG|gene|STDERR|Pval', colnames(data))]
  
  # get ride of 'N' 'T'
  data <- data[ , !grepl('N.AVG|N.BEAD|N.Detection|T.AVG|T.BEAD|T.Detection', colnames(data))]
  
  # transpose
  data <- as.data.frame(t(data), stringsAsFactors = F)
  colnames(data) <- data[1,]
  data <- data[-1,]
  
  # make numeric
  data[, 1:ncol(data)] <- apply(data[, 1:ncol(data)], 2, as.numeric)
  
  for( i in seq(1, nrow(data), by=3) ) {
    
    temp.2 <- data[ i:(i+2), ] 
    
    # subset by sig array
    pval <- temp.2[3,]
    pval_sig_ind <- pval < 0.05
    temp.3 <- temp.2[, pval_sig_ind]

    data[ i:(i+2), ] <- temp.3
    
    print(i)
    
  }
  
  # keep only rows that have avg in it 
  data$model_id <- rownames(data)
  
  data <- data[grepl('AVG', data$model_id),]
  features <- colnames(data)[(1:ncol(data) - 1)]
  data <- data[, c('model_id', features)]
  
  # normalize data
  data[, 2:ncol(data)] <- scale( data[, 2:ncol(data)])
  
  # # apply log2 
  # data[, 2:ncol(data)] <- log2(data[, 2:ncol(data)])
  # 
  # # transpose for normlization
  # data <- t(data)
  # data <- data[-1,]
  # 
  # # get rownmaes, colnames, id
  # row_names <- rownames(data)
  # col_names <- colnames(data)
  # 
  # # make numeric
  # data <- apply(data, 2, function(x) as.numeric(x))
  # 
  # # normalize
  # data <- normalize.quantiles(data)
  # 
  # # add back colnames and rownames
  # rownames(data) <- row_names
  # colnames(data) <- col_names
  # 
  # data <- t(data)
  # 
  datList[[dat]] <- data
  
  print(dat)
  
}


rm(pval, pval_sig_ind, temp.2, temp.3, data)

save.image('/home/benbrew/Desktop/temp_lung.RData')
# load('/home/benbrew/Desktop/temp_lung.RData')

### going forward
# two data sets - dasl and all together
# get intersection of features

##########
# Loop through and get dasl_int and full_int
##########

# fill col_list_full with all column names  from list elements 1-10
col_list_full <- list()
for (dat in 1:length(datList)) {
  data <- datList[[dat]]
  col_list_full[[dat]] <- colnames(data)
}

# fill col_list_dasl with all column names from list elements 2-10
col_list_dasl <- list()
for (dat in 2:length(datList)) {
  data <- datList[[dat]]
  col_list_dasl[[dat]] <- colnames(data)
}


# get intersection of both lists
full_int <- Reduce(intersect, col_list_full)

# remove first element 
col_list_dasl[[1]] <- NULL

dasl_int <- Reduce(intersect, col_list_dasl)

##########
# loop through and subset each data type by 
##########

getSubInt <- function(data_list, full) 
{
  
  if (full == T) {
    
    for (dat in 1:length(data_list)) {
      
      data <- data_list[[dat]]
      
      data <- data[, full_int]
      
      data_list[[dat]] <- data
      
    }
      
  } else {
    
    for (dat in 2:length(data_list)) {
      
      data <- data_list[[dat]]
      
      data <- data[, dasl_int]
      
      data_list[[dat]] <- data
      
      
      
    }
    data_list[[1]] <- NULL
    
  }
  return(data_list)
}
  
# get full data
data_list_full <- getSubInt(datList, full = T)

data_list_dasl <- getSubInt(datList, full = F)

##########
# collapse lists into data frames 
##########
data_full <- do.call(rbind, data_list_full)

data_dasl <- do.call(rbind, data_list_dasl)


##########
# get gene avg
##########
getGene <- function(data, full) 
{
  
  # transpose data
  data <- as.data.frame(t(data), stringsAsFactors = F)
  
  # keep only signal1 column
  data <- data[, grepl('_Signal', colnames(data))]
  
  # get features 
  features <- colnames(data)
  
  # make gene a column
  data$gene <- rownames(data)
  
  # remove every .*
  new_gene <- strsplit(data$gene, '.', fixed = T)
  data$gene <- as.character(lapply(new_gene, function(x) x[1]))

  dat_group <- data %>% 
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  # transpose again
  data <- as.data.frame(t(dat_group), stringsAsFactors = F)
  
  # make first rown colnames and remove
  colnames(data) <- as.character(data[1,])
  data <- data[-1,]
  
  # make column for id 
  data$id <- rownames(data)
  
  # rearrange 
  features <- colnames(data[, (1:ncol(data) -1)])
  
  data <- data[ , c('id', features)]
  
  ids <- strsplit(data$id, '.', fixed = T)
  data$id <- unlist(lapply(ids, function(x) x[1]))
  
  
  return(data)
}

data_dasl <- getGene(data_dasl)
data_full <- getGene(data_full)

#########
# create new colum in id_map that can join both dasl and hybrid
#########

# grap first 42 of array_barcode
array_bar <- id_map$array_barcode[1:42]

# put first 42 of array_barcod into array_label
id_map$array_label[1:42] <- array_bar

#########
# clean id column and join data with id_map
#########
cleanMap <- function(data) 
{
  features <- colnames(data)[2:ncol(data)]
  data$id <- as.factor(data$id)
  id_map$array_label <- as.factor(id_map$array_label)
  data <- inner_join(data, id_map, c("id" = "array_label") )
  data <- data[, c("id" ,"model_id", "array_barcode", "passage" ,"project_title" , features)]
  
  return(data)
}

data_dasl <- cleanMap(data_dasl)
data_full <- cleanMap(data_full)

#########
# drop duplicates
#########
treatDups <- function(data, full)
{
  
  data$array_barcode <- NULL
  data$id <- NULL
  # first make numeric 
  data[,4:ncol(data)] <- apply(data[,4:ncol(data)], 2, as.numeric)
  
  # make group by variables factors
  data$model_id <- as.character(data$model_id)
  data$passage <- as.character(data$passage)
  data$project_title <- as.character(data$project_title)
  
  # make factors 
  if(full) {
    batch_factor <- c("a", "b", "c", "d", "e", 
                      "f", "g", "h", "i", "j")
    
    data$project_title <- batch_factor[as.factor(data$project_title)]
  } else {
    batch_factor <- c("b", "c", "d", "e", 
                      "f", "g", "h", "i", "j")
    
    data$project_title <- batch_factor[as.factor(data$project_title)]
  }
  
  # sort data
  data <- data[order(data$model_id, decreasing = F),]
  
  # avg genes that have same model, passage, and project_title
  data_avg <- data %>% #144
    group_by(model_id, passage, project_title) %>%
    summarise_each(funs(mean))
  
  # remove perfect duplicates (probably first few coulmns)
  data_dup <- data_avg[,4:100]
  data_dup <- data_avg[!duplicated(data_dup),]
  
  # remove !complete.cases
  data_dup <- data_dup[complete.cases(data_dup),]
  data_dup <- as.data.frame(data_dup)
  
  return(data_dup)
  
}

data_dasl <- treatDups(data_dasl, full = F)
data_full <- treatDups(data_full, full = T)



# 

##########
# pca 
##########

getPCA <- function(pca_data, 
                   column_name, 
                   gene_start, 
                   name, 
                   full) 
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
  #fill in factors with colors 
  if (full) {
    pca_data[, column_name] <- ifelse(pca_data[, column_name] == 'a', 'a', 'b')
    col_vec <- c('red', 'green')
  } else {
    col_vec <- c('red', 'green', 'blue', 'brown', 'black', 
                 'orange', 'yellow', 'grey', 'lightblue', 'lightgreen')
  }
  
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

getPCA(data_full, 'project_title', 4, 'full data', full = T)
getPCA(data_dasl, 'project_title', 4, 'dasl data', full = F)


##########
# batch 
##########
getCombat <- function(data, full)
{
  
  if (full) {
    data$project_title <- ifelse(data$project_title == 'a', 'a', 'b')
    
    dup_index <- duplicated(data$model_id)
    data$model_id[dup_index] <- paste0(data$model_id[dup_index], '_dup')
    
    dup_index <- duplicated(data$model_id)
    data$model_id[dup_index] <- paste0(data$model_id[dup_index], '_dup2')
    
    dup_index <- duplicated(data$model_id)
    data$model_id[dup_index] <- paste0(data$model_id[dup_index], '_dup3')
    
    dup_index <- duplicated(data$model_id)
    data$model_id[dup_index] <- paste0(data$model_id[dup_index], '_dup4')
    
    dup_index <- duplicated(data$model_id)
    data$model_id[dup_index] <- paste0(data$model_id[dup_index], '_dup5')
    
  } else {
    # add unique indicator to duplicated model ids - multiple duplicates - this will do for now
    dup_index <- duplicated(data$model_id)
    data$model_id[dup_index] <- paste0(data$model_id[dup_index], '_dup')
  }
 
  
  # get passage to combine later
  passage <- data$passage
  
  # put model model_ids in rownames and remove columns
  rownames(data) <- data$model_id
  mat_data <- data[, 4:ncol(data)]
  # get features 
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # # full
  # if(full) {
  #   data$batch <- ifelse(data$)
  # }
  
  # get intercept
  modcombat <- model.matrix(~1, data = data)
  batch <- as.factor(data$project_title)
  combat <- ComBat(dat = mat_data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$model_id <- rownames(final_dat)
  final_dat$project_title <- batch
  final_dat$passage <- passage
  final_dat <- final_dat[, c('model_id', 'passage', 'project_title', features)]
  final_dat$project_title <- as.character(final_dat$project_title)
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
}

data_dasl <- getCombat(data_dasl, full = F)

getPCA(data_dasl, 'project_title', 4, 'dasl data', full = F)


#########
# put new batch corrected dasl into full data
#########

getFull <- function(dasl, full) 
{
  # sort each by project_title (batch)
  dasl <- dasl[order(dasl$project_title, decreasing = F),]
  full <- full[order(full$project_title, decreasing = F),]
  
  # get overlapping intersection
  int_feat <- Reduce(intersect, list(colnames(dasl), colnames(full)))
  
  # subset both 
  dasl <- dasl[, int_feat]
  full <- full[, int_feat]
  
  # check to see if dasl is same in both
  stopifnot(all(dasl$project_title  == full$project_title[full$project_title != 'a']))
  
  # put dasl in where project_title is not a
  full[full$project_title != 'a',] <- dasl
  
  return(full)
  
}

data_full <- getFull(data_dasl, data_full)

#########
# apply combat to full data with 
#########

# first pca
getPCA(data_full, 'project_title',4, 'full data, dasl corrected', full = T)

# apply combat
data_full <- getCombat(data_full, full = T)

# get pca again
getPCA(data_full, 'project_title',4, 'full data, dasl corrected', full = T)


#########
# save both data sets
#########
saveRDS(data_full, paste0(data_folder, '/data_full.rds'))
saveRDS(data_dasl, paste0(data_folder, '/data_dasl.rds'))


