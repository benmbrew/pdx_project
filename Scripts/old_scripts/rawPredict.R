##########################
# this script will read in lung pdx data and run SNF
# look at class
# different responses
# check novartis unique models
##########
# initiate libraries
##########
library(dplyr)
library(lumi)
library(SNFtool)
library(caret)
library(randomForest)
library(Metrics)
library(sva)
library(randomForest)
library(kernlab)
library(pROC)
library(doParallel)
library(nnet)

##########
# initiate folder
##########
home_folder = '/home/benbrew/Documents/'
project_folder = paste0(home_folder, 'pdxSNF/')
data_folder = paste0(project_folder, 'Data/')
micro_array_data = paste0(data_folder, 'microarray/')
code_folder = paste0(project_folder, 'Code/')

##########
# source functions
##########
source(paste0(code_folder, 'functions.R'))

##########
# read in id map
##########
id_map <- read.csv(paste0(data_folder, '/pdx_microarray_annotation.csv'), stringsAsFactors = F)
# recall that a.* is a different technology and we gotta do batch correction when combined with other data.

# remove white spaces from array_label
id_map$array_label <- trimws(id_map$array_label, which = "both")

##########
# read in data - raw data and novartis data
##########
quan_lumi_norm <- readRDS(paste0(data_folder, 'quan.lumi.rda'))

a.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 10.09.17 (DirectHyb lungXeno 2 Raphael)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)


b.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 11.11.16 Set1 (DASL lungXeno)',
                              '/sample_probe_nonnorm_FinalReport.csv'), 
                       stringsAsFactors = F)


c.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 11.11.16 Set2 (DASL lungXeno)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)


d.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 12.01.05 set1 (DASL 9 Notch 3Lung 12Pan)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)


e.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 12.01.05 set2  (DASL 24 lung Xeno)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

f.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 12.06.26 (DASL Part1 LungXeno)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

g.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 12.06.26 Part2 (DSSL 16 LungXeno + 8 Niki)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

h.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 12.08.15 (DASL Lisa Notch Xeno)',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

i.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 12.11.30_Lisa',
                              '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)


j.lumi_raw <- read.csv(paste0(micro_array_data, 
                              'Tsao 13.02.07_PHLC164_ORF',
                              '/sample_probe_nonnorm_FinalReport.csv'), 
                       stringsAsFactors = F)



# read in novartis data
mrna_seq <- read.csv(paste0(data_folder, '/novartis', '/rna_seq.csv'))

# # read in house keeping genes
# house_keep <- read.csv(paste0(data_folder, '/HK_exons.csv'))
# 
# # get number of hk genes
# length(unique(house_keep$Gene.Name))
# 
# # get unique gene
# hk_genes <- unique(house_keep$Gene.Name)
  

##########
# function that cleans raw data
##########

sepScaleData <- function(data)
{
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
  
  num_models <- nrow(data)/3
  
  # make numeric
  data[, 1:ncol(data)] <- apply(data[, 1:ncol(data)], 2, as.numeric)
  
  data_list <- list()
  for( i in seq(1, nrow(data), by=3) ) {
    data_list[[i]] <- data[ i:(i+2), ] 
    temp.2 <- data_list[[i]]
    # subset by sig array
    pval <- temp.2[3,]
    pval_sig_ind <- pval < 0.05
    temp.3 <- temp.2[, pval_sig_ind]
    mean_sig <- mean(as.numeric(temp.3[1,]))
    mean_values <- sweep(temp.3[1,1:ncol(temp.3)], 2, mean_sig, `-`)    
    scaled_values <- mean_values/temp.3[2,]
    temp.4 <- rbind(temp.3, scaled_values)
    
    # temp.4 <- as.data.frame(t(temp.3))
    model_id <- strsplit(rownames(temp.4)[1], '.', fixed = T )
    temp.4$model_id <- unlist(model_id)[[1]]
    features <- colnames(temp.4)[(1:ncol(temp.4) - 1)]
    temp.4 <- temp.4[, c('model_id', features)]
    data_list[[i]] <- temp.4
  
  }
  
  # 
  # for (ids in 1:num_groups){
  #   # get num_ids data sets
  #   d <- data[, ]
  # }
  # 
  # 
  # # get only significant genes
  # data <- data
  # 
  # remove variance 0 (maybe after fkpm)
  
  # normalize to rkpm/fkpm
  
  # scale data
  
  
  
  return(data_list)
}

# a.lumi_sep
a.lumi_sep <- sepScaleData(a.lumi_raw)
b.lumi_sep <- sepScaleData(b.lumi_raw)
c.lumi_sep <- sepScaleData(c.lumi_raw)
d.lumi_sep <- sepScaleData(d.lumi_raw)
e.lumi_sep <- sepScaleData(e.lumi_raw)
f.lumi_sep <- sepScaleData(f.lumi_raw)
g.lumi_sep <- sepScaleData(g.lumi_raw)
h.lumi_sep <- sepScaleData(h.lumi_raw)
i.lumi_sep <- sepScaleData(i.lumi_raw)
j.lumi_sep <- sepScaleData(j.lumi_raw)

rm(a.lumi_raw, b.lumi_raw, c.lumi_raw,
   d.lumi_raw, e.lumi_raw, f.lumi_raw,
   g.lumi_raw, i.lumi_raw, j.lumi_raw,
   h.lumi_raw)

   

##########
# function to subset intersectio features and combine data
##########
combineData <- function(sep_list) 
{
  features_list <- list()
  
  # first remove nulls from sep_list
  new_list <- Filter(Negate(function(x) is.null(unlist(x))), sep_list)
  
  # nod get intersection
  for(dat in 1:length(new_list)) {
    temp.1 <- new_list[[dat]]
    features_list[[dat]] <- colnames(temp.1)[2:ncol(temp.1)]
  }
  

  intersecting_features <- Reduce(intersect, features_list)
   
  
  for(dat in 1:length(new_list)) {
     temp.1 <- new_list[[dat]]
     temp.2 <- temp.1[, c('model_id', intersecting_features)]
     new_list[[dat]] <- temp.2
  }
     
  combined_data <- do.call(rbind, new_list)
   
  return(combined_data)
}

a.lumi_com <- combineData(a.lumi_sep)
b.lumi_com <- combineData(b.lumi_sep)
c.lumi_com <- combineData(c.lumi_sep)
d.lumi_com <- combineData(d.lumi_sep)
e.lumi_com <- combineData(e.lumi_sep)
f.lumi_com <- combineData(f.lumi_sep)
g.lumi_com <- combineData(g.lumi_sep)
i.lumi_com <- combineData(i.lumi_sep)
j.lumi_com <- combineData(j.lumi_sep)

#########
# transpose, get gene column, and avg over genes
#########
arrangeData <- function(data, dasl) 
{
  
  # transpose data
  data <- as.data.frame(t(data), stringsAsFactors = F)
  
  # keep only signal1 column
  data <- data[, grepl('_Signal1', colnames(data))]
  
  # make gene a column
  data$gene <- rownames(data)
  
  # make first rown colnames and remove
  colnames(data) <- as.character(data[1,])
  data <- data[-1,]
  
  # remove every .*
  new_gene <- strsplit(data$model_id, '.', fixed = T)
  data$gene <- as.character(lapply(new_gene, function(x) x[1]))
  data$model_id <- NULL
  
  # make data frame and numeric
  data <- as.data.frame(data)
  data[, 1:(ncol(data)- 1)] <- as.data.frame(apply(data[, 1:(ncol(data)- 1)], 2, as.numeric))
  
  # group by gene and summarise
  data <- data %>%
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  # transpose
  data <- as.data.frame(t(data), stringsAsFactors= F)
  
  # colnames first row and remove first row
  colnames(data) <- data[1,]
  data <- data[-1,]
  
  
  
  return(data)
}
  

a.lumi_fin <- arrangeData(a.lumi_com, dasl = F)
b.lumi_fin <- arrangeData(b.lumi_com, dasl = T)
c.lumi_fin <- arrangeData(c.lumi_com, dasl = T)
d.lumi_fin <- arrangeData(d.lumi_com, dasl = T)
e.lumi_fin <- arrangeData(e.lumi_com, dasl = T)
f.lumi_fin <- arrangeData(f.lumi_com, dasl = T)
g.lumi_fin <- arrangeData(g.lumi_com, dasl = T)
i.lumi_fin <- arrangeData(i.lumi_com, dasl = T)
j.lumi_fin <- arrangeData(j.lumi_com, dasl = T)

##########
# combine data: two sets: normal, full
##########
aggData <- function(full)
{
  lumi_intersection <- Reduce(intersect, list(colnames(b.lumi_fin)[5:ncol(b.lumi_fin)], 
                                              colnames(c.lumi_fin)[5:ncol(c.lumi_fin)],
                                              colnames(d.lumi_fin)[5:ncol(d.lumi_fin)],
                                              colnames(e.lumi_fin)[5:ncol(e.lumi_fin)],
                                              colnames(f.lumi_fin)[5:ncol(f.lumi_fin)],
                                              colnames(g.lumi_fin)[5:ncol(g.lumi_fin)],
                                              colnames(i.lumi_fin)[5:ncol(i.lumi_fin)],
                                              colnames(j.lumi_fin)[5:ncol(j.lumi_fin)]))

  b.lumi <- b.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  c.lumi <- c.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  d.lumi <- d.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  e.lumi <- e.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  f.lumi <- f.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  g.lumi <- g.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  # h.lumi_fin <- h.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  i.lumi <- i.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  j.lumi <- j.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
  
  lumi <- rbind(b.lumi, c.lumi,d.lumi,
                e.lumi, f.lumi, g.lumi, 
                i.lumi, j.lumi)
  
  if (full) {
    lumi_full_intersection <- Reduce(intersect, list(colnames(a.lumi_fin)[5:ncol(a.lumi_fin)], 
                                                colnames(lumi)[5:ncol(lumi)]))
    
    a.lumi <- a.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    b.lumi <- b.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    c.lumi <- c.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    d.lumi <- d.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    e.lumi <- e.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    f.lumi <- f.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    g.lumi <- g.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    # h.lumi_fin <- h.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_intersection)]
    i.lumi <- i.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    j.lumi <- j.lumi_fin[, c('array_barcode', 'model_id', 'passage', 'project_title', lumi_full_intersection)]
    
    lumi <- rbind(a.lumi, b.lumi, c.lumi,
                  d.lumi, e.lumi, f.lumi,
                  g.lumi, i.lumi, j.lumi)
  }
  return(lumi)
                              
}

lumi <- aggData(full = F)
lumi_full <- aggData(full = T)

##########
# Treat for dups
##########
treatDups <- function(data, full)
{
  
  data$array_barcode <- NULL
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
  data_avg <- data %>%
    group_by(model_id, passage, project_title) %>%
    summarise_each(funs(mean))
  
  # remove perfect duplicates (probably first few coulmns)
  data_dup <- data_avg[4:20]
  data_dup <- data_avg[!duplicated(data_dup),]
  
  # remove !complete.cases
  data_dup <- data_dup[complete.cases(data_dup),]
  data_dup <- as.data.frame(data_dup)
  
  return(data_dup)
  
}

lumi <- treatDups(lumi, full = F) #formally 103, 7985
lumi_full <- treatDups(lumi_full, full = T) # 12, 1265

getPCA(lumi, column_name = 'project_title', gene_start = 4, name = 'DASL', full = F)
getPCA(lumi_full, column_name = 'project_title', gene_start = 4, name = 'DASL + DirectHyb', full = T)

##########
# batch correction
##########
lumi_correct <- getCombat(lumi, full = F)
lumi_full_correct <- getCombat(lumi_correct, full = T, a_lumi = a.lumi_fin)

# check pca
getPCA(lumi_correct, column_name = 'batch', gene_start = 4, name = 'DASL Corrected', full = F)
getPCA(lumi_full_correct, column_name = 'batch', gene_start = 3, name = 'DASL + DirectHyb Corrected', full = T)


##########
# clean novartis data
##########
# transpose
novartis <- as.data.frame(t(mrna_seq), stringsAsFactors = F)

# put first row in colnames then delete first row
colnames(novartis) <- novartis[1,]
novartis <- novartis[-1,]

# make id column
novartis$id <- rownames(novartis)

# rearrange 
features <- colnames(novartis[, (1:ncol(novartis) - 1)])
novartis <- novartis[, c('id', features)]
rownames(novartis) <- NULL

##########
# apply log2 to novaris data
##########
# make numeric
novartis[, 2:ncol(novartis)] <- apply(novartis[, 2:ncol(novartis)], 2, as.numeric)

# zero var
ids <- novartis$id
temp.nov <- novartis[, 2:ncol(novartis)]
near_zero <- nearZeroVar(temp.nov)
remove_these <- colnames(temp.nov)[near_zero]
keep_these <- !colnames(temp.nov) %in%  remove_these
temp.nov_2 <- temp.nov[, keep_these]
novartis <- temp.nov_2
# log2
novartis[, 1:ncol(novartis)] <- log(novartis[, 1:ncol(novartis)] + .1, base=2)

# scale novartis 
for (i in 1:ncol(novartis)) {
  sub_col <- novartis[, i]
  mean_col <- mean(sub_col)
  sd_col <- sd(sub_col)
  novartis[, i] <- (sub_col - mean_col)/sd_col
  print(i)
}

# add back id
novartis_features <- colnames(novartis)
novartis$id <- ids
novartis <- novartis[, c('id', novartis_features)]

# combine novartis and mings data 
##########

commonFeat <- function(lumi_data, novartis_data) 
{
  # get common features
  lumi_features <- colnames(lumi_data)[4:ncol(lumi_data)]
  novartis_features <- colnames(novartis_data)[2:ncol(novartis_data)]
  intersect_feat <- intersect(lumi_features, novartis_features)
  
  # fix lumi, subset
  lumi_data$passage =NULL
  lumi_data$project_title = NULL
  names(lumi_data)[1] <- 'id'
  lumi_data <- lumi_data[, c('id', intersect_feat)]
  
  # fix novartis, subset
  novartis_data <- novartis_data[, c('id', intersect_feat)]
  
  return(list(lumi_data, novartis_data))
  
  
}

# for normal
com_feat_norm <- commonFeat(lumi_correct, novartis)
lumi <- com_feat_norm[[1]]
novartis <- com_feat_norm[[2]]

dat_norm <- rbind(lumi, novartis)

# for full
com_feat_full <- commonFeat(lumi_full_correct, novartis)
lumi_full <- com_feat_full[[1]]
novartis_full <- com_feat_full[[2]]

dat_full <- rbind(lumi_full, novartis_full)

##########
# PCA
##########

# create batches 
dat_norm$batch <- ifelse(grepl('X.', dat_norm$id), 'b', 'a')
dat_full$batch <- ifelse(grepl('X.', dat_full$id), 'b', 'a')

# rearrange
dat_norm <- dat_norm[, c('id', 'batch', colnames(dat_norm[2:(ncol(dat_norm) - 1)]))]
dat_full <- dat_full[, c('id', 'batch', colnames(dat_full[2:(ncol(dat_full) - 1)]))]

getPCA(dat_norm, column_name = 'batch', gene_start = 3, name = 'DASL + Novartis', full = F)
getPCA(dat_full, column_name = 'batch', gene_start = 3, name = 'DASL + DirectHyb + Novartis', full = T)

# batch
dat_norm_cor <- getCombatNov(dat_norm)
dat_full_cor <- getCombatNov(dat_full)

removeOutlier <- function(data, full) 
{
  if(full) {
    data <- data[data$id != '66',]
    data <- data[data$id != '58',]
    data <- data[data$id != '110',]
    
  } else {
    data <- data[data$id != '303',]
    data <- data[data$id != '274',]
    data <- data[data$id != '110',]
    
  }

  return(data)
}


dat_norm_cor <- removeOutlier(dat_norm_cor, full = F)
dat_full_cor <- removeOutlier(dat_full_cor, full = T)


# PCA
getPCA(dat_norm_cor, column_name = 'batch', gene_start = 3, name = 'DASL + Novartis Corrrected', full = F)
getPCA(dat_full_cor, column_name = 'batch', gene_start = 3, name = 'DASL + DirectHyb + Novartis Corrected', full = T)


##########
# split data back up
##########
novartis <- dat_norm_cor[grepl('X.', dat_norm_cor$id),]
novartis_full <- dat_full_cor[grepl('X.', dat_full_cor$id),]


# save.image('/home/benbrew/Desktop/temp.RData')
# load('/home/benbrew/Desktop/temp.RData')

##########

# data <- dat_full_cor
# list_result <- list()
# for(i in 3:ncol(data)) {
#   list_result[[i]] <- data$id[max(data[,i], na.rm = T)]
#   print(i)
# }
# 
# temp <- as.data.frame(do.call(rbind, list_result))
# result <- temp %>% group_by(V1) %>% summarise(counts = n())
# result <- result[order(result$counts, decreasing = T),]
# ##########
# read novartis raw outcome
##########
# raw_outcome <- read.csv(paste0(data_folder, '/raw_outcome.csv'), stringsAsFactors = F)
# 
# # subset by young cancer NSCLC
# raw_nsclc <- raw_outcome[raw_outcome$Tumor.Type == 'NSCLC',]
# 
# # define outcome as max - min for grouped model, treatment
# raw_avg_outcome <- raw_nsclc %>%
#   group_by(Model, Treatment) %>%
#   summarise(max_size = max(Volume..mm3.),
#             min_size = min(Volume..mm3.))
# 
# # define outcome variable as difference between max and min
# raw_avg_outcome$response <- raw_avg_outcome$max_size - raw_avg_outcome$min_size
# 
# ##########
# # combine Novartis data 
# ##########
# # recode raw_avg_outcome from '-' to '.' for join
# raw_avg_outcome$Model <- gsub('-', '.', raw_avg_outcome$Model,  fixed = T)
# 
# # recode to id 
# names(raw_avg_outcome)[1] <- 'id'
# 
# # join novartis on id
# outcome_inner <- inner_join(novartis, raw_avg_outcome, by = 'id')
# 
# # rearrange 
# features <- colnames(outcome_inner)[2:(ncol(outcome_inner) - 4)]
# 
# outcome <- outcome_inner[, c('id', 'Treatment', 'max_size', 'min_size', 'response', features)]
# 
# # order outcome
# nov_mod_data <- outcome[order(outcome$Treatment),]

# ##########
# # find intersection of novartis and mings data
# ##########
# # get features
# nov_feat <- colnames(nov_mod_data[, 6:ncol(nov_mod_data)])
# lumi_feat <- colnames(lumi_correct[, 2:ncol(lumi_correct)])
# lumi_full_feat <- colnames(lumi_correct[, 2:ncol(lumi_full_correct)])
# 
# 
# # intersect
# mod_feat <- intersect(nov_feat, lumi_feat)
# mod_full_feat <- intersect(nov_feat, lumi_full_feat)
# 
# 
# # subset model_data
# model_data <- nov_mod_data[,c('id', 'Treatment', 'max_size', 'min_size', 'response', mod_feat) ]
# model_data_full <- nov_mod_data[,c('id', 'Treatment', 'max_size', 'min_size', 'response', mod_feat) ]
# 
# 
# ##########
# # for presentation
# ##########
library(reshape2)
nov_drugs <- mod_norm[, 1:63]
nov_drugs_melt <- melt(nov_drugs, id.vars = 'id')
# recode valu
nov_drugs_melt$value <-  ifelse(nov_drugs_melt$value == 'CR', 'yes', 
                                ifelse(nov_drugs_melt$value == 'PR', 'yes', 'no'))
nov_drugs_melt <- nov_drugs_melt[complete.cases(nov_drugs_melt),]

# group by id and drug and get counts

drug_counts <- nov_drugs_melt %>%
  group_by(variable) %>%
  summarise(counts = n())

drug_counts_res <- nov_drugs_melt %>%
  group_by(variable) %>%
  summarise(response = sum(value == 'yes'),
            no_response = sum(value == 'no'))


##########
# read in mRECIST for classification
##########
nov_outcome <- read.csv(paste0(data_folder, 'pdxe_outcome.csv'))

# transpose 
nov_outcome <- as.data.frame(t(nov_outcome), stringsAsFactors = F)

# get colnames
colnames(nov_outcome) <- nov_outcome[1,]
nov_outcome <- nov_outcome[-1,]
nov_feat <- colnames(nov_outcome)
nov_outcome$id <- rownames(nov_outcome)
nov_outcome <- nov_outcome[, c('id', nov_feat)]

# join with novartis normal and novartis full
mod_norm <- inner_join(nov_outcome, novartis, by = 'id')
mod_full <- inner_join(nov_outcome, novartis_full,  by = 'id')

# list of drugs 
drug_list <- colnames(mod_norm)[ 2:63]

##########
# predict on novartis classifier
##########
# CR - complete response
# PR - partial response,
# PD - progressive disease
# SD - stable disease
# Random forest model 
rfPredictFac <- function(model_data,
                         drug_name,
                         cutoff,
                         iterations) 
{
  
  model <- list()
  best_features <- list()
  importance <- list()
  train.predictions <- list()
  test.predictions <- list()
  train.ground_truth <- list()
  test.ground_truth <- list()
  test_acc <- list()
  test_stats  <- list()
  
  selected_features <- names(model_data)[65:ncol(model_data)]
  
  model_data <- model_data[, c('id', drug_name, selected_features )]
  model_data <- model_data[complete.cases(model_data),]
  model_data[, drug_name] <- ifelse(model_data[, drug_name] == 'CR', 'yes', 
                                    ifelse(model_data[, drug_name] == 'PR', 'yes', 'no'))
  
  # balance class
  neg <- length(which(model_data[, drug_name] == 'no'))
  pos <- length(which(model_data[, drug_name] == 'yes'))
  stopifnot(pos >10)
  great <- pos > neg
  if(great == F){
    diff <- neg - pos
    remove_index <- which(model_data[, drug_name] == 'no') 
    amount_yes <- length(which(model_data[, drug_name] == 'yes'))
    amount_no <- length(which(model_data[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_no - amount_yes) - 10, replace = F )
    model_data <- model_data[-remove_index,]
  } else {
    diff <- pos - neg
    remove_index <- which(model_data[, drug_name] == 'yes') 
    amount_yes <- length(which(model_data[, drug_name] == 'yes'))
    amount_no <- length(which(model_data[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_yes - amount_no)+5, replace = F )
    model_data <- model_data[-remove_index,]
    
  }
  
  new_neg <- length(which(model_data[, drug_name] == 'no'))
  new_pos <- length(which(model_data[, drug_name] == 'yes'))

  dims <- dim(model_data)
  
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    # determines how you train the model.
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),
      classProbs = TRUE,
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = twoClassSummary
      
    )
    
    y = as.factor(model_data[,drug_name][train_index])
  
    mtry <- sqrt(ncol(model_data[train_index, selected_features]))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = model_data[train_index, selected_features]
                        , y = y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model[[i]])[[1]]
    importance[[i]] <- cbind(rownames(temp), temp$X1)
    
    test.predictions[[i]] <- predict(model[[i]] 
                                     , newdata = model_data[-train_index, selected_features]
                                     , type = "prob")
    
    train.predictions[[i]] <- predict(model[[i]] 
                                      , newdata = model_data[train_index, selected_features]
                                      ,type = "prob")
    
    
    
    
    
    train.ground_truth[[i]] <- as.factor(make.names(model_data[, drug_name][train_index]))
    test.ground_truth[[i]] <- as.factor(make.names(model_data[, drug_name][-train_index]))
    
    # For age of diagnosis
    # Accuracy
    test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(test.predictions[[i]])[1]
    # Confustion Matrix
    test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    
    
    print(i)
  }
  
  return(list(train.predictions, test.predictions, train.ground_truth, test.ground_truth,
              test_acc, test_stats, model, importance, dims, new_pos, new_neg))
  
  
}

##########################
drug <- print(drug_list[[2]])
# "binimetinib"
bini_norm <- rfPredictFac(mod_norm, 
                          drug, 
                          cutoff = 0.7, 
                          iterations = 10)

bini_tab_norm <- print(conMatrix(bini_norm)[[1]])
bini_acc_norm <- print(conMatrix(bini_norm)[[2]])

bini_norm[[9]] #dims
bini_norm[[10]] #pos
bini_norm[[11]] #neg
 
bini_full <- rfPredictFac(mod_full, 
                          drug, 
                          cutoff = 0.7, 
                          iterations = 10)
bini_tab_full <- print(conMatrix(bini_full)[[1]])
bini_acc_full <- print(conMatrix(bini_full)[[2]])

bini_full[[9]] # dims
bini_full[[10]] # pos
bini_full[[11]]# neg
#this one
##########################
drug <- print(drug_list[[3]])
# "BYL719_LJM716"
BYL719_LJM716_norm <- rfPredictFac(mod_norm, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
BYL719_LJM716_tab_norm <- print(conMatrix(BYL719_LJM716_norm)[[1]])
BYL719_LJM716_acc_norm <- print(conMatrix(BYL719_LJM716_norm)[[2]])

BYL719_LJM716_norm[[9]] # dims 
BYL719_LJM716_norm[[10]] # pos
BYL719_LJM716_norm[[11]] # neg

BYL719_LJM716_full <- rfPredictFac(mod_full, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
BYL719_LJM716_tab_full <- print(conMatrix(BYL719_LJM716_full)[[1]])
BYL719_LJM716_acc_full <- print(conMatrix(BYL719_LJM716_full)[[2]])

BYL719_LJM716_full[[9]] # dims
BYL719_LJM716_full[[10]] # pos
BYL719_LJM716_full[[11]] # neg

##########################
drug <- print(drug_list[[6]])
# "BYL719_+_LJM716"
BYL719_LJM716_norm <- rfPredictFac(mod_norm, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
BYL719_LJM716_tab_norm <- print(conMatrix(BYL719_LJM716_norm)[[1]])
BYL719_LJM716_acc_norm <- print(conMatrix(BYL719_LJM716_norm)[[2]])

BYL719_LJM716_norm[[9]] # dims 
BYL719_LJM716_norm[[10]] # pos
BYL719_LJM716_norm[[11]] # neg

BYL719_LJM716_full <- rfPredictFac(mod_full, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
BYL719_LJM716_tab_full <- print(conMatrix(BYL719_LJM716_full)[[1]])
BYL719_LJM716_acc_full <- print(conMatrix(BYL719_LJM716_full)[[2]])

BYL719_LJM716_full[[9]] # dims
BYL719_LJM716_full[[10]] # pos
BYL719_LJM716_full[[11]] # neg

##########################
drug <- print(drug_list[[7]])
# "CLR457"
CLR457_norm <- rfPredictFac(mod_norm, 
                                   drug, 
                                   cutoff = 0.7, 
                                   iterations = 10)
CLR457_tab_norm <- print(conMatrix(CLR457_norm)[[1]])
CLR457_acc_norm <- print(conMatrix(CLR457_norm)[[2]])

CLR457_norm[[9]] # dims 
CLR457_norm[[10]] # pos
CLR457_norm[[11]] # neg

CLR457_full <- rfPredictFac(mod_full, 
                                   drug, 
                                   cutoff = 0.7, 
                                   iterations = 10)
CLR457_tab_full <- print(conMatrix(CLR457_full)[[1]])
CLR457_acc_full <- print(conMatrix(CLR457_full)[[2]])

CLR457_full[[9]] # dims
CLR457_full[[10]] # pos
CLR457_full[[11]] # neg


##########################
drug <- print(drug_list[[14]])
# "LEE011"
LEE011_norm <- rfPredictFac(mod_norm, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
LEE011_tab_norm <- print(conMatrix(LEE011_norm)[[1]])
LEE011_acc_norm <- print(conMatrix(LEE011_norm)[[2]])

LEE011_norm[[9]] # dims 
LEE011_norm[[10]] # pos
LEE011_norm[[11]] # neg

LEE011_full <- rfPredictFac(mod_full, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
LEE011_tab_full <- print(conMatrix(LEE011_full)[[1]])
LEE011_acc_full <- orint(conMatrix(LEE011_full)[[2]])

LEE011_full[[9]] # dims
LEE011_full[[10]] # pos
LEE011_full[[11]] # neg


##########################
drug <- print(drug_list[[15]])
# "LEE011_everolimus"
LEE011_everolimus_norm <- rfPredictFac(mod_norm, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
LEE011_everolimus_tab_norm <- print(conMatrix(LEE011_everolimus_norm)[[1]])
LEE011_everolimus_acc_norm <- print(conMatrix(LEE011_everolimus_norm)[[2]])

LEE011_everolimus_norm[[9]] # dims 
LEE011_everolimus_norm[[10]] # pos
LEE011_everolimus_norm[[11]] # neg

LEE011_everolimus_full <- rfPredictFac(mod_full, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
LEE011_everolimus_tab_full <- print(conMatrix(LEE011_everolimus_full)[[1]])
LEE011_everolimus_acc_full <- print(conMatrix(LEE011_everolimus_full)[[2]])

LEE011_everolimus_full[[9]] # dims
LEE011_everolimus_full[[10]] # pos
LEE011_everolimus_full[[11]] # neg
#this one

##########################
drug <- print(drug_list[[23]])
# "BYL719_binimetinib"
BYL719_binimetinib_norm <- rfPredictFac(mod_norm, 
                                       drug, 
                                       cutoff = 0.7, 
                                       iterations = 10)
BYL719_binimetinib_tab_norm <- print(conMatrix(BYL719_binimetinib_norm)[[1]])
BYL719_binimetinib_acc_norm <- print(conMatrix(BYL719_binimetinib_norm)[[2]])

BYL719_binimetinib_norm[[9]] # dims 
BYL719_binimetinib_norm[[10]] # pos
BYL719_binimetinib_norm[[11]] # neg

BYL719_binimetinib_full <- rfPredictFac(mod_full, 
                                       drug, 
                                       cutoff = 0.7, 
                                       iterations = 10)
BYL719_binimetinib_tab_full <- print(conMatrix(BYL719_binimetinib_full)[[1]])
BYL719_binimetinib_acc_full <- print(conMatrix(BYL719_binimetinib_full)[[2]])

BYL719_binimetinib_full[[9]] # dims
BYL719_binimetinib_full[[10]] # pos
BYL719_binimetinib_full[[11]] # neg
# this one

#########################
drug <- print(drug_list[[37]])
# "BKM120_binimetinib"
BKM120_binimetinib_norm <- rfPredictFac(mod_norm, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
BKM120_binimetinib_tab_norm <- print(conMatrix(BKM120_binimetinib_norm)[[1]])
BKM120_binimetinib_acc_norm <- print(conMatrix(BKM120_binimetinib_norm)[[2]])

BKM120_binimetinib_norm[[9]] # dims 
BKM120_binimetinib_norm[[10]] # pos
BKM120_binimetinib_norm[[11]] # neg

BKM120_binimetinib_full <- rfPredictFac(mod_full, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
BKM120_binimetinib_tab_full <- print(conMatrix(BKM120_binimetinib_full)[[1]])
BKM120_binimetinib_acc_full <- print(conMatrix(BKM120_binimetinib_full)[[2]])

BKM120_binimetinib_full[[9]] # dims
BKM120_binimetinib_full[[10]] # pos
BKM120_binimetinib_full[[11]] # neg

#########################
drug <- print(drug_list[[40]])
# gemcitabine-50mpk
gemcitabine_50mpk_norm <- rfPredictFac(mod_norm, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
gemcitabine_50mpk_tab_norm <- print(conMatrix(gemcitabine_50mpk_norm)[[1]])
gemcitabine_50mpk_acc_norm <- print(conMatrix(gemcitabine_50mpk_norm)[[2]])

gemcitabine_50mpk_norm[[9]] # dims 
gemcitabine_50mpk_norm[[10]] # pos
gemcitabine_50mpk_norm[[11]] # neg

gemcitabine_50mpk_full <- rfPredictFac(mod_full, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
gemcitabine_50mpk_tab_full <- print(conMatrix(gemcitabine_50mpk_full)[[1]])
gemcitabine_50mpk_acc_full <- print(conMatrix(gemcitabine_50mpk_full)[[2]])

gemcitabine_50mpk_full[[9]] # dims
gemcitabine_50mpk_full[[10]] # pos
gemcitabine_50mpk_full[[11]] # neg

#########
# predict on ming pdx
##########
# Random forest model 
rfPredictTest <- function(model_data,
                         drug_name) 
{
  
  dims <- dim(model_data)
  novartis <- model_data[grepl('X.', model_data$id),]
  pdxe <-model_data[!grepl('X.', model_data$id),]
  mod_norm <- inner_join(nov_outcome, novartis, by = 'id')
  
  selected_features <- names(mod_norm)[65:ncol(mod_norm)]
  
  mod_norm <- mod_norm[, c('id', drug_name, selected_features )]
  mod_norm <- mod_norm[complete.cases(mod_norm),]
  mod_norm[, drug_name] <- ifelse(mod_norm[, drug_name] == 'CR', 'yes', 
                                    ifelse(mod_norm[, drug_name] == 'PR', 'yes', 'no'))
  
  # balance class
  neg <- length(which(mod_norm[, drug_name] == 'no'))
  pos <- length(which(mod_norm[, drug_name] == 'yes'))
  stopifnot(pos >10)
  great <- pos > neg
  if(great == F){
    diff <- neg - pos
    remove_index <- which(mod_norm[, drug_name] == 'no') 
    amount_yes <- length(which(mod_norm[, drug_name] == 'yes'))
    amount_no <- length(which(mod_norm[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_no - amount_yes) - 10, replace = F )
    mod_norm <- mod_norm[-remove_index,]
  } else {
    diff <- pos - neg
    remove_index <- which(mod_norm[, drug_name] == 'yes') 
    amount_yes <- length(which(mod_norm[, drug_name] == 'yes'))
    amount_no <- length(which(mod_norm[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_yes - amount_no)+5, replace = F )
    mod_norm <- mod_norm[-remove_index,]
    
  }
  
  new_neg <- length(which(mod_norm[, drug_name] == 'no'))
  new_pos <- length(which(mod_norm[, drug_name] == 'yes'))
  
  dims <- dim(mod_norm)

  # determines how you train the model.
  NFOLDS <- 2
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary
    
  )
  
  y = as.factor(mod_norm[,drug_name])
  
  mtry <- sqrt(ncol(mod_norm))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = mod_norm[, selected_features]
                      , y = y
                      , method = "rf"
                      , trControl = fitControl
                      , tuneGrid = tunegrid
                      , importance = T
                      , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$X1)
  
  test.predictions <- predict(model 
                                   , newdata = pdxe[, selected_features]
                                   , type = "prob")
  

  
  return(list(test.predictions, model, importance, dims, new_pos, new_neg))
  
}

#2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[2]])
# "binimetinib"

binimetinib_norm <- rfPredictTest(dat_norm_cor, 
                      drug)[[1]]
binimetinib_full <- rfPredictTest(dat_full_cor, 
                      drug)[[1]]

result <- binimetinib_norm[order(binimetinib_norm$yes), ]


#2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[3]])
# BKM120

BKM120_norm <- rfPredictTest(dat_norm_cor, 
                                  drug)
BKM120_full <- rfPredictTest(dat_full_cor, 
                                  drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[6]])
# "BYL719_LJM716"

BYL719_LJM716_norm <- rfPredictTest(dat_norm_cor, 
                             drug)
BYL719_LJM716_full <- rfPredictTest(dat_full_cor, 
                             drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[7]])
# "CLR457"

CLR457_norm <- rfPredictTest(dat_norm_cor, 
                                    drug)
CLR457_full <- rfPredictTest(dat_full_cor, 
                                    drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[14])
# "LEE011"

LEE011_norm <- rfPredictTest(dat_norm_cor, 
                             drug)
LEE011_full <- rfPredictTest(dat_full_cor, 
                             drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[23])
# "BYL719_binimetinib"

BYL719_binimetinib_norm <- rfPredictTest(dat_norm_cor, 
                             drug)
LBYL719_binimetinib_full <- rfPredictTest(dat_full_cor, 
                             drug)

# 2,3,6,7,14,23,37,40
############################
drug <- print(drug_list[37])
# "BKM120_binimetinib"

BKM120_binimetinib_norm <- rfPredictTest(dat_norm_cor, 
                                         drug)
BKM120_binimetinib_full <- rfPredictTest(dat_full_cor, 
                                          drug)

# 2,3,6,7,14,23,37,40
############################
drug <- print(drug_list[40])
# "gemcitabine-50mpk"

gemcitabine_50mpk_norm <- rfPredictTest(dat_norm_cor, 
                                         drug)
gemcitabine_50mpk_full <- rfPredictTest(dat_full_cor, 
                                         drug)


