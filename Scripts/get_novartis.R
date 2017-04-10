##########################
# this is the 1st step in the pipeline - read in novartis data and clean


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
project_folder <- paste0(home_folder, 'pdxSNF/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in novartis data
##########
rna_seq <- read.csv(paste0(data_folder, '/novartis', '/rna_seq.csv'))

# make row one the row names then remove 
rownames(rna_seq) <- rna_seq$Sample
rna_seq$Sample <- NULL

# transpose data 
rna_seq <- t(rna_seq)

##########
# clean data
##########

# remove zero variance genes
near_zero <- nearZeroVar(rna_seq, saveMetrics = T)

keep_these <- rownames(near_zero)[near_zero$nzv != T]

rna_seq <- rna_seq[, keep_these]

# remove highly correlated genes
x_mat <- cor(as.matrix(rna_seq))

highlyCorDescr <- findCorrelation(x_mat, cutoff = .80)

dim(rna_seq)
rna_seq <- rna_seq[, -highlyCorDescr]
dim(rna_seq)


# remove old objects
rm(near_zero, x_mat, keep_these)

# ##########
# # normalize 
# ##########
# rna_seq_2 <- scale(rna_seq, scale = F)
# rna_seq <- scale(rna_seq)


##########
# log2 and quantile 
##########

# first replace zeros with 1s
rna_seq[rna_seq == 0] <- 1

# apply log2 transformation
rna_seq <- log2(rna_seq)

# get row and col names
row_names <- rownames(rna_seq)
col_names <- colnames(rna_seq)


# apply qunatile normilzation
rna_norm <- normalize.quantiles(rna_seq)
rownames(rna_norm) <- row_names
colnames(rna_norm) <- col_names


# save data 
saveRDS(rna_norm, paste0(data_folder, '/rna_normalized.rda'))


#######################



