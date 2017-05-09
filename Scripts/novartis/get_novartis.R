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
project_folder <- paste0(home_folder, 'pdx_project/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in novartis data
##########
rna_seq <- read.csv(paste0(data_folder, '/novartis', '/rna_seq.csv'), stringsAsFactors = F)

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
dim(x_mat)

highlyCorDescr <- findCorrelation(x_mat, cutoff = .80)

dim(rna_seq)
rna_seq <- rna_seq[, -highlyCorDescr]
dim(rna_seq)


# remove old objects
rm(near_zero, x_mat, keep_these)

###########
# subset to lung 
###########

# read in pdxe info
nov_map <- read.csv(paste0(data_folder, '/novartis', '/pdxe_model.csv'), stringsAsFactors = F)

# get summary of cancer types 
summary(as.factor(nov_map$tumor.type))

# subet to NSCLC
nov_map <- nov_map[nov_map$tumor.type == 'NSCLC',]

# check NAs
length(which(is.na(nov_map)))

##########
# join nov_map and rna seq to grab only lung cancer
##########
# make rna_seq a data frame
rna_seq <- as.data.frame(rna_seq, strinsAsFactors = F)
# first get id colmn for rna_seq
rna_seq$patient.id <- rownames(rna_seq)

# check NAs
length(which(is.na(rna_seq)))

# replace '-' with '.' in patient id from nov map
nov_map$patient.id <- gsub('-', '.', nov_map$patient.id, fixed = T)

# check uniquness of joining variable on both data sets
length(unique(nov_map$patient.id))
length(unique(rna_seq$patient.id))

# join data 
nov_dat <- inner_join(nov_map, rna_seq, by = 'patient.id')

# check NAs
length(which(is.na(nov_dat)))
length(which(is.na(rna_seq)))


# get number of unique models
length(unique(nov_dat$patient.id))

# remove duplicates
nov_dat <- nov_dat[!duplicated(nov_dat$patient.id),]

# ##########
# scale data
# ##########

# first get clin dat 
clin_dat <- nov_dat[, 1:6]

# now remove clin part of nov_dat and convert to matrix for scale
nov_dat <- as.matrix(nov_dat[, -c(1:6)])

# scale nov_dat
nov_dat <- apply(nov_dat, 2, function(y) (y - mean(y)) / sd(y) ^ as.logical(sd(y)))

# cbind nov_dat and clin_dat
nov_dat <- cbind(clin_dat, nov_dat)


# ##########
# # log2 and quantile 
# ##########
# 
# # first replace zeros with 1s
# rna_seq[rna_seq == 0] <- 1
# 
# # apply log2 transformation
# rna_seq <- log2(rna_seq)
# 
# # get row and col names
# row_names <- rownames(rna_seq)
# col_names <- colnames(rna_seq)
# 
# 
# # apply qunatile normilzation
# rna_norm <- normalize.quantiles(rna_seq)
# rownames(rna_norm) <- row_names
# colnames(rna_norm) <- col_names
# 

# save data 
saveRDS(nov_dat, paste0(data_folder, '/nov_dat_scaled.rda'))


#######################



