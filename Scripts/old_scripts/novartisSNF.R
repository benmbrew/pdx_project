##########################
# this script will run SNF on the data

##########
# initiate libraries
##########
library(dplyr)
library(SNFtool)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdxSNF/')
data_folder <- paste0(project_folder, 'Data')
code_folder <- paste0(project_folder, 'Code')

##########
# source functions
##########
source(paste0(code_folder, '/functions.R'))

##########
# read in data and transpose it
##########
rna_seq <- read.csv(paste0(data_folder, '/novartis', '/rna_seq.csv'))
copy_number <- read.csv(paste0(data_folder, '/novartis', '/copy_number.csv'))

# make row one the row names then remove 
rownames(rna_seq) <- rna_seq$Sample
rna_seq$Sample <- NULL

# if duplicated paste a .x
rownames_for_cn <- make.names(copy_number$Sample, unique = T)
rownames(copy_number) <- rownames_for_cn
copy_number$Sample <- NULL

##########
# combine data sets and take column intersection
##########
cases <- list()
cases[[1]] <- rna_seq
cases[[2]] <- copy_number


##########
# get column intersection and sorted ids
##########
cases <- columnIntersection(cases)[[1]]
sorted_ids <- columnIntersection(cases)[[2]]

##########
# normalize data
##########
stat <- rowStatistics(cases)
cases <- normalizeData(cases, stat)
temp <- cases[[1]]
##########
# run SNF
##########

# samples ids are in columns
sampleRows <- FALSE

# run snf
labels <- SNFClustering(cases, 5, sampleRows)

# get id for labels
labels <- as.data.frame(cbind(labels, id = sorted_ids))

# save labels
write.csv(labels, (paste0(data_folder, '/novartis_labels.csv')))
