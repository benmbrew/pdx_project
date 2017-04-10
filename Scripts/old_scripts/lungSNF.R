##########################
# this script will read in lung pdx data and run SNF

##########
# initiate libraries
##########
library(dplyr)
library(lumi)
library(SNFtool)

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
# read in data
##########
quan = readRDS(paste0(data_folder, 'quan_new.rda'))
rsn = readRDS(paste0(data_folder, 'rsn_new.rda'))
ssn = readRDS(paste0(data_folder, 'ssn_new.rda'))
loess = readRDS(paste0(data_folder, 'loess_new.rda'))
rank = readRDS(paste0(data_folder, 'rank_new.rda'))

##########
# get sorted id
##########
quan.sorted_id = sort(quan$model_id)
rsn.sorted_id = sort(rsn$model_id)
ssn.sorted_id = sort(ssn$model_id)
loess.sorted_id = sort(loess$model_id)
rank.sorted_id = sort(rank$model_id)

##########
# put data in SNF format
##########
quan = prepareSNF(quan)
rsn = prepareSNF(rsn)
ssn = prepareSNF(ssn)
loess = prepareSNF(loess)
rank = prepareSNF(rank)

##########
# Run SNF
##########

# quan
quan_labels = lungSNf(data = quan, 
                       numClust = 5, 
                       sampleRows = F, 
                       sorted_ids = quan.sorted_id)
# rsn
rsn_labels = lungSNf(data = rsn, 
                      numClust = 5, 
                      sampleRows = F, 
                      sorted_ids = rsn.sorted_id)

# ssn
ssn_labels = lungSNf(data = ssn, 
                      numClust = 5, 
                      sampleRows = F, 
                      sorted_ids = ssn.sorted_id)

# loess
loess_labels = lungSNf(data = loess, 
                        numClust = 5, 
                        sampleRows = F, 
                        sorted_ids = loess.sorted_id)

# rank
rank_labels = lungSNf(data = rank, 
                       numClust = 5, 
                       sampleRows = F, 
                       sorted_ids = rank.sorted_id)


# save labels
write.csv(quan_labels, (paste0(data_folder, '/quan_labels.csv')))
write.csv(rsn_labels, (paste0(data_folder, '/rsn_labels.csv')))
write.csv(ssn_labels, (paste0(data_folder, '/ssn_labels.csv')))
write.csv(loess_labels, (paste0(data_folder, '/loess_labels.csv')))
write.csv(rank_labels, (paste0(data_folder, '/rank_labels.csv')))
