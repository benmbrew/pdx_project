
##########
# initialize libraries
##########
library(caret)
library(glmnet)
library(randomForest)
library(kernlab)
library(pROC)
library(Metrics)
library(doParallel)
library(nnet)
library(dplyr)
library(bumphunter)
library(sqldf)

registerDoParallel(1)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdx_project/')
data_folder <- paste0(project_folder, 'Data/')
champ_onc_folder <- paste0(data_folder, 'champ_onc/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in data
##########
rna_mod <- readRDS(paste0(champ_onc_folder, '/rna_mod.rda'))

length(unique(rna_mod$treatment))
# 259 treatments

# get drug counts 
temp <- rna_mod %>% group_by(treatment) %>% summarise(counts = n())
temp <- temp[order(temp$counts, decreasing = T),]

