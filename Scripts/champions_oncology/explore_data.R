##################################################################################
# this script will explore champions oncology data

##########
# initialize libraries
##########
library(dplyr)
library(lumi)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdx_project/')
data_folder <- paste0(project_folder, 'Data/')
champ_onc_folder <- paste0(data_folder, 'champ_onc/')
code_folder <- paste0(project_folder, 'Code/')

##########
# Read in data
##########

# rna_seq_expression
rna <- read.csv(paste0(champ_onc_folder, '/rna_seq_expression.csv'))

# cnv_export
cnv <- read.csv(paste0(champ_onc_folder, '/cnv_export.csv'))

# drugs
drugs <- read.csv(paste0(champ_onc_folder, '/drugs.csv'))

# fusions
fusions <- read.csv(paste0(champ_onc_folder, '/fusions.csv'))

# mutations_export
mut <- read.csv(paste0(champ_onc_folder, '/mutations_export.csv'))

# treatment_history
treatment <- read.csv(paste0(champ_onc_folder, '/treatment_history.csv'))

