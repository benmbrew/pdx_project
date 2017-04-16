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
rna <- read.csv(paste0(champ_onc_folder, '/rna_seq_expression.csv'), stringsAsFactors = F)

# cnv_export
cnv <- read.csv(paste0(champ_onc_folder, '/cnv_export.csv'), stringsAsFactors = F)

# drugs
drugs <- read.csv(paste0(champ_onc_folder, '/drugs.csv'), stringsAsFactors = F)

# fusions
fusions <- read.csv(paste0(champ_onc_folder, '/fusions.csv'), stringsAsFactors = F)

# mutations_export
mut <- read.csv(paste0(champ_onc_folder, '/mutations_export.csv'), stringsAsFactors = F)

# treatment_history
treatment <- read.csv(paste0(champ_onc_folder, '/treatment_history.csv'), stringsAsFactors = F)

# readin model overview
mod_summary <- read.csv(paste0(champ_onc_folder, '/mod_overview.csv'), stringsAsFactors = F)

# read in s3 file info
s3_files <- read.csv(paste0(champ_onc_folder, '/s3_files.csv'), stringsAsFactors = F)

##########
# explore mod_summary and s3_files 
##########

length(unique(mod_summary$Model))
# 971 unique models 

mod_summary %>% group_by(Tumor.status) %>% summarise(counts = n())
# 4 cat for met- local, primary, metastic, NA

mod_summary %>% group_by(Cancer.type) %>% summarise(counts = n())
# 51 unique cancers

mod_summary %>% group_by(Model.status) %>% summarise(counts = n())
# 4 cat: early banking, engrafter, established, late banking.

mod_summary %>% group_by(Histology) %>% summarise(counts = n())
# histology is the tissue (study of microscopic tissues)

mod_summary %>% group_by(Harvest.site) %>% summarise(counts = n())
# lots of cats

mod_summary %>% group_by(Disease.stage) %>% summarise(counts = n())
# 5 cats: I, II, III, IV, NA

mod_summary %>% group_by(Tumor.grade) %>% summarise(counts = n())
# 8 cats: well differentiated, moderately diff, moderately well, undifferentiated, etc

mod_summary %>% group_by(Diagnosis) %>% summarise(counts = n())
# 7 cat: relapsed, refractory, recurrent, Na, fist diagnosis, de novo, and blank

mod_summary %>% group_by(Treatment.history) %>% summarise(counts = n())
# 3 cats naive, NA, pretreated

summary(as.numeric(mod_summary$Age))
# mean 55, min 1, max 93, 99 NAs

mod_summary %>% group_by(Gender) %>% summarise(counts = n())
#440 female, 474 male

mod_summary %>% group_by(Ethnicity) %>% summarise(counts = n())
# 593 white

mod_summary %>% group_by(Exome.seq.status) %>% summarise(counts = n())
#no 444, yes 527

mod_summary %>% group_by(RNA.seq.status) %>% summarise(counts = n())
# 449 no, 522 yes

mod_summary %>% group_by(Characterization) %>% summarise(counts = n())
# 522 available, 449 blank (NA)

##########
# look at rna file and get to model by gene format
##########

# rename
names(rna)

colnames(rna) <- c('model_id', 'tumor_type', 'gene', 'gene_id', 'log_rpkm', 'z', 'mean_log_rpkm', 'sd', 'fold', 'trans')

# how many unique models 
length(unique(rna$model_id))
# 515

# how many unqiue genes
length(unique(rna$gene))

# how many z stats greater than + 1.65 or less than -1.65
length(which(rna$z > 1.65))
# 5179
length(which(rna$z < -1.65))
#3401

# .05 sig
rna_sig <- rna[rna$z > 1.65 | rna$z < -1.65,]

# group by model_id and gene and get mean(log_rpkm)
rna_sig_by_gene <- rna_sig %>% group_by(model_id, gene) %>% summarise(mean_log_rpkm = mean(mean_log_rpkm))

##########
# rna with mean log rpkm 
##########

# how many genes per model 
temp <- rna_sig_by_gene %>% group_by(model_id) %>% summarise(counts = n())
mean(temp$counts)
# on avg each model has only 16 genes

# put in wide format(maybe use loop)

