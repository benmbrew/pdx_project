##################################################################################
# this script will explore champions oncology data

##########
# initialize libraries
##########
library(plyr)
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
# explore mod_summary
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

# recode mod_summary colnames
colnames(mod_summary) <- tolower(colnames(mod_summary))
colnames(mod_summary) <- gsub('.', '_', colnames(mod_summary), fixed = T)
colnames(mod_summary)[1] <- 'model_id'

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

model <- 'CTG-0011'
# put in wide format(maybe use loop)
result_list <- list()
for(model in unique(rna_sig_by_gene$model_id)){
  
  sub_dat <- as.data.frame(t(rna_sig_by_gene[rna_sig_by_gene$model_id == model,]), stringsAsFactors = F)
  sub_dat$model_id <- unique(sub_dat[1,1])
  colnames(sub_dat) <- sub_dat[2,] 
  result_list[[model]] <- sub_dat[-c(1:2),]
  print(model)
  
}

# what is the minimun
temp <- lapply(result_list, function(x) unlist(colnames(x)))

# histogram of gene counts for each model. 
hist(unlist(lapply(temp, function(x) length(x))), breaks = 20)

# how many models have over 20 genes
length(which(unlist(lapply(temp, function(x) length(x))) > 5))

# # create column in each sub data frame that indicates how many genes we have for that model
# length(result_list) 
# for(model in 1:length(result_list)) {
#   sub_dat <- result_list[[model]]
#   sub_dat$num_genes <- (length(colnames(sub_dat)) - 1)
#   result_list[[model]] <- sub_dat
#   print(model)
# }

##########
# collapse data
##########

# first make last column for model_id in each data frame 
for(i in 1:length(result_list)) {
  sub_dat <- result_list[[i]]
  colnames(sub_dat)[length(sub_dat)] <- 'model_id'
  feats <- colnames(sub_dat)[(1:length(colnames(sub_dat))- 1)]
  sub_dat <- sub_dat[, c('model_id', feats)]
  result_list[[i]] <- sub_dat
  print(i)
}

rna_data <- rbind.fill(result_list)

###########
# we now have a data frame with 511 unique models and 195 unique genes
###########

# avg amount of missing values per gene 
mean(apply(rna_data[, 2:ncol(rna_data)], 2,function(x) length(which(!is.na(x)))))

# # get gene names and corresponding counts 
# temp_counts <- as.data.frame(apply(rna_data[, 2:ncol(rna_data)], 2,function(x) length(which(!is.na(x)))))
# colnames(temp_counts) <- 'counts'
# temp_counts$gene <- rownames(temp_counts)
# temp_counts <- temp_counts[order(temp_counts$counts, decreasing = T),]

##########
# merge mod_summary with rna_data, and then merge that with treatment (has outcome)
##########
rna_mod <- left_join(treatment, rna_data, by = c('model_id' = 'model_id'))

# subset to "Responded" and "No response"
rna_mod <- rna_mod[grepl('No response|Responded', rna_mod$outcome),]

###########
# save rna_mod
###########
saveRDS(rna_mod, paste0(champ_onc_folder, '/rna_mod.rda'))


