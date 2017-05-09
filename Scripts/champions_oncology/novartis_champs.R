###############################################
# find overlapping drugs-tissue pairs between novartis and 

##########
# initialize libraries
##########
library(plyr)
library(dplyr)
library(lumi)
library(reshape2)
library(tidyr)

##########
# initialize folders
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdx_project/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
novartis_folder <- paste0(data_folder, 'novartis/')
champ_onc_folder <- paste0(data_folder, 'champ_onc/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in mRECIST for classification
##########
# read in outcome
nov_outcome <- read.csv(paste0(novartis_folder, 'pdxe_outcome.csv'))

# transpose 
nov_outcome <- as.data.frame(t(nov_outcome), stringsAsFactors = F)

# get colnames
colnames(nov_outcome) <- nov_outcome[1,]
nov_outcome <- nov_outcome[-1,]
nov_feat <- colnames(nov_outcome)
nov_outcome$id <- rownames(nov_outcome)
nov_outcome <- nov_outcome[, c('id', nov_feat)]

nov <- gather(nov_outcome, id)

# rename nov
colnames(nov) <- c('patient.id', 'drug', 'value')
nov %>% group_by(patient.id,drug)
rm(nov_outcome)

# read in novartis pdx map
pdx_nov <- read.csv(paste0(novartis_folder, 'pdxe_model.csv'))

# recode patient.id gsub '-' to '.'
pdx_nov$patient.id <- gsub('-', '.', pdx_nov$patient.id, fixed = T)

##########
# join pdx nov and nov by patient.id and id
##########

# merge patient.id and tumor type from pdx nov
pdx_nov$id_tumor <- paste(pdx_nov$patient.id, pdx_nov$tumor.type.name, sep = '_')

# now get unique id_tumor
temp.1 <- unique(pdx_nov$id_tumor)

# now seperate again and make data frame with unique id and corresponding tissue
temp.2 <- as.data.frame(do.call(rbind, strsplit(temp.1, split = '_')))
names(temp.2) <- c('patient.id', 'tissue')

# now left join nov onto temp.2
drug_tiss_set <- left_join(temp.2, nov, by = 'patient.id')

# now merge drug and tissue
temp.3 <- paste(drug_tiss_set$tissue, drug_tiss_set$drug, sep= '_')

temp.4 <- unique(temp.3)

temp.5 <- as.data.frame(do.call(rbind, strsplit(temp.4, split = '_')))

names(temp.5) <- c('tissue', 'drug')

full_nov <- left_join(temp.5, drug_tiss_set)
full_nov$value <- NULL

##########
# read in champs onclology
##########

# treatment_history
treatment <- read.csv(paste0(champ_onc_folder, '/treatment_history.csv'), stringsAsFactors = F)

# readin model overview
mod_summary <- read.csv(paste0(champ_onc_folder, '/mod_overview.csv'), stringsAsFactors = F)

# left join treatment and mod_summary by model id
# mod summary has 971 uniue models and corresponding cancer
# treatment has all the drugs for every model. 

champs_dat <- left_join(treatment, mod_summary, by= c('model_id' = 'Model'))
names(champs_dat)
# combine model_id and Cancer.type
champs_dat$new_var <- paste(champs_dat$treatment, champs_dat$Cancer.type, sep = '_')

# get unique
temp_1 <- unique(champs_dat$new_var)

# now seperate again and make data frame with unique id and corresponding tissue
temp_2 <- as.data.frame(do.call(rbind, strsplit(temp_1, split = '_')))
names(temp_2) <- c('treatment', 'Cancer.type')

# now left join champs_dat into temp2
full_champs <- left_join(temp_2, champs_dat)

# keep only treatment, cancer.type, and id
full_champs <- full_champs[, c('Cancer.type','treatment',  'model_id')]
names(full_champs) <- c('tissue', 'drug', 'patient.id')

##########
# full_nov and full_champs
##########

# look at tissue columns and homogenize
summary(full_nov$tissue)
# Breast Cancer, Colorectal Cancer, Cutaneous Melanoma, Gastric Cancer, 
# Non-small Cell Lung Carcinoma, Pancreatic Ductal Carcinoma 

any(grepl('lung', summary(full_champs$tissue)))
# Breast, Colorectal, Melanoma, Gastric, Pancreatic,

# recode full_nov 
full_nov$tissue <- ifelse(grepl('Breast', full_nov$tissue), 'Breast',
                 ifelse(grepl('Colorectal', full_nov$tissue), 'Colorectal',  
                        ifelse(grepl('Melanoma', full_nov$tissue), 'Melanoma',
                               ifelse(grepl('Gastric', full_nov$tissue), 'Gastric', 
                                      ifelse(grepl('Pancreatic', full_nov$tissue), 'Pancreatic', 'Lung')))))

# now tissue is homogneized
# homogenize drugs
summary(as.factor(full_nov$drug))
summary(as.factor(full_champs$drug))

# gsub '+' for '/'
full_nov$drug <- gsub(' + ', '/', full_nov$drug, fixed = T)

# make both entirely lower case and seperated by 
full_nov$drug <- tolower(full_nov$drug)
full_champs$drug <- tolower(full_champs$drug)
# use these later

short_nov <- full_nov %>% group_by(tissue, drug) %>% summarise(counts = n())
short_champs <- full_champs %>% group_by(tissue, drug) %>% summarise(counts = n())

# sort by drug
short_nov <- short_nov[order(short_nov$drug),]
short_champs <- short_champs[order(short_champs$drug),]

# recode short_champs$drug 
short_champs$drug <- gsub('gemcitabine', 'gemcitabine-50mpk', short_champs$drug)
short_champs$drug <- gsub('5-fluorouracil', '5fu', short_champs$drug)

# loop through tissues and match on that, much smaller set (full_nov and full_champs)
intersecting_tissue <- intersect(short_nov$tissue, short_champs$tissue)
int_drugs <- list()

for (i in intersecting_tissue) {
  
  #subet both by cancer
  sub_nov <- short_nov[short_nov$tissue == i,]
  sub_champ <- short_champs[short_champs$tissue == i,]
  
  if(i != 'Melanoma') {
    # join by drug
    sub_dat <- inner_join(sub_nov, sub_champ, by = 'drug')

    int_drugs[[i]] <- as.data.frame(sub_dat)
  }
  
 

  
  print(i)
}

temp_results <- do.call(rbind, int_drugs)

# remove columns
results <- temp_results[, c("tissue.x", "drug", "counts.x", "counts.y")]
names(results) <- c('Tissue', 'Drug', 'nov_counts', 'champs_counts')

write.csv(results, '/home/benbrew/Desktop/results.csv')
