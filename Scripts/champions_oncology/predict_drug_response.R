
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
library(tidyr)
library(reshape2)

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
# 259 unique drug/drug combinations
hist(temp$counts)
# most unique drug/drug combinations have only less than 5 observations

##########
# look at number of drugs used and outcome variable
##########
treat_outcome <- rna_mod[, c('treatment', 'outcome')]

# strsplit 
treat_outcome$num_drugs_used <- unlist(lapply(strsplit(treat_outcome$treatment, '/'), function(x) length(x)))
str(treat_outcome)
hist(treat_outcome$num_drugs_used)
summary(as.factor(treat_outcome$outcome))

# get rid of treatment variable
treat_out <- treat_outcome[, c('outcome', 'num_drugs_used')]

# group by outcome and get mean
mean_drugs <- treat_out %>% group_by(outcome) %>% summarise(mean_drugs = mean(num_drugs_used))

boxplot(num_drugs_used ~ outcome, data = treat_out)


treat_out$outcome <- as.factor(treat_out$outcome)
glm_mod <- glm(outcome ~ num_drugs_used, data = treat_out, family=binomial(link="logit"))
summary(glm_mod)
coef(glm_mod)
exp(cbind(Odds=coef(glm_mod), confint(glm_mod)))

# Odds     2.5 %    97.5 %
#   (Intercept)    0.6336092 0.4390488 0.9094927
# num_drugs_used 1.6138151 1.3540473 1.9349707

##########
# which drugs are associated with good outcomes (point system)
##########

# splitrows (get other function)
splitRows <- function(data, duplicate_table){
  
  for (i in 1:nrow(data)) {
    sub_data <- data[i,]
    
    if (grepl('/', sub_data$blood_dna_malkin_lab_) | 
        grepl('/', sub_data$age_sample_collection)) {
      
      split_malkin <- strsplit(as.character(sub_data$blood_dna_malkin_lab_), '/')
      split_age <- strsplit(as.character(sub_data$age_sample_collection), '/')
      split_malkin <- cbind(unlist(split_malkin))
      split_age <- cbind(unlist(split_age))
      duplicate <- cbind(split_malkin, split_age, sub_data)
      duplicate_table <- rbind(duplicate_table, duplicate)
      
    }
  }
  
  data <- data[!grepl('/', data$blood_dna_malkin_lab_),]
  data <- data[!grepl('/', data$age_sample_collection),]
  
  duplicate_table$blood_dna_malkin_lab_ <- NULL
  duplicate_table$age_sample_collection <- NULL
  colnames(duplicate_table)[1:2] <- c('blood_dna_malkin_lab_', 'age_sample_collection')
  duplicate_table <- duplicate_table[, colnames(data)]
  data <- rbind(duplicate_table, data)
  return(data)
  
}

empty_table <- data.frame(matrix(ncol = ncol(clin), nrow = 0))
clin <- splitRows(clin, empty_table)


# some way of weighting for number of drugs


##########
# visualize our universe of genes per model
##########

# ##########
# # model outcome with just drugs
# ##########
# 
# # first group by treatment (need unqiue) and get counts for response and no response
# treat_counts <- treat_mod %>% group_by(treatment) %>% summarise(response = sum(outcome == 'Responded'),
#                                                                 no_response = sum(outcome == 'No response'))
# 
# # make column for total and percent responded
# treat_counts$total <- treat_counts$response + treat_counts$no_response
# treat_counts$per_response <- round((treat_counts$response/treat_counts$total)*100,2)
# 
# # make data frame
# treat_counts <- as.data.frame(treat_counts)
# 
# # remove respone, no responde, total
# treat_counts$response <- treat_counts$no_response <- treat_counts$total <- NULL
# temp <- treat_counts %>%
#   gather(treatment, drug) %>%
#   mutate(present = 1) %>%
#   select(-treatment) %>%
#   spread(treatment,present,fill = 0)
# 
# temp <- treat_counts %>% spread(treatment, treatment)
# 
# # fill with dummy 1 or zero
# temp[, 2:ncol(temp)][!is.na(temp[2:ncol(temp)])] <- 1
# temp[, 2:ncol(temp)][is.na(temp[2:ncol(temp)])] <- 0

##########
# for each drug or combination predict outcome
##########