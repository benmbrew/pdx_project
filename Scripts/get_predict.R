##########################
# this is the 4th step in the pipeline - predict


##########
# initialize libraries
##########
library(dplyr)
library(lumi)
library(preprocessCore)
library(caret)
library(sva)
library(nnet)
##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdxSNF/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# read in mRECIST for classification
##########
nov_outcome <- read.csv(paste0(data_folder, 'pdxe_outcome.csv'))

# transpose 
nov_outcome <- as.data.frame(t(nov_outcome), stringsAsFactors = F)

# get colnames
colnames(nov_outcome) <- nov_outcome[1,]
nov_outcome <- nov_outcome[-1,]
nov_feat <- colnames(nov_outcome)
nov_outcome$id <- rownames(nov_outcome)
nov_outcome <- nov_outcome[, c('id', nov_feat)]

#########
# read in dasl and full
#########
data_full <- readRDS(paste0(data_folder, '/data_full.rds'))
data_dasl <- readRDS(paste0(data_folder, '/data_dasl.rds'))

#########
# seperate data
#########
dasl <- data_dasl[data_dasl$batch != 'nov',]
dasl_nov <- data_dasl[data_dasl$batch == 'nov',]

full <- data_full[data_full$batch != 'nov',]
full_nov <- data_full[data_full$batch == 'nov',]

# join novartis with outcome
dasl_nov <- inner_join(nov_outcome, dasl_nov ,by = c('id' = 'model_id'))
full_nov <- inner_join(nov_outcome, full_nov ,by = c('id' = 'model_id'))

# library(reshape2)
# nov_drugs <- dasl_nov[, 1:64]
# nov_drugs_melt <- melt(nov_drugs, id.vars = 'id')
# # recode valu
# nov_drugs_melt <- nov_drugs_melt[complete.cases(nov_drugs_melt),]
# 
# nov_drugs_melt$value <-  ifelse(grepl('CR|PR', nov_drugs_melt$value), 'yes', 'no')
# 
# # # group by id and drug and get counts
# # 
# # drug_counts <- nov_drugs_melt %>%
# #   group_by(variable) %>%
# #   summarise(counts = n())
# # 
# # drug_counts_res <- nov_drugs_melt %>%
# #   group_by(variable) %>%
# #   summarise(response = sum(value == 'yes'),
# #             no_response = sum(value == 'no'))
# 
# 

# # list of drugs 
# drug_list <- colnames(mod_norm)[ 2:63]


##########
# predict on novartis classifier
##########
# CR - complete response
# PR - partial response,
# PD - progressive disease
# SD - stable disease
# Random forest model 

rfPredictFac <- function(model_data,
                         drug_name,
                         cutoff,
                         iterations) 
{
  
  model <- list()
  best_features <- list()
  importance <- list()
  train.predictions <- list()
  test.predictions <- list()
  train.ground_truth <- list()
  test.ground_truth <- list()
  test_acc <- list()
  test_stats  <- list()
  
  selected_features <- names(model_data)[65:ncol(model_data)]
  
  model_data <- model_data[, c('id', drug_name, selected_features )]
  model_data <- model_data[complete.cases(model_data),]
  model_data[, drug_name] <- ifelse(grepl('CR|PR', model_data[, drug_name]), 'yes', 'no')
  
  # balance class
  neg <- length(which(model_data[, drug_name] == 'no'))
  pos <- length(which(model_data[, drug_name] == 'yes'))
  stopifnot(pos >10)
  great <- pos > neg
  if(great == F){
    diff <- neg - pos
    remove_index <- which(model_data[, drug_name] == 'no') 
    amount_yes <- length(which(model_data[, drug_name] == 'yes'))
    amount_no <- length(which(model_data[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_no - amount_yes) - 10, replace = F )
    model_data <- model_data[-remove_index,]
  } else {
    diff <- pos - neg
    remove_index <- which(model_data[, drug_name] == 'yes') 
    amount_yes <- length(which(model_data[, drug_name] == 'yes'))
    amount_no <- length(which(model_data[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_yes - amount_no)+5, replace = F )
    model_data <- model_data[-remove_index,]
    
  }
  
  new_neg <- length(which(model_data[, drug_name] == 'no'))
  new_pos <- length(which(model_data[, drug_name] == 'yes'))
  
  dims <- dim(model_data)
  
  
  for (i in 1:iterations) {
    
    set.seed(i)
    
    train_index <- sample(nrow(model_data), nrow(model_data) *cutoff, replace = F)
    
    # determines how you train the model.
    NFOLDS <- 2
    fitControl <- trainControl( 
      method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
      number = min(10, NFOLDS),
      classProbs = TRUE,
      repeats = 1,
      allowParallel = TRUE,
      summaryFunction = twoClassSummary
      
    )
    
    y = as.factor(model_data[,drug_name][train_index])
    
    mtry <- sqrt(ncol(model_data[train_index, selected_features]))
    tunegrid <- expand.grid(.mtry=mtry)
    
    model[[i]] <- train(x = model_data[train_index, selected_features]
                        , y = y
                        , method = "rf"
                        , trControl = fitControl
                        , tuneGrid = tunegrid
                        , importance = T
                        , verbose = FALSE)
    
    temp <- varImp(model[[i]])[[1]]
    importance[[i]] <- cbind(rownames(temp), temp$X1)
    
    test.predictions[[i]] <- predict(model[[i]] 
                                     , newdata = model_data[-train_index, selected_features]
                                     , type = "prob")
    
    train.predictions[[i]] <- predict(model[[i]] 
                                      , newdata = model_data[train_index, selected_features]
                                      ,type = "prob")
    
    
    
    
    
    train.ground_truth[[i]] <- as.factor(make.names(model_data[, drug_name][train_index]))
    test.ground_truth[[i]] <- as.factor(make.names(model_data[, drug_name][-train_index]))
    
    # For age of diagnosis
    # Accuracy
    test_acc[[i]] <- sum(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)] == test.ground_truth[[i]]) / dim(test.predictions[[i]])[1]
    # Confustion Matrix
    test_stats[[i]] <- confusionMatrix(levels(test.ground_truth[[i]])[apply(test.predictions[[i]], 1, which.is.max)], test.ground_truth[[i]])
    
    
    print(i)
  }
  
  return(list(train.predictions, test.predictions, train.ground_truth, test.ground_truth,
              test_acc, test_stats, model, importance, dims, new_pos, new_neg))
  
  
}
#########
# function for conMatrix
#########
conMatrix <- function(results) 
{
  
  # test acc for age of diagnosis
  acc <- mean(unlist(results[[5]]))
  
  
  # confustion matrix age of diagnosis 10
  iterations <- 10
  temp <- list()
  for (i in 1:10){
    temp[[i]] <- results[[6]][[i]]$table
  }
  mat <- unlist(temp)
  new_mat <- matrix(, 2, 2)
  
  mat_index <- seq(1, length(mat), 4)
  
  new_mat[1,1] <- sum(mat[mat_index])/iterations
  new_mat[2,1] <- sum(mat[mat_index + 1])/iterations
  new_mat[1,2] <- sum(mat[mat_index + 2])/iterations
  new_mat[2,2] <- sum(mat[mat_index + 3])/iterations
  
  return(list(new_mat, acc))
  
}


# list of drugs 
drug_list <- colnames(dasl_nov)[2:63]

##########################
drug <- print(drug_list[[2]])
# "binimetinib"
bini_norm <- rfPredictFac(dasl_nov, 
                          drug, 
                          cutoff = 0.7, 
                          iterations = 10)

bini_tab_norm <- print(conMatrix(bini_norm)[[1]])
bini_acc_norm <- print(conMatrix(bini_norm)[[2]])

bini_norm[[9]] #dims
bini_norm[[10]] #pos
bini_norm[[11]] #neg

bini_full <- rfPredictFac(full_nov, 
                          drug, 
                          cutoff = 0.7, 
                          iterations = 10)
bini_tab_full <- print(conMatrix(bini_full)[[1]])
bini_acc_full <- print(conMatrix(bini_full)[[2]])

bini_full[[9]] # dims
bini_full[[10]] # pos
bini_full[[11]]# neg
#this one
##########################
drug <- print(drug_list[[3]])
# "BYL719_LJM716"
BYL719_LJM716_norm <- rfPredictFac(dasl_nov, 
                                   drug, 
                                   cutoff = 0.7, 
                                   iterations = 10)
BYL719_LJM716_tab_norm <- print(conMatrix(BYL719_LJM716_norm)[[1]])
BYL719_LJM716_acc_norm <- print(conMatrix(BYL719_LJM716_norm)[[2]])

BYL719_LJM716_norm[[9]] # dims 
BYL719_LJM716_norm[[10]] # pos
BYL719_LJM716_norm[[11]] # neg

BYL719_LJM716_full <- rfPredictFac(full_nov, 
                                   drug, 
                                   cutoff = 0.7, 
                                   iterations = 10)
BYL719_LJM716_tab_full <- print(conMatrix(BYL719_LJM716_full)[[1]])
BYL719_LJM716_acc_full <- print(conMatrix(BYL719_LJM716_full)[[2]])

BYL719_LJM716_full[[9]] # dims
BYL719_LJM716_full[[10]] # pos
BYL719_LJM716_full[[11]] # neg

##########################
drug <- print(drug_list[[6]])
# "BYL719_+_LJM716"
BYL719_LJM716_norm <- rfPredictFac(dasl_nov, 
                                   drug, 
                                   cutoff = 0.7, 
                                   iterations = 10)
BYL719_LJM716_tab_norm <- print(conMatrix(BYL719_LJM716_norm)[[1]])
BYL719_LJM716_acc_norm <- print(conMatrix(BYL719_LJM716_norm)[[2]])

BYL719_LJM716_norm[[9]] # dims 
BYL719_LJM716_norm[[10]] # pos
BYL719_LJM716_norm[[11]] # neg

BYL719_LJM716_full <- rfPredictFac(full_nov, 
                                   drug, 
                                   cutoff = 0.7, 
                                   iterations = 10)
BYL719_LJM716_tab_full <- print(conMatrix(BYL719_LJM716_full)[[1]])
BYL719_LJM716_acc_full <- print(conMatrix(BYL719_LJM716_full)[[2]])

BYL719_LJM716_full[[9]] # dims
BYL719_LJM716_full[[10]] # pos
BYL719_LJM716_full[[11]] # neg

##########################
drug <- print(drug_list[[7]])
# "CLR457"
CLR457_norm <- rfPredictFac(dasl_nov, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
CLR457_tab_norm <- print(conMatrix(CLR457_norm)[[1]])
CLR457_acc_norm <- print(conMatrix(CLR457_norm)[[2]])

CLR457_norm[[9]] # dims 
CLR457_norm[[10]] # pos
CLR457_norm[[11]] # neg

CLR457_full <- rfPredictFac(full_nov, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
CLR457_tab_full <- print(conMatrix(CLR457_full)[[1]])
CLR457_acc_full <- print(conMatrix(CLR457_full)[[2]])

CLR457_full[[9]] # dims
CLR457_full[[10]] # pos
CLR457_full[[11]] # neg


##########################
drug <- print(drug_list[[14]])
# "LEE011"
LEE011_norm <- rfPredictFac(dasl_nov, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
LEE011_tab_norm <- print(conMatrix(LEE011_norm)[[1]])
LEE011_acc_norm <- print(conMatrix(LEE011_norm)[[2]])

LEE011_norm[[9]] # dims 
LEE011_norm[[10]] # pos
LEE011_norm[[11]] # neg

LEE011_full <- rfPredictFac(full_nov, 
                            drug, 
                            cutoff = 0.7, 
                            iterations = 10)
LEE011_tab_full <- print(conMatrix(LEE011_full)[[1]])
LEE011_acc_full <- orint(conMatrix(LEE011_full)[[2]])

LEE011_full[[9]] # dims
LEE011_full[[10]] # pos
LEE011_full[[11]] # neg


##########################
drug <- print(drug_list[[15]])
# "LEE011_everolimus"
LEE011_everolimus_norm <- rfPredictFac(dasl_nov, 
                                       drug, 
                                       cutoff = 0.7, 
                                       iterations = 10)
LEE011_everolimus_tab_norm <- print(conMatrix(LEE011_everolimus_norm)[[1]])
LEE011_everolimus_acc_norm <- print(conMatrix(LEE011_everolimus_norm)[[2]])

LEE011_everolimus_norm[[9]] # dims 
LEE011_everolimus_norm[[10]] # pos
LEE011_everolimus_norm[[11]] # neg

LEE011_everolimus_full <- rfPredictFac(full_nov, 
                                       drug, 
                                       cutoff = 0.7, 
                                       iterations = 10)
LEE011_everolimus_tab_full <- print(conMatrix(LEE011_everolimus_full)[[1]])
LEE011_everolimus_acc_full <- print(conMatrix(LEE011_everolimus_full)[[2]])

LEE011_everolimus_full[[9]] # dims
LEE011_everolimus_full[[10]] # pos
LEE011_everolimus_full[[11]] # neg
#this one

##########################
drug <- print(drug_list[[23]])
# "BYL719_binimetinib"
BYL719_binimetinib_norm <- rfPredictFac(dasl_nov, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
BYL719_binimetinib_tab_norm <- print(conMatrix(BYL719_binimetinib_norm)[[1]])
BYL719_binimetinib_acc_norm <- print(conMatrix(BYL719_binimetinib_norm)[[2]])

BYL719_binimetinib_norm[[9]] # dims 
BYL719_binimetinib_norm[[10]] # pos
BYL719_binimetinib_norm[[11]] # neg

BYL719_binimetinib_full <- rfPredictFac(full_nov, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
BYL719_binimetinib_tab_full <- print(conMatrix(BYL719_binimetinib_full)[[1]])
BYL719_binimetinib_acc_full <- print(conMatrix(BYL719_binimetinib_full)[[2]])

BYL719_binimetinib_full[[9]] # dims
BYL719_binimetinib_full[[10]] # pos
BYL719_binimetinib_full[[11]] # neg
# this one

#########################
drug <- print(drug_list[[37]])
# "BKM120_binimetinib"
BKM120_binimetinib_norm <- rfPredictFac(dasl_nov, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
BKM120_binimetinib_tab_norm <- print(conMatrix(BKM120_binimetinib_norm)[[1]])
BKM120_binimetinib_acc_norm <- print(conMatrix(BKM120_binimetinib_norm)[[2]])

BKM120_binimetinib_norm[[9]] # dims 
BKM120_binimetinib_norm[[10]] # pos
BKM120_binimetinib_norm[[11]] # neg

BKM120_binimetinib_full <- rfPredictFac(full_nov, 
                                        drug, 
                                        cutoff = 0.7, 
                                        iterations = 10)
BKM120_binimetinib_tab_full <- print(conMatrix(BKM120_binimetinib_full)[[1]])
BKM120_binimetinib_acc_full <- print(conMatrix(BKM120_binimetinib_full)[[2]])

BKM120_binimetinib_full[[9]] # dims
BKM120_binimetinib_full[[10]] # pos
BKM120_binimetinib_full[[11]] # neg

#########################
drug <- print(drug_list[[40]])
# gemcitabine-50mpk
gemcitabine_50mpk_norm <- rfPredictFac(dasl_nov, 
                                       drug, 
                                       cutoff = 0.7, 
                                       iterations = 10)
gemcitabine_50mpk_tab_norm <- print(conMatrix(gemcitabine_50mpk_norm)[[1]])
gemcitabine_50mpk_acc_norm <- print(conMatrix(gemcitabine_50mpk_norm)[[2]])

gemcitabine_50mpk_norm[[9]] # dims 
gemcitabine_50mpk_norm[[10]] # pos
gemcitabine_50mpk_norm[[11]] # neg

gemcitabine_50mpk_full <- rfPredictFac(full_nov, 
                                       drug, 
                                       cutoff = 0.7, 
                                       iterations = 10)
gemcitabine_50mpk_tab_full <- print(conMatrix(gemcitabine_50mpk_full)[[1]])
gemcitabine_50mpk_acc_full <- print(conMatrix(gemcitabine_50mpk_full)[[2]])

gemcitabine_50mpk_full[[9]] # dims
gemcitabine_50mpk_full[[10]] # pos
gemcitabine_50mpk_full[[11]] # neg

#########
# predict on ming pdx
##########
# Random forest model 
rfPredictTest <- function(model_data,
                          drug_name) 
{
  
  dims <- dim(model_data)
  novartis <- model_data[grepl('X.', model_data$id),]
  pdxe <-model_data[!grepl('X.', model_data$id),]
  dasl_nov <- inner_join(nov_outcome, novartis, by = 'id')
  
  selected_features <- names(dasl_nov)[65:ncol(dasl_nov)]
  
  dasl_nov <- dasl_nov[, c('id', drug_name, selected_features )]
  dasl_nov <- dasl_nov[complete.cases(dasl_nov),]
  dasl_nov[, drug_name] <- ifelse(dasl_nov[, drug_name] == 'CR', 'yes', 
                                  ifelse(dasl_nov[, drug_name] == 'PR', 'yes', 'no'))
  
  # balance class
  neg <- length(which(dasl_nov[, drug_name] == 'no'))
  pos <- length(which(dasl_nov[, drug_name] == 'yes'))
  stopifnot(pos >10)
  great <- pos > neg
  if(great == F){
    diff <- neg - pos
    remove_index <- which(dasl_nov[, drug_name] == 'no') 
    amount_yes <- length(which(dasl_nov[, drug_name] == 'yes'))
    amount_no <- length(which(dasl_nov[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_no - amount_yes) - 10, replace = F )
    dasl_nov <- dasl_nov[-remove_index,]
  } else {
    diff <- pos - neg
    remove_index <- which(dasl_nov[, drug_name] == 'yes') 
    amount_yes <- length(which(dasl_nov[, drug_name] == 'yes'))
    amount_no <- length(which(dasl_nov[, drug_name] == 'no'))
    
    
    remove_index <- sample(remove_index, (amount_yes - amount_no)+5, replace = F )
    dasl_nov <- dasl_nov[-remove_index,]
    
  }
  
  new_neg <- length(which(dasl_nov[, drug_name] == 'no'))
  new_pos <- length(which(dasl_nov[, drug_name] == 'yes'))
  
  dims <- dim(dasl_nov)
  
  # determines how you train the model.
  NFOLDS <- 2
  fitControl <- trainControl( 
    method = "repeatedcv",  # could train on boostrap resample, here use repeated cross validation.
    number = min(10, NFOLDS),
    classProbs = TRUE,
    repeats = 1,
    allowParallel = TRUE,
    summaryFunction = twoClassSummary
    
  )
  
  y = as.factor(dasl_nov[,drug_name])
  
  mtry <- sqrt(ncol(dasl_nov))
  tunegrid <- expand.grid(.mtry=mtry)
  
  model <- train(x = dasl_nov[, selected_features]
                 , y = y
                 , method = "rf"
                 , trControl = fitControl
                 , tuneGrid = tunegrid
                 , importance = T
                 , verbose = FALSE)
  
  temp <- varImp(model)[[1]]
  importance <- cbind(rownames(temp), temp$X1)
  
  test.predictions <- predict(model 
                              , newdata = pdxe[, selected_features]
                              , type = "prob")
  
  
  
  return(list(test.predictions, model, importance, dims, new_pos, new_neg))
  
}

#2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[2]])
# "binimetinib"

binimetinib_norm <- rfPredictTest(dat_norm_cor, 
                                  drug)[[1]]
binimetinib_full <- rfPredictTest(dat_full_cor, 
                                  drug)[[1]]

result <- binimetinib_norm[order(binimetinib_norm$yes), ]


#2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[3]])
# BKM120

BKM120_norm <- rfPredictTest(dat_norm_cor, 
                             drug)
BKM120_full <- rfPredictTest(dat_full_cor, 
                             drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[6]])
# "BYL719_LJM716"

BYL719_LJM716_norm <- rfPredictTest(dat_norm_cor, 
                                    drug)
BYL719_LJM716_full <- rfPredictTest(dat_full_cor, 
                                    drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[[7]])
# "CLR457"

CLR457_norm <- rfPredictTest(dat_norm_cor, 
                             drug)
CLR457_full <- rfPredictTest(dat_full_cor, 
                             drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[14])
# "LEE011"

LEE011_norm <- rfPredictTest(dat_norm_cor, 
                             drug)
LEE011_full <- rfPredictTest(dat_full_cor, 
                             drug)

# 2,3,6,7,14,23,27,40
############################
drug <- print(drug_list[23])
# "BYL719_binimetinib"

BYL719_binimetinib_norm <- rfPredictTest(dat_norm_cor, 
                                         drug)
LBYL719_binimetinib_full <- rfPredictTest(dat_full_cor, 
                                          drug)

# 2,3,6,7,14,23,37,40
############################
drug <- print(drug_list[37])
# "BKM120_binimetinib"

BKM120_binimetinib_norm <- rfPredictTest(dat_norm_cor, 
                                         drug)
BKM120_binimetinib_full <- rfPredictTest(dat_full_cor, 
                                         drug)

# 2,3,6,7,14,23,37,40
############################
drug <- print(drug_list[40])
# "gemcitabine-50mpk"

gemcitabine_50mpk_norm <- rfPredictTest(dat_norm_cor, 
                                        drug)
gemcitabine_50mpk_full <- rfPredictTest(dat_full_cor, 
                                        drug)



