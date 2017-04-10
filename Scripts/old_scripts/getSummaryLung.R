###################
# this script will examine and get summary stats for mings lung data

##########
# initiate libraries
##########
library(dplyr)
library(lumi)
library(sva)

##########
# initiate folder
##########
home_folder <- '/home/benbrew/Documents/'
project_folder <- paste0(home_folder, 'pdxSNF/')
data_folder <- paste0(project_folder, 'Data/')
micro_array_data <- paste0(data_folder, 'microarray/')
code_folder <- paste0(project_folder, 'Code/')

##########
# source functions
##########
source(paste0(code_folder, 'functions.R'))

##########
# read in quantile normalized data since basic clinical attributes are similar across methods
##########

# read data
quan <- readRDS(paste0(data_folder, 'quan.lumi.rda'))
quan_full <- readRDS(paste0(data_folder, 'quan.lumi_full.rda'))

##########
# create new factor for batches
##########
getBatch <- function(data)
{
  batch_factor <- c("a", "b", "c", "d", "e", 
                    "f", "g", "h", "i", "j")
  
  data$project_title <- batch_factor[as.factor(data$project_title)]
  
  
  return(data)
}

quan <- getBatch(quan)
quan_full <- getBatch(quan_full)

# sort data 
quan_full <- quan_full[order(quan_full$model_id, decreasing = F),]

##########
# get relevant columns
##########
features <- colnames(quan_full[, 6:ncol(quan_full)])
quan_clin <- quan_full[, c('model_id', 'passage', 'project_title', features)]

# get unique
length(unique(quan_clin$model_id))

##########
# Grab controls data and map them (also treat model data for duplicates)
##########

# first make numeric 
quan_clin[,4:ncol(quan_clin)] <- apply(quan_clin[,4:ncol(quan_clin)], 2, as.numeric)

# make group by variables factors
quan_clin$model_id <- as.character(quan_clin$model_id)
quan_clin$passage <- as.character(quan_clin$passage)
quan_clin$project_title <- as.character(quan_clin$project_title)


# avg genes that have same model, passage, and project_title
quan_avg <- quan_clin %>%
  group_by(model_id, passage, project_title) %>%
  summarise_each(funs(mean))

# remove perfect duplicates (probably first few coulmns)
quan_dup <- quan_avg[c("7A5", "A1BG", "A1CF", "A26C3", "A2BP1")]
quan_dup <- quan_avg[!duplicated(quan_dup),]

# remove !complete.cases
quan_dup <- quan_dup[complete.cases(quan_dup),]
quan_dup <- as.data.frame(quan_dup)
#(dup pipeline end)

# get same model id and batch, but different project title
duplicates <- quan_dup[c('model_id', 'passage')]
temp <- duplicates[duplicated(duplicates, fromLast = T),]
temp.1 <- paste(unique(temp$model_id), collapse = '|')

# remove project title
con_list <- list()
for (i in 1:nrow(quan_dup)){
  sub_dat <- quan_dup[i,]
  if(any(grepl(sub_dat$model_id, temp.1))) {
    con_list[[i]] <- sub_dat
  }
  print(i)
  
}

controls <- do.call(rbind, con_list)

# remove model_id = 8
controls <- controls[controls$model_id != '8',]

##########
# create function that takes model id, and plots controls against each other
##########
controls_data <- controls
column_name <- '116'
plotControls <- function(controls_data, column_name)
{
  plot_data <- controls_data[controls_data$model_id == column_name,]
  x_axis <- log(as.numeric(plot_data[1, 4:ncol(plot_data)]))
  y_axis <- as.numeric(plot_data[2, 4:ncol(plot_data)])
  smoothScatter(x_axis, y_axis, xlab = "Time 1", ylab = "Time 2")
}
# plot for first patient
orig_1 <- as.numeric(beta_orig[1, 3:ncol(beta_orig)])
control_1 <- as.numeric(beta_control[1, 3:ncol(beta_orig)])

smoothScatter(orig_1, control_1, main = '1st Sample',
              xlab = 'Original', ylab = 'Control')


##########
# recode project.title
##########
quan_clin$project_title <- ifelse(quan_clin$project_title == 'a', 'a', 'b')

##########
# group by passage and get counts 
##########
passages <- quan_clin %>%
  filter(!is.na(passage)) %>%
  group_by(passage) %>%
  summarise(counts = n())

##########
# group by project title and get counts 
##########
batch <- quan_clin %>%
  filter(!is.na(project_title)) %>%
  group_by(project_title) %>%
  summarise(counts = n())

##########
# group by project_title.1 and passage 
##########
passage_batch <- quan_clin %>%
  filter(!is.na(passage)) %>%
  filter(!is.na(project_title)) %>%
  group_by(project_title, passage) %>%
  summarise(counts = n())

##########
# group by project_title.1, passage, model_id
##########
passage_batch_model <- quan_clin %>%
  filter(!is.na(passage)) %>%
  filter(!is.na(project_title)) %>%
  filter(!is.na(model_id)) %>%
  group_by(project_title, passage, model_id) %>%
  summarise(counts = n())

#check unique now
length(unique(passage_batch_model$model_id))

# duplicates here will likely have pefect matches. sort data
passage_batch_model <- passage_batch_model[order(passage_batch_model$model_id, decreasing = F),]
