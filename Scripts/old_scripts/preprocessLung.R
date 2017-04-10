##########################
# this script will read in and explore pdx data

##########
# initiate libraries
##########
library(dplyr)
library(lumi)
library(preprocessCore)

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
# read in id map
##########
id_map <- read.csv(paste0(data_folder, '/pdx_microarray_annotation.csv'), stringsAsFactors = F)
# recall that a.* is a different technology and we gotta do batch correction when combined with other data.

# remove white spaces from array_label
id_map$array_label <- trimws(id_map$array_label, which = "both")

##########
# Read in data from all seperate folders
##########

# Tsao_10.09.1_(DirectHyb_lungXeno)
a.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 10.09.17 (DirectHyb lungXeno 2 Raphael)',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

a.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 10.09.17 (DirectHyb lungXeno 2 Raphael)',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 11.11.16 Set1 (DASL lungXeno)
b.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 11.11.16 Set1 (DASL lungXeno)',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

b.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 11.11.16 Set1 (DASL lungXeno)',
                       '/sample_probe_nonnorm_FinalReport.csv'), 
                       stringsAsFactors = F)


# Tsao 11.11.16 Set2 (DASL lungXeno)
c.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 11.11.16 Set2 (DASL lungXeno)',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

c.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 11.11.16 Set2 (DASL lungXeno)',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 12.01.05 set1 (DASL 9 Notch 3Lung 12Pan)
d.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 12.01.05 set1 (DASL 9 Notch 3Lung 12Pan)',
                       '/sample-probe-nonnorm-FinalReport.txt' ))

d.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 12.01.05 set1 (DASL 9 Notch 3Lung 12Pan)',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 12.01.05 set2  (DASL 24 lung Xeno)
e.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 12.01.05 set2  (DASL 24 lung Xeno)',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

e.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 12.01.05 set2  (DASL 24 lung Xeno)',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 12.06.26 (DASL Part1 LungXeno)
f.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 12.06.26 (DASL Part1 LungXeno)',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

f.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 12.06.26 (DASL Part1 LungXeno)',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 12.06.26 Part2 (DSSL 16 LungXeno + 8 Niki)
g.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 12.06.26 Part2 (DSSL 16 LungXeno + 8 Niki)',
                       '/sample-Probe_nonnorm_FinalReport.txt' ))

g.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 12.06.26 Part2 (DSSL 16 LungXeno + 8 Niki)',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 12.08.15 (DASL Lisa Notch Xeno)
h.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 12.08.15 (DASL Lisa Notch Xeno)',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

h.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 12.08.15 (DASL Lisa Notch Xeno)',
                       '/sample_probe_nonnorm_FinalReport1.csv' ), 
                       stringsAsFactors = F)

#Tsao_12.11.30_Lisa
i.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 12.11.30_Lisa',
                       '/sample_probe_nonnorm_FinalReport.txt' ))

i.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 12.11.30_Lisa',
                       '/sample_probe_nonnorm_FinalReport.csv' ), 
                       stringsAsFactors = F)

# Tsao 13.02.07_PHLC164_ORF
j.lumi <- lumiR(paste0(micro_array_data, 
                       'Tsao 13.02.07_PHLC164_ORF',
                       '/sample_probe_nonnorm_FinalReport.txt' ))


j.lumi_raw <- read.csv(paste0(micro_array_data, 
                       'Tsao 13.02.07_PHLC164_ORF',
                       '/sample_probe_nonnorm_FinalReport.csv'), 
                       stringsAsFactors = F)

# ##########
# # function that cleans raw data puts genes in columns and applies log2 and scales using qunaitle normalization
# ##########
# 
# data <- b.lumi_raw
# 
# # remove first 7 rows
# data <- data[-(1:7),]
# 
# # colnames first row, delete first row
# colnames(data) <- data[1,]
# data <- data[-1,]
# 
# # recode Target 
# colnames(data)[1] <- 'gene'
# 
# # keep data that has colnames 'AVG'
# data <- data[ , grepl('AVG|gene|STDERR', colnames(data))]
# 
# # get ride of 'N' 'T'
# data <- data[ , !grepl('N.AVG|N.BEAD|N.Detection|T.AVG|T.BEAD|T.Detection', colnames(data))]
# 
# # put genes in rownames
# gene_names <- make.names(data$gene,unique = T)
# rownames(data) <- gene_names
# 
# # remove gene column
# data$gene <- NULL
# 
# for(i in 1:ncol(data)){
#   data[,i] <- as.numeric(data[,i])
#   print(i)
# }
# 
# # data <-new("ExpressionSet", exprs=as.matrix(data))
# # 
# # data <- lumiT(data, method = 'vst')
# 
# # apply vst
# for(i in seq(1, ncol(data), by=2) ) {
#   data_cols <- data[, i:(i+1) ] 
#   data_results <- vst(data_cols[,1], std = data_cols[,2])
#   data[, i:(i+1) ] <- data_results
# }
# 
# # remove standard deviation
# data <- data[, !grepl('STDER', colnames(data))]
# 
# # get colnames and rownames
# data_cols <- colnames(data)
# data_rows <- rownames(data)
# 
# data_quan <- lumiN(as.matrix(data), method = 'quantile')
# 
# # # apply quantile 
# # data_quan <- normalize.quantiles(as.matrix(data))
# 
# # get names 
# colnames(data_quan) <- data_cols
# rownames(data_quan) <- data_rows
# 
# # make data frame
# data_quan <- as.data.frame(data_quan)


##########
# function that takes a largebatch lumi object and raw data  and 
##########

getData <- function(lumi_object, raw_object, dasl)
{
  # get probe gene map
  gene_probe <- getGene(raw_object)
  
  # get expression matrices based on normalization method
  exp_quan <- getExp(lumi_object, method = 'quantile')
  exp_rsn <- getExp(lumi_object, method = 'rsn')
  exp_ssn <- getExp(lumi_object, method = 'ssn')
  exp_loess <- getExp(lumi_object, method = 'loess')
  exp_rank <- getExp(lumi_object, method = 'rankinvariant')
  
  # get gene column from gene_probe
  for(probe in gene_probe$ProbeID) {
    
    exp_quan$gene[which(row.names(exp_quan) == probe)] <- gene_probe$TargetID[which(gene_probe$ProbeID == probe)]
    exp_rsn$gene[which(row.names(exp_rsn) == probe)] <- gene_probe$TargetID[which(gene_probe$ProbeID == probe)]
    exp_ssn$gene[which(row.names(exp_ssn) == probe)] <- gene_probe$TargetID[which(gene_probe$ProbeID == probe)]
    exp_loess$gene[which(row.names(exp_loess) == probe)] <- gene_probe$TargetID[which(gene_probe$ProbeID == probe)]
    exp_rank$gene[which(row.names(exp_rank) == probe)] <- gene_probe$TargetID[which(gene_probe$ProbeID == probe)]
    
  }
  

  # group by gene and get mean 
  exp_quan <- exp_quan %>%
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  exp_rsn <- exp_rsn %>%
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  exp_ssn <- exp_ssn %>%
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  exp_loess <- exp_loess %>%
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  exp_rank <- exp_rank %>%
    group_by(gene) %>%
    summarise_each(funs(mean))
  
  # Transpose data 
  exp_quan <- tPose(exp_quan)
  exp_rsn <- tPose(exp_rsn)
  exp_ssn <- tPose(exp_ssn)
  exp_loess <- tPose(exp_loess)
  exp_rank <- tPose(exp_rank)
  
  # merge data 
  if (dasl) {
    exp_quan <- getMerger(exp_quan, dasl = T)
    exp_rsn <- getMerger(exp_rsn, dasl = T)
    exp_ssn <- getMerger(exp_ssn, dasl = T)
    exp_loess <- getMerger(exp_loess, dasl = T)
    exp_rank <- getMerger(exp_rank, dasl = T)
    
  } else {
    exp_quan <- getMerger(exp_quan, dasl = F)
    exp_rsn <- getMerger(exp_rsn, dasl = F)
    exp_ssn <- getMerger(exp_ssn, dasl = F)
    exp_loess <- getMerger(exp_loess, dasl = F)
    exp_rank <- getMerger(exp_rank, dasl = F)
  }
  
  return(list(exp_quan, exp_rsn, exp_ssn, exp_loess, exp_rank))
  
}

##########
# apply to data and get objects within list
##########

# a_lumi
a_lumi <- getData(a.lumi, a.lumi_raw, dasl = F)
a.lumi_quan <- a_lumi[[1]]
a.lumi_rsn <- a_lumi[[2]]
a.lumi_ssn <- a_lumi[[3]]
a.lumi_loess <- a_lumi[[4]]
a.lumi_rank <- a_lumi[[5]]

# b_lumi
b_lumi <- getData(b.lumi, b.lumi_raw, dasl = T)
b.lumi_quan <- b_lumi[[1]]
b.lumi_rsn <- b_lumi[[2]]
b.lumi_ssn <- b_lumi[[3]]
b.lumi_loess <- b_lumi[[4]]
b.lumi_rank <- b_lumi[[5]]

# c_lumi
c_lumi <- getData(c.lumi, c.lumi_raw, dasl = T)
c.lumi_quan <- c_lumi[[1]]
c.lumi_rsn <- c_lumi[[2]]
c.lumi_ssn <- c_lumi[[3]]
c.lumi_loess <- c_lumi[[4]]
c.lumi_rank <- c_lumi[[5]]

# d_lumi
d_lumi <- getData(d.lumi, d.lumi_raw, dasl = T)
d.lumi_quan <- d_lumi[[1]]
d.lumi_rsn <- d_lumi[[2]]
d.lumi_ssn <- d_lumi[[3]]
d.lumi_loess <- d_lumi[[4]]
d.lumi_rank <- d_lumi[[5]]

# e_lumi
e_lumi <- getData(e.lumi, e.lumi_raw, dasl = T)
e.lumi_quan <- e_lumi[[1]]
e.lumi_rsn <- e_lumi[[2]]
e.lumi_ssn <- e_lumi[[3]]
e.lumi_loess <- e_lumi[[4]]
e.lumi_rank <- e_lumi[[5]]

# f_lumi
f_lumi <- getData(f.lumi, f.lumi_raw, dasl = T)
f.lumi_quan <- f_lumi[[1]]
f.lumi_rsn <- f_lumi[[2]]
f.lumi_ssn <- f_lumi[[3]]
f.lumi_loess <- f_lumi[[4]]
f.lumi_rank <- f_lumi[[5]]

# g_lumi
g_lumi <- getData(g.lumi, g.lumi_raw, dasl = T)
g.lumi_quan <- g_lumi[[1]]
g.lumi_rsn <- g_lumi[[2]]
g.lumi_ssn <- g_lumi[[3]]
g.lumi_loess <- g_lumi[[4]]
g.lumi_rank <- g_lumi[[5]]

# h_lumi
h_lumi <- getData(h.lumi, h.lumi_raw, dasl = T)
h.lumi_quan <- h_lumi[[1]]
h.lumi_rsn <- h_lumi[[2]]
h.lumi_ssn <- h_lumi[[3]]
h.lumi_loess <- h_lumi[[4]]
h.lumi_rank <- h_lumi[[5]]

# i_lumi
i_lumi <- getData(i.lumi, i.lumi_raw, dasl = T)
i.lumi_quan <- i_lumi[[1]]
i.lumi_rsn <- i_lumi[[2]]
i.lumi_ssn <- i_lumi[[3]]
i.lumi_loess <- i_lumi[[4]]
i.lumi_rank <- i_lumi[[5]]

# j_lumi
j_lumi <- getData(j.lumi, j.lumi_raw, dasl = T)
j.lumi_quan <- j_lumi[[1]]
j.lumi_rsn <- j_lumi[[2]]
j.lumi_ssn <- j_lumi[[3]]
j.lumi_loess <- j_lumi[[4]]
j.lumi_rank <- j_lumi[[5]]

##########
# combine data a-i by preprocessing method
##########

# quan
quan.lumi <- combineLumi(a.lumi_quan, b.lumi_quan, c.lumi_quan, d.lumi_quan,
                         e.lumi_quan, f.lumi_quan, g.lumi_quan,
                         h.lumi_quan, i.lumi_quan, j.lumi_quan, combined = F)

quan.lumi_full <- combineLumi(a.lumi_quan, b.lumi_quan, c.lumi_quan, d.lumi_quan,
                              e.lumi_quan, f.lumi_quan, g.lumi_quan,
                              h.lumi_quan, i.lumi_quan, j.lumi_quan, combined = T)

# rsn
rsn.lumi <- combineLumi(a.lumi_rsn, b.lumi_rsn, c.lumi_rsn, d.lumi_rsn,
                         e.lumi_rsn, f.lumi_rsn, g.lumi_rsn,
                         h.lumi_rsn, i.lumi_rsn, j.lumi_rsn, combined = F)

rsn.lumi_full <- combineLumi(a.lumi_rsn, b.lumi_rsn, c.lumi_rsn, d.lumi_rsn,
                              e.lumi_rsn, f.lumi_rsn, g.lumi_rsn,
                              h.lumi_rsn, i.lumi_rsn, j.lumi_rsn, combined = T)

# ssn
ssn.lumi <- combineLumi(a.lumi_ssn, b.lumi_ssn, c.lumi_ssn, d.lumi_ssn,
                         e.lumi_ssn, f.lumi_ssn, g.lumi_ssn,
                         h.lumi_ssn, i.lumi_ssn, j.lumi_ssn, combined = F)

ssn.lumi_full <- combineLumi(a.lumi_ssn, b.lumi_ssn, c.lumi_ssn, d.lumi_ssn,
                              e.lumi_ssn, f.lumi_ssn, g.lumi_ssn,
                              h.lumi_ssn, i.lumi_ssn, j.lumi_ssn, combined = T)

# loess
loess.lumi <- combineLumi(a.lumi_loess, b.lumi_loess, c.lumi_loess, d.lumi_loess,
                         e.lumi_loess, f.lumi_loess, g.lumi_loess,
                         h.lumi_loess, i.lumi_loess, j.lumi_loess, combined = F)

loess.lumi_full <- combineLumi(a.lumi_loess, b.lumi_loess, c.lumi_loess, d.lumi_loess,
                              e.lumi_loess, f.lumi_loess, g.lumi_loess,
                              h.lumi_loess, i.lumi_loess, j.lumi_loess, combined = T)

# rank
rank.lumi <- combineLumi(a.lumi_rank, b.lumi_rank, c.lumi_rank, d.lumi_rank,
                         e.lumi_rank, f.lumi_rank, g.lumi_rank,
                         h.lumi_rank, i.lumi_rank, j.lumi_rank, combined = F)

rank.lumi_full <- combineLumi(a.lumi_rank, b.lumi_rank, c.lumi_rank, d.lumi_rank,
                              e.lumi_rank, f.lumi_rank, g.lumi_rank,
                              h.lumi_rank, i.lumi_rank, j.lumi_rank, combined = T)



# log2 raw and do quantile manual and see if it is the same as this pipeline.
# check to see if there are any weird ones (h)

# what is quantile doing 

# mircroarray vs gene expression

# save data
saveRDS(quan.lumi, paste0(data_folder, 'quan.lumi.rda'))
saveRDS(quan.lumi_full, paste0(data_folder, 'quan.lumi_full.rda'))
saveRDS(rsn.lumi, paste0(data_folder, 'rsn.lumi.rda'))
saveRDS(rsn.lumi_full, paste0(data_folder, 'rsn.lumi_full.rda'))
saveRDS(ssn.lumi, paste0(data_folder, 'ssn.lumi.rda'))
saveRDS(ssn.lumi_full, paste0(data_folder, 'ssn.lumi_full.rda'))
saveRDS(loess.lumi, paste0(data_folder, 'loess.lumi.rda'))
saveRDS(loess.lumi_full, paste0(data_folder, 'loess.lumi_full.rda'))
saveRDS(rank.lumi, paste0(data_folder, 'rank.lumi.rda'))
saveRDS(rank.lumi_full, paste0(data_folder, 'rank.lumi_full.rda'))
