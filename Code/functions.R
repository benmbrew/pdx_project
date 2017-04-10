# ##########
# # this function will find the intersection of samples across modalities 
# ##########
# 
# # column intersection 
# columnIntersection <- function(L)
# {
#   n <- length(L)
#   cnames <- colnames(L[[1]])
#   
#   for (i in 1:n) {
#     cnames <- intersect(cnames, colnames(L[[i]]))
#   }
#   
#   if (length(cnames) == 0) {
#     return(NULL)
#   }
#   
#   cnames <- sort(cnames)
#   instersectL <- vector("list", n)
#   
#   for (i in 1:n) {
#     instersectL[[i]] <- L[[i]][, match(cnames, colnames(L[[i]])), drop=FALSE]
#   }
#   
#   return(list(instersectL, cnames))
# }
# 
# 
# ##########
# # normalize 
# ##########
# rowStatistics <- function(cases) 
# {
#   numViews <- length(cases)
#   rowStats <- vector("list", numViews)
#   
#   for (v in 1:numViews) {
#     # Calculate the row means and standard deviations
#     rowMean <- apply(cases[[v]], 1, mean, na.rm=TRUE)
#     rowSd <- apply(cases[[v]], 1, sd, na.rm=TRUE)
#     constantInd <- rowSd==0
#     rowSd[constantInd] <- 1
#     rowStats[[v]] <- list(mean=rowMean, sd=rowSd, ind=constantInd)
#   }
#   
#   return(rowStats)
# }
# 
# normalizeData <- function(data, stat) {
#   for (v in 1:length(data)) {
#     data[[v]] <- (data[[v]] - stat[[v]]$mean) / stat[[v]]$sd
#     data[[v]] <- data[[v]][!stat[[v]]$ind, ]
#   }
#   
#   return(data)
# }
# 
# 
# #########
# # this will perform SNF clustering 
# #########
# 
# SNFClustering <- function(data, numClus, sampleRows) 
# {
#   # Transpose data for the dist function
#   # The dist function takes distances between rows
#   sampleRowsRequired <- TRUE
#   transposeData <- sampleRows != sampleRowsRequired
#   if (transposeData) {
#     data <- lapply(data, t)
#   }
#   
#   # Calculate the distance between samples
#   distances <- lapply(data, function(x) as.matrix(dist(x)))
#   
#   # Convert the distances to affinities
#   affinities <- lapply(distances, affinityMatrix)
#   
#   # Fuse the affinity matrices
#   fusedMatrix <- SNF(affinities)
#   
#   # Cluster the fused matrix
#   labels <- spectralClustering(fusedMatrix, numClus)
#   
#   return(labels)
# }

##########
# this function will get gene info
#########

getGene <- function(raw_data) 
{
  raw_data <- raw_data[-c(1:7),]
  colnames(raw_data) <- raw_data[1, ]
  temp <- raw_data[1,]
  temp <- as.vector(temp)
  colnames(raw_data) <- temp
  raw_data <- raw_data[-1, ] 
  gene_probe <- raw_data[,1:2]
  return(gene_probe)
}

#########
# this function will an expression matrix from lumi object
#########
getExp <- function(lumi_object, method) 
{
  # ## Quality control based on the raw data
  # lumi.Q <- lumiQ(lumi_object) # (optional)
  ## Do default VST variance stabilizing transform
  lumi.T <- lumiT(lumi_object, method = 'log2')
  ## Do lumi between microarray normaliazation
  lumi.N <- lumiN(lumi.T, method = method)
  # ## Quality control after normalization
  # lumi.N.Q <- lumiQ(lumi.N) # (optional)
  # # 5.2 Identify differentiate genes
  # Identify the differentiated genes based on moderated t-test using limma.
  # Retrieve the normalized data
  dataMatrix <- exprs(lumi.N)
  dataMatrix <- as.data.frame(dataMatrix)
  
  return(dataMatrix)
  
}

##########
# function that transposes data and sets columns, creates id variable
##########

tPose <- function(data)
{
  data <- as.data.frame(t(data), stringsAsFactors = F)
  colnames(data) <- data[1,]
  data <- data[-1,]
  data$id <- rownames(data)
  return(data)
}

##########
# function that mergers data with id_map
##########
getMerger <- function(data, dasl) 
{
  features <- colnames(data[, (1:ncol(data) - 1)])
  data <- data[, c("id", features)]
  
  if (dasl) {
    data <- left_join(data, id_map, c("id" = "array_label"))
    data <- data[, c("model_id" , "array_barcode", "passage" ,"project_title" , features)]
    # data <- data[!duplicated(data$id),]
    data <- data[!is.na(data$array_barcode),]
    
  } else {
    data <- left_join(data, id_map, c("id" = "array_barcode")) 
    data <- data[, c("model_id", "id","passage" , "project_title", features)]
    names(data)[2] <- 'array_barcode'
    data <- data[!is.na(data$array_barcode),]
    # data <- data[!duplicated(data$id),]
    
  }
  
  return(data)
}

##########
# takes a-i and combines them for each preprocessing method
##########

combineLumi <- function(a.lumi, b.lumi, c.lumi, d.lumi, e.lumi,
                        f.lumi, g.lumi, h.lumi, i.lumi, j.lumi, combined) 
{
  if (combined) {
    
    feature_intersection <- Reduce(intersect, list(
      colnames(a.lumi)[4:ncol(a.lumi)], 
      colnames(b.lumi)[4:ncol(b.lumi)],
      colnames(c.lumi)[4:ncol(c.lumi)],
      colnames(d.lumi)[4:ncol(d.lumi)],
      colnames(e.lumi)[4:ncol(e.lumi)],
      colnames(f.lumi)[4:ncol(f.lumi)],
      colnames(g.lumi)[4:ncol(g.lumi)],
      colnames(h.lumi)[4:ncol(h.lumi)],
      colnames(i.lumi)[4:ncol(i.lumi)],
      colnames(j.lumi)[4:ncol(j.lumi)])
    )
    
    a.lumi <- a.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    b.lumi <- b.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    c.lumi <- c.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    d.lumi <- d.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    e.lumi <- e.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    f.lumi <- f.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    g.lumi <- g.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    h.lumi <- h.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    i.lumi <- i.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    j.lumi <- j.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    
    lumi <- rbind(a.lumi, b.lumi, c.lumi, d.lumi,
                  e.lumi, f.lumi, g.lumi, h.lumi,
                  i.lumi, j.lumi)
    
  } else {
    
    feature_intersection <- Reduce(intersect, list(
      colnames(b.lumi)[4:ncol(b.lumi)],
      colnames(c.lumi)[4:ncol(c.lumi)],
      colnames(d.lumi)[4:ncol(d.lumi)],
      colnames(e.lumi)[4:ncol(e.lumi)],
      colnames(f.lumi)[4:ncol(f.lumi)],
      colnames(g.lumi)[4:ncol(g.lumi)],
      colnames(h.lumi)[4:ncol(h.lumi)],
      colnames(i.lumi)[4:ncol(i.lumi)],
      colnames(j.lumi)[4:ncol(j.lumi)])
    )
    
    b.lumi <- b.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    c.lumi <- c.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    d.lumi <- d.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    e.lumi <- e.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    f.lumi <- f.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    g.lumi <- g.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    h.lumi <- h.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    i.lumi <- i.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    j.lumi <- j.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    
    lumi <- rbind(b.lumi, c.lumi, d.lumi,
                  e.lumi, f.lumi, g.lumi, 
                  h.lumi, i.lumi, j.lumi)
  }
  
  return(lumi)
  
}


##########
# takes a-i and combines them for each preprocessing method
##########

combineLumiRaw <- function(a.lumi, b.lumi, c.lumi, d.lumi, e.lumi,
                           f.lumi, g.lumi, h.lumi, i.lumi, j.lumi, combined) 
{
  
  if (combined) {
    
    feature_intersection <- Reduce(intersect, list(
      colnames(a.lumi)[4:ncol(a.lumi)], 
      colnames(b.lumi)[4:ncol(b.lumi)],
      colnames(c.lumi)[4:ncol(c.lumi)],
      colnames(d.lumi)[4:ncol(d.lumi)],
      colnames(e.lumi)[4:ncol(e.lumi)],
      colnames(f.lumi)[4:ncol(f.lumi)],
      colnames(g.lumi)[4:ncol(g.lumi)],
      colnames(h.lumi)[4:ncol(h.lumi)],
      colnames(i.lumi)[4:ncol(i.lumi)],
      colnames(j.lumi)[4:ncol(j.lumi)])
    )
    
    a.lumi <- a.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    b.lumi <- b.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    c.lumi <- c.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    d.lumi <- d.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    e.lumi <- e.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    f.lumi <- f.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    g.lumi <- g.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    h.lumi <- h.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    i.lumi <- i.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    j.lumi <- j.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    
    lumi <- rbind(a.lumi, b.lumi, c.lumi, d.lumi,
                  e.lumi, f.lumi, g.lumi, h.lumi,
                  i.lumi, j.lumi)
    
  } else {
    
    feature_intersection <- Reduce(intersect, list(
      colnames(b.lumi)[4:ncol(b.lumi)],
      colnames(c.lumi)[4:ncol(c.lumi)],
      colnames(d.lumi)[4:ncol(d.lumi)],
      colnames(e.lumi)[4:ncol(e.lumi)],
      colnames(f.lumi)[4:ncol(f.lumi)],
      colnames(g.lumi)[4:ncol(g.lumi)],
      colnames(h.lumi)[4:ncol(h.lumi)],
      colnames(i.lumi)[4:ncol(i.lumi)],
      colnames(j.lumi)[4:ncol(j.lumi)])
    )
    
    b.lumi <- b.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    c.lumi <- c.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    d.lumi <- d.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    e.lumi <- e.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    f.lumi <- f.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    g.lumi <- g.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    h.lumi <- h.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    i.lumi <- i.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    j.lumi <- j.lumi[, c('array_barcode', 'model_id', 'passage', 'project_title', feature_intersection)]
    
    lumi <- rbind(b.lumi, c.lumi, d.lumi,
                  e.lumi, f.lumi, g.lumi, 
                  h.lumi, i.lumi, j.lumi)
  }
  
  return(lumi)
  
}

##########
# function that corrects for batch effects of full data
# ##########
# full_data <- lumi_correct
# full = T
# a_lumi <- a.lumi_fin
getCombat <- function(full_data, full, a_lumi = NULL)
{
  # remove rows with NA in model_id
  full_data <- full_data[!is.na(full_data$model_id),]
  
  if(full) {
    feature_intersection <- Reduce(intersect, list(
      colnames(a_lumi)[5:ncol(a_lumi)], 
      colnames(full_data)[2:ncol(full_data)]))
    full_data <- full_data[, c('model_id', 'batch', feature_intersection)]
    a_lumi <- a_lumi[, c('model_id', 'project_title', feature_intersection)]
    names(full_data)[2] <-'project_title'
    
    # remove columns from a.lumi
    a_lumi$array_barcode <- NULL
    a_lumi$passage <- NULL
    
    full_data <- rbind(a_lumi, full_data)
    # recode project title to two batches 
    full_data$project_title <- ifelse(grepl('Tsao', full_data$project_title), 'a', full_data$project_title)
    
    full_data <- full_data[complete.cases(full_data),]
  }
  
  
  # remove unnecessary columns
  full_data$array_barcode <- full_data$passage <- full_data$project_title.1 <- NULL
  
  # first make numeric data numeric
  full_data[, 3:ncol(full_data)] <- apply(full_data[, 3:ncol(full_data)], 2, as.numeric)
  # add unique indicator to duplicated model ids - multiple duplicates - this will do for now
  dup_index <- duplicated(full_data$model_id)
  full_data$model_id[dup_index] <- paste0(full_data$model_id[dup_index], '_dup')
  
  dup_index <- duplicated(full_data$model_id)
  full_data$model_id[dup_index] <- paste0(full_data$model_id[dup_index], '_dup2')
  
  dup_index <- duplicated(full_data$model_id)
  full_data$model_id[dup_index] <- paste0(full_data$model_id[dup_index], '_dup3')
  
  dup_index <- duplicated(full_data$model_id)
  full_data$model_id[dup_index] <- paste0(full_data$model_id[dup_index], '_dup4')
  
  dup_index <- duplicated(full_data$model_id)
  full_data$model_id[dup_index] <- paste0(full_data$model_id[dup_index], '_dup5')
  
  # put model ids in rownames and remove columns
  rownames(full_data) <- full_data$model_id
  mat_data <- full_data[, 3:ncol(full_data)]
  # get features 
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # make full just two batches
  if(full) {
    full_data$project_title <- ifelse(full_data$project_title == 'a', 'a', 'b')
  }
  
  # get intercept
  modcombat <- model.matrix(~1, data = full_data)
  batch <- as.factor(full_data$project_title)
  combat <- ComBat(dat = mat_data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$model_id <- rownames(final_dat)
  final_dat$batch <- batch
  final_dat <- final_dat[, c('model_id', 'batch', features)]
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
}


##########
# function that corrects for batch effects of full data
# ##########
# full_data <- lumi_correct
# full = T
# a_lumi <- a.lumi_fin
getCombatNov <- function(data, full)
{
  
  # add unique indicator to duplicated model ids - multiple duplicates - this will do for now
  dup_index <- duplicated(data$id)
  data$id[dup_index] <- paste0(data$id[dup_index], '_dup')
  
  dup_index <- duplicated(data$id)
  data$id[dup_index] <- paste0(data$id[dup_index], '_dup2')
  
  dup_index <- duplicated(data$id)
  data$id[dup_index] <- paste0(data$id[dup_index], '_dup3')
  
  dup_index <- duplicated(data$id)
  data$id[dup_index] <- paste0(data$id[dup_index], '_dup4')
  
  dup_index <- duplicated(data$id)
  data$id[dup_index] <- paste0(data$id[dup_index], '_dup5')
  
  # put model ids in rownames and remove columns
  rownames(data) <- data$id
  mat_data <- data[, 3:ncol(data)]
  # get features 
  features <- colnames(mat_data)
  mat_data <- t(mat_data)
  
  # # full
  # if(full) {
  #   data$batch <- ifelse(data$)
  # }
  
  # get intercept
  modcombat <- model.matrix(~1, data = data)
  batch <- as.factor(data$batch)
  combat <- ComBat(dat = mat_data, batch = batch, mod = modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  # transpose and add back columns
  final_dat <- as.data.frame(t(combat))
  final_dat$id <- rownames(final_dat)
  final_dat$batch <- batch
  final_dat <- final_dat[, c('id', 'batch', features)]
  rownames(final_dat) <- NULL
  
  return(final_dat)
  
}

##########
# function to prepare data for snf
##########
prepareSNF <- function(data_new) 
{
  # sort the model_ids 
  data_new <- data_new[order(data_new$model_id, decreasing = F),]
  
  # put model id in rownames
  rownames(data_new) <- data_new$model_id
  
  # subet by 3rd column
  data_new <- data_new[, 3:ncol(data_new)]
  
  # transpose
  data_final <- t(data_new)
  
  return(data_final)
}

##########
# get snf for lung object
#########
# data <- quan
# sorted_ids <- quan.sorted_id
# sampleRows <- F
# numClust <- 5
lungSNf <- function(data, numClust, sampleRows, sorted_ids)
{
  # Transpose data for the dist function
  # The dist function takes distances between rows
  sampleRowsRequired <- TRUE
  transposeData <- sampleRows != sampleRowsRequired
  if (transposeData) {
    data <- t(data)
  }
  
  # Calculate the distance between samples
  distances <- as.matrix(dist(data))
  
  # Convert the distances to affinities
  affinities <- affinityMatrix(distances)
  
  # Cluster the fused matrix
  labels <- spectralClustering(affinities, numClust)
  
  # get id for labels
  labels <- as.data.frame(cbind(labels, id = sorted_ids))
  
  return(labels)
  
}

##########
# PCA of data 
##########

getPCA <- function(pca_data, column_name, gene_start, name, full) 
{
  # get features sites
  genes <- colnames(pca_data)[gene_start:ncol(pca_data)]
  
  # if (full){
  #   pca_data <- pca_data[!is.na(pca_data$model_id),]
  # }
  
  # put column name with genes
  pca_data <- pca_data[ ,c(column_name, genes)]
  
  # make data numeric 
  pca <- prcomp(pca_data[,2:ncol(pca_data)])
  
  # plot data
  #fill in factors with colors 
  if (full) {
    pca_data[, column_name] <- ifelse(pca_data[, column_name] == 'a', 'a', 'b')
    col_vec <- c('red', 'green')
  } else {
    col_vec <- c('red', 'green', 'blue', 'brown', 'black', 
                 'orange', 'yellow', 'grey', 'lightblue', 'lightgreen')
  }
  
  colors <- col_vec[as.factor(pca_data[, column_name])]

  
  plot <- plot(pca$x[,1], 
               pca$x[,2],
               xlab = 'pca 1',
               ylab = 'pca 2',
               cex = 1,
               main = name,
               pch = 16,
               col = adjustcolor(colors, alpha.f = 0.5)
  )
  abline(v = c(0,0),
         h = c(0,0))
  
}



plotResult <- function(drug_result, main) 
{
  x_max <- max(unlist(drug_result[[4]]))
  y_max <- max(unlist(drug_result[[6]]))
  
  # plot predictions against ground truth
  plot(unlist(drug_result[[4]]), unlist(drug_result[[6]]), 
       xlim = c(0, x_max),
       ylim = c(0, y_max),
       bty = 'n',
       col = adjustcolor('blue', alpha.f = 0.6),
       pch = 16,
       xlab = 'Predictions',
       ylab = 'Real Values',
       main = main)
  abline(0,1, lty = 3)
  corr <- round(cor(unlist(drug_result[[4]]), unlist(drug_result[[6]])), 2)
  legend("topleft", legend = paste0('correlation = ', corr), cex = 1, bty = 'n')
  legend("bottomright", legend = paste0('# samples = ', 
                                        drug_result[[11]][1]))
  
  
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





