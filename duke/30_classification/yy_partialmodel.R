#!/usr/bin/env Rscript --vanilla

# System Arguments
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)
nthreads <- 24

#' ## Setup
#' 
#' This has packages, paths, and other things to setup for running the code.
#+ setup
cat("Setup\n")
# packages
library(plyr)
library(ggplot2)
library(RColorBrewer)
library(caret)
library(glmnet)
library(Matrix)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(nthreads)
# file paths
base       <- "/data1/ffg05"
roidir     <- file.path(base, "analysis/rois")
scriptdir  <- file.path(base, "scripts/repsim")
taskdir    <- file.path(base, "analysis/task_activity")
timedir    <- file.path(base, "notes")
# variable information
subjects   <- as.character(read.table(file.path(scriptdir, "sublist.txt"))[,1])
# rois
roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", "R")
roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", "M1", "M1", "A1", "A1")
roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
                            fullnames=paste(roi.hemis, roi.names))
roi.df        <- roi.df[-1,]
roi.names.title <- roi.df$fullnames
roi.names     <- sub(" ", "_", tolower(roi.names.title))
# categories
categories <- c("faces", "fruits", "letters", "vehicles")
categories.title <- c("Faces", "Fruits", "Letters", "Vehicles")

cat("Source\n")
wrap.model <- NULL; save.vars <- NULL
source("/data1/ffg05/scripts/distrep/duke/30_classification/yy_partialmodel_funs.R")

#' ## Load Data
#' 
#' We get the brain activity (beta-series) for each category.
#' This is done for each subject and combined together into a list.
#' 
#+ load-data
cat("Load\n")
lst.dats <- llply(subjects, load.data, .progress="text")
names(lst.dats) <- subjects

#' Load the processed data to get the best regions
#' 
#' We want to order our regions based on how often they were kept in the final
#' reduced model across subjects. We will then add each region to our model 
#' based on how consistently it was found across subjects.
#' 
#+ load-results
rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")
rns <- llply(rets, function(x) {
  as.character(roi.df$fullnames[roi.df$vals %in% x$res$alpha0.5$partial$rois])
})
print(sort(table(unlist(rns)), decreasing=T))
sorted.roi.names <- names(sort(table(unlist(rns)), decreasing=T))




#' Sample run
# loop through rois
# Get the roi name and roi number that we will be dealing with
# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
all.rns <- sroi.names
cur.urois <- c()
mat.accs <- matrix(0, length(subjects), length(roi.names))
best.acc <- 0

# run each subject with full model
cat("Full Model\n")
prev.accuracies <- laply(lst.dats, function(dat) {
  sdat    <- scale(dat$dat)
  rois    <- dat$grouping
  yfactor <- dat$yfactor
  runs    <- dat$runs
  
  # Combine the ROIs that are bilateral into one ROI
  #old.rois <- rois 
  tab <- table(roi.df$names)
  for (x in names(tab[tab>1])) {
    inds <- which(roi.df$names==x)
    rois[rois %in% inds] <- inds[1] # set everything to the first one
  }
  # get the unique numbers
  urois <- sort(unique(rois))
  
  # Get the cross-folds
  nrepeats <- 1
  uruns    <- sort(unique(runs))
  nruns    <- length(uruns)
  runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
  foldsI   <- lapply(runFolds, function(x) {
    runs %in% x
  })
  
  # Model
  res <- wrap.model(sdat, yfactor, foldsI, alpha=0.5, nlambda=100)
  res <- save.vars(res, rois)
  
  return(res$acc)
}, .progress="text")

cat("Reduced Model\n")
cur.urois <- which(!duplicated(roi.df$names))
for (iter in 1:(length(cur.urois)-1)) {
#for (iter in 2:(length(cur.urois)-1)) {
  cat("+++\nIteration:", iter, "\n")
  cat("fitting\n")
  
  # Run each possible model with one of the regions removed
  roi.combns <- llply(1:length(cur.urois), function(i) {
    cur.urois[-i]
  })
  
  # look through subjects
  cur.accuracies <- laply(lst.dats, function(dat) {
    sdat    <- scale(dat$dat)
    rois    <- dat$grouping
    yfactor <- dat$yfactor
    runs    <- dat$runs
    
    # Combine the ROIs that are bilateral into one ROI
    #old.rois <- rois 
    tab <- table(roi.df$names)
    for (x in names(tab[tab>1])) {
      inds <- which(roi.df$names==x)
      rois[rois %in% inds] <- inds[1] # set everything to the first one
    }
    # get the unique numbers
    urois <- sort(unique(rois))
    
    # Get the cross-folds
    nrepeats <- 1
    uruns    <- sort(unique(runs))
    nruns    <- length(uruns)
    runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
    foldsI   <- lapply(runFolds, function(x) {
      runs %in% x
    })
    
    # Model
    res <- laply(roi.combns, function(reduced.urois) {
      x <- sdat[, rois %in% reduced.urois] # take subset of data
      res <- wrap.model(x, yfactor, foldsI, alpha=0.5, nlambda=100)
      res <- save.vars(res, rois[rois %in% reduced.urois])
      res$acc
    })
    
    # Return the results
    return(res)
  }, .parallel=FALSE, .progress="text")
  
  # Determine the model with the largest positive or least negative change 
  # in accuracy relative to the best model accuracy in the previous iteration
  cat("determining best\n")
  mean.prev.acc <- mean(prev.accuracies)
  mean.cur.accs <- colMeans(cur.accuracies)
  best.ind      <- which.max(mean.cur.accs - mean.prev.acc)
  best.diff     <- max(mean.cur.accs - mean.prev.acc)
  
  # If there's no improvement (<0) in accuracy with the best model, then end
  if (best.diff < 0) {
    cat("nothing to remove, stopping!!!\n")
    break
  } else {
    disp <- roi.names.title[roi.combns[[best.ind]]]
    disp <- paste(disp, collapse=', ')
    cat('removed:', as.character(roi.names.title)[cur.urois[best.ind]], "\n")
    cat("new model is:", disp, "\n")
  }
  
  # Save the best accuracies and rois
  prev.accuracies <- cur.accuracies[,best.ind]
  cur.urois       <- roi.combns[[best.ind]]
}

# The result here is a reduced model including the FFA, PPA, LOC, VWF, and V1

## redo roi.df
#roi.df <- roi.df[!duplicated(roi.df$names),]
#roi.df$fullnames <- sroi.names

save()
