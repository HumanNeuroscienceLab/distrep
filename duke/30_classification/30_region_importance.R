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
categories <- c("faces", "letter", "fruits", "vehicles")
categories.title <- c("Faces", "Letters", "Fruits", "Vehicles")
# new roi names
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")

outdir <- "/data1/ffg05/analysis/duke/classification"

cat("Source\n")
wrap.model <- NULL; save.vars <- NULL
source("/data1/ffg05/scripts/distrep/duke/30_classification/yy_partialmodel_funs.R")


###
# Determine importance by number of times region is found in reduced model
###

#' Load the processed data to get the best regions
#' 
#' We want to order our regions based on how often they were kept in the final
#' reduced model across subjects. We will then add each region to our model 
#' based on how consistently it was found across subjects.
#' 
#+ load-results
# Get all the data
subs.rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")
# Collect roi names selected by each subject with the elastic-net
rns <- llply(subs.rets, function(x) {
  as.character(roi.df$fullnames[roi.df$vals %in% x$res$alpha0.5$partial$rois])
})
# Print out those names as a table
print(table(unlist(rns)))
print(sort(table(unlist(rns)), decreasing=T))
# Save
outdir <- "/data1/ffg05/analysis/duke/classification"
outfile <- file.path(outdir, "region_importance_reduced_model.csv")
saveret <- sort(table(unlist(rns)), decreasing=T)
names(saveret) <- sroi.names
write.csv(t(saveret), row.names=F, file=outfile)


###
# Determine importance by change in accuracy by leave-one-region-out
###

#' ## Load Data
#' 
#' We get the brain activity (beta-series) for each category.
#' This is done for each subject and combined together into a list.
#' 
#+ load-data
cat("Load\n")
lst.dats <- llply(subjects, load.data, .progress="text")
names(lst.dats) <- subjects


# Let's first get the accuracies of the full model before leaving any region out
cat("Full Model\n")
full.accuracies <- laply(lst.dats, function(dat) {
  sdat    <- scale(dat$dat)
  rois    <- dat$grouping
  yfactor <- dat$yfactor
  runs    <- dat$runs
  
  # Combine the ROIs that are bilateral into one ROI
  tab <- table(roi.df$names)
  for (x in names(tab[tab>1])) {
    inds <- which(roi.df$names==x)
    rois[rois %in% inds] <- inds[1] # set everything to the first one
  }
  
  # Get the cross-folds
  nrepeats <- 1
  uruns    <- sort(unique(runs))
  nruns    <- length(uruns)
  runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
  foldsI   <- lapply(runFolds, function(x) {
    runs %in% x
  })
  
  # Model
  res <- wrap.model(sdat, yfactor, foldsI, alpha=1, nlambda=100)
  res <- save.vars(res, rois)
  
  return(res$acc)
}, .progress="text")


# Now we leave each region out in each subject and see the drop/increase in 
# classification accuracy
cat("Leave One Region Out Model\n")

# Run each possible model with one of the regions removed
cur.urois <- which(!duplicated(roi.df$names))
roi.combns <- llply(1:length(cur.urois), function(i) {
  cur.urois[-i]
})

# look through subjects
reduced.accuracies <- laply(lst.dats, function(dat) {
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
    res <- wrap.model(x, yfactor, foldsI, alpha=1, nlambda=100)
    res <- save.vars(res, rois[rois %in% reduced.urois])
    res$acc
  })
  
  # Return the results
  return(res)
}, .parallel=FALSE, .progress="text")

# Combine the accuracies
colnames(reduced.accuracies) <- sroi.names
accuracies <- cbind(Full=full.accuracies, reduced.accuracies)
head(accuracies)

# Save
outfile <- file.path(outdir, "region_importance_leave_one_out.csv")
write.csv(accuracies, row.names=F, col.names=T, quote=F, file=outfile)




# Get the whole classification accurayc


# TODO: permutation test

head(tmp2)

###
# Determine the importance of a region based on it's individual accruacy
###

#' ## Load Data
#' 
#' Load the previously run classification analyses from each subject.
#' 
#+ load-data
sub.rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")

#' ## Combine Accuracies
#' 
#' We are only looking at the lasso results. Note that we also combine the 
#' pvalues.
#' 
#+ combine-individual-accs

indiv.df <- ldply(1:length(sub.rets), function(si) {
  indiv.res <- sub.rets[[si]]$res$alpha1$individual
  data.frame(
    subject=si, 
    roi=sroi.names, 
    accuracy=indiv.res$acc, 
    pvalue=indiv.res$pval
  )
})
head(indiv.df)

#' ## Save
#' 
#+ individual-save
outfile <- file.path(outdir, "region_importance_individual_accuracy.csv")
write.csv(indiv.df, row.names=F, col.names=T, quote=T, file=outfile)
