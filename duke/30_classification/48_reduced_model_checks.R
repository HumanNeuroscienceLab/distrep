
###
# Group Average
###

#' This page will deal with getting the output from the group-level analyses
#' in particular with a reduced 4 or 5 region model.
#' 
#' ## Setup
#' 
#' First, we load all the packages
#' 
#+ load-packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)
library(plyr)
library(RColorBrewer)
library(corrplot)
library(ggplot2)
library(ggthemes)
library(caret)
library(glmnet)
library(Matrix)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)


#' Let's load some needed functions including
#' 
#' - `load.data`: Does what it says, loads the data for one subject
#' - `get_stats`: Used to get classification performance
#' - `category_stats`: Used to get classification performance for each category
#' 
#+ load-funs
load.data.std <- NULL; wrap.model <- NULL
source("/data1/ffg05/scripts/repsim/scratch/stats_lab_funs.R")

#' Now we gather useful variables like the list of categories for the 
#' different stimuli and the subject
#' 
#+ vars
roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", "R")
roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", "M1", "M1", "A1", "A1")
roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
                            fullnames=paste(roi.hemis, roi.names))
roi.df        <- roi.df[-1,]
roi.names.title <- roi.df$fullnames
roi.names     <- sub(" ", "_", tolower(roi.names.title))
# categories
categories <- c("faces", "letters", 'fruits', "vehicles")
categories.title <- c("Faces", "Letters", "Fruits", "Vehicles")
# cols
dat.cols <- rev(brewer.pal(11, "RdYlBu"))
cat.cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # Faces, Letters, Fruits, & Vehicles
#roi.cols <- brewer.pal(8, "Set3")[-c(1,7)]
roi.cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)]
# subjects
subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])
## output for saving
outdir <- "/data1/ffg05/analysis/duke/classification"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/duke/classification/plots"
dir.create(plotdir, showWarnings = FALSE)

#' Let's get the rois, grouping, etc information
#' 
#+ vars2 
tmp        <- load.data.std(subjects[1])
mask       <- tmp$orig.rois
grouping   <- tmp$grouping
rois       <- tmp$rois
ys         <- tmp$yfactor
categories <- tmp$categories # this overwrites the old categories variable
roi.df$nvoxs <- table(grouping)


#' ## Load Data
#'  
#' Load in all the subject data
#' 
#+ data-load
dat.subs.std <- laply(subjects, function(subject) {
  cat(subject, "\n")
  sublst       <- load.data.std(subject)
  scale(sublst$dat)
})


#' Combine the ROIs that are bilateral into one ROI
#' 
#+ data-select
# create a new roi grouping vector
old.grouping <- grouping 
old.mask     <- mask
mask[mask==10] <- 0; mask[mask==12] <- 1; mask[mask==40] <- 6
tab <- table(roi.df$names)
for (x in names(tab[tab>1])) {
  inds <- which(roi.df$names==x)
  grouping[grouping %in% inds] <- inds[1] # set everything to the first one
  vals <- roi.df$vals[inds]
  mask[mask %in% vals] <- inds[1]
}
# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
# get the unique numbers
urois <- sort(unique(grouping))

#' Average the data across subjects
#' 
#+ data-average
grpdat <- apply(dat.subs.std, c(2,3), mean, na.rm=T)
# weird NA for 6th subject at voxel 1074...exclude from average


#' ## Classification
#' 
#' ### Cross-Validation Setup
#' 
#' We will be doing a leave-one run out cross-validation. Since we have 10 runs
#' of data for each subject, this could also be considered a 10-fold cross-
#' validation. Presentation of each stimulus from each category is more or less
#' evenly distributed across the runs. Here we will get the indices for each 
#' cross-validation fold.
#' 
#+ cross-validation-setup
nrepeats <- 1
runs     <- tmp$runs # since all subjects the same can use runs from 1 subject
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
folds    <- lapply(runFolds, function(x) {
  which(runs %in% x)
})
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})


#' ### Run the Four Models At The Group Level
#' 
#' Now we run the ridge, elastic-net, and lasso on the group model. Here we run
#' the full, partial, and individual region models...just cuz.
#' 
#+ run-models
sgrpdat         <- scale(grpdat)
subset.urois    <- urois[1:4]
subset.grouping <- grouping[grouping %in% subset.urois]
subset.grpdat   <- sgrpdat[,grouping %in% subset.urois]
subset.roi.orders <- roi.orders[grouping %in% subset.urois]
alpha <- 1

#rgrp.res4 <- wrap.rcfe(subset.grpdat, ys, foldsI, subset.grouping, 
#                       alpha=alpha, nlambda=100)
grp.res4 <- wrap.model(subset.grpdat, ys, foldsI, alpha=alpha, nlambda=100)

#' ### Beta Coefficient
#' 
#' Let's get the beta-coefficients
#' 
#+ std-betas
sublst       <- load.data.std(subjects[2])
mean.betas <- laply(1:20, function(si) {
  inds <- grouping < 7
  tmp <- apply(dat.subs.std[si,,inds], 2, function(x) tapply(x, sublst$yfactor, mean))
  tmp2 <- apply(tmp, 1, tapply, grouping[inds], mean)
  tmp2
})
mean.betas <- apply(mean.betas, 2:3, mean)
rownames(mean.betas) <- sroi.names[1:4]
round(mean.betas, 2)
# for group data
tmp <- apply(subset.grpdat, 2, function(x) tapply(x, ys, mean))
tmp2 <- apply(tmp, 1, tapply, subset.grouping, mean)
rownames(tmp2) <- sroi.names[1:4]
round(tmp2, 2)


#' ### Classification
#' 
#' Let's try the classification again here
#' 
#+ run-classifiers
lst.sres <- llply(1:length(subjects), function(si) {
  #cat(si, "\n")
  
  # Loop through each alpha
  inds     <- grouping < 7
  sdat     <- dat.subs.std[si,,inds]
  obs      <- sublst$yfactor
  
  # Run the model!!!
  res      <- wrap.model(sdat, obs, foldsI, alpha=1, nlambda=100, parallel=T)
  
  res
}, .parallel=F, .progress="text")

#' Summarize the location of non-zero voxels
#' 
#+ 
sret <- laply(1:20, function(si) {
  sres      <- lst.sres[[si]]
  
  inds     <- grouping < 7
  grp      <- grouping[inds]
  betas    <- sres$best$betas[,inds]
  
  fac      <- factor(grp, levels=c(1,2,4,6), labels=c("FFA", "PPA", "LOC", "VWFA"))
  
  ret <- laply(1:4, function(fi) {
    cbind(pos=table(fac[betas[fi,]>0]), neg=table(fac[betas[fi,]<0]))
  })
  dimnames(ret)[[1]] <- categories.title
  
  ret
})
round(apply(sret, 2:4, mean), 1)



###
# Subject
###

# Load the data
dat.subs <- NULL
infile <- file.path(outdir, "subjects_4region_models.rda")
# ran this: save(lst.res, dat.subs, roi.df, categories.title, file=outfile)
load(infile)

#' ### Beta Coefficient
#' 
#' Let's do it here for the native space data!
#' 
#+ native-betas
mean.betas.nat <- laply(1:20, function(si) {
  tmp  <- apply(dat.subs[[si]]$dat, 2, function(x) tapply(x, dat.subs[[si]]$yfactor, mean))
  tmp2 <- apply(tmp, 1, tapply, dat.subs[[si]]$grouping, mean)
  tmp2
})
mean.betas.nat <- apply(mean.betas.nat, 2:3, mean)
rownames(mean.betas.nat) <- sroi.names[1:4]
round(mean.betas.nat, 2)
# Scaled
mean.betas.nat2 <- laply(1:20, function(si) {
  tmp  <- apply(scale(dat.subs[[si]]$dat), 2, function(x) tapply(x, dat.subs[[si]]$yfactor, mean))
  tmp2 <- apply(tmp, 1, tapply, dat.subs[[si]]$grouping, mean)
  tmp2
})
mean.betas.nat2 <- apply(mean.betas.nat2, 2:3, mean)
rownames(mean.betas.nat2) <- sroi.names[1:4]
round(mean.betas.nat2, 2)

#' ### Classification
#' 
#' I'm even re-running the classification.
#' 
#+ native-classification
lst.res2 <- llply(1:length(subjects), function(si) {
  #cat(si, "\n")
  
  # Loop through each alpha
  sdat     <- dat.subs[[si]]$dat
  obs      <- dat.subs[[si]]$yfactor
  
  # Run the model!!!
  res      <- wrap.model(sdat, obs, foldsI, alpha=1, nlambda=100, parallel=T)
  
  res
}, .parallel=F, .progress="text")
nret <- laply(1:20, function(si) {
  res      <- lst.res2[[si]]
  betas    <- res$best$betas
  fac      <- factor(dat.subs[[si]]$grouping, levels=c(1:4), labels=c("FFA", "PPA", "LOC", "VWFA"))
  
  ret <- laply(1:4, function(fi) {
    cbind(pos=table(fac[betas[fi,]>0]), neg=table(fac[betas[fi,]<0]))
  })
  dimnames(ret)[[1]] <- categories.title
  
  ret
})
round(apply(nret, 2:4, mean), 1)

