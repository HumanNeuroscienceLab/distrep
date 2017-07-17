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
# FOR EACH CATEGORY:
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

#' ## Analysis
#' 
#' ### Full Model
#' Let's first get the accuracies of the full model before leaving any region out
#' 
#+ full-model-analysis
cat("Full Model\n")
lst.full <- llply(lst.dats, function(dat) {
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
  res <- wrap.model(sdat, yfactor, foldsI, alpha=1, nlambda=100, parallel=T)
  res <- save.vars(res, rois)
  
  return(res)
}, .parallel=F, .progress="text")
names(lst.full) <- subjects

#' ### Reduced Model
#' 
#' Remove one region at a time.
#' 
#+ reduced-model-analysis
cur.urois <- which(!duplicated(roi.df$names))
roi.combns <- llply(1:length(cur.urois), function(i) {
  cur.urois[-i]
})
lst.reduced <- llply(lst.dats, function(dat) {
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
  lst.res <- llply(roi.combns, function(reduced.urois) {
    x <- sdat[, rois %in% reduced.urois] # take subset of data
    res <- wrap.model(x, yfactor, foldsI, alpha=1, nlambda=100, parallel=T)
    res <- save.vars(res, rois[rois %in% reduced.urois])
    res
  }, .parallel=T)
  names(lst.res) <- sub("[RL]\ ", "", sroi.names)
  
  # Return the results
  return(lst.res)
}, .parallel=F, .progress="text")
names(lst.reduced) <- subjects

#' ### Individual Regions
#' 
#+ individual-model-analysis
lst.individual <- llply(1:length(lst.dats), function(i) {
  cat("subject", i, "\n")
  dat     <- lst.dats[[i]]
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
  lst.res <- llply(urois, function(uroi) {
    # take subset of data
    x   <- sdat[, rois == uroi]
    res <- wrap.model(x, yfactor, foldsI, alpha=1, nlambda=100, parallel=F)
    res <- save.vars(res, rois[rois == uroi])
    res
  }, .parallel=T)
  names(lst.res) <- sub("[RL]\ ", "", sroi.names[1:4])
  
  # Return the results
  return(lst.res)
}, .parallel=F, .progress="text")
names(lst.individual) <- subjects


#' ### Category Accuracies
#' 
#' We find the best accuracies in each category.
#' 
#+ category-accuracies
# TODO: make sure average accuarcy works and stuff
get.best.accs <- function(res) {
  ave.acc  <- as.numeric(as.character(res$best$stats$Accuracy))
  cat.accs <- t(sapply(1:100, function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    acc <- category_stats(res$obs, pred, prob)$Accuracy
    as.numeric(as.character(acc))
  }))
  best.accs <- apply(cat.accs, 2, max)
  names(best.accs) <- categories.title
  best.accs <- c(average=ave.acc, best.accs)
  return(best.accs)
}
full.accs    <- laply(lst.full, function(x) get.best.accs(x$res), .parallel=T)
reduced.accs <- laply(1:7, function(ri) {
  laply(lst.reduced, function(x) get.best.accs(x[[ri]]$res), .parallel=T)
}, .progress="text")
dimnames(reduced.accs)[[1]] <- sub("[RL]\ ", "", sroi.names) # Remove R and L 
# subtract each of the fulls from the reduced
titles <- c("Average", categories.title)
diff.accs <- aaply(reduced.accs, .(1), function(x) round(full.accs*100 - x*100, 1))
dimnames(diff.accs) <- list(
  roi=names(reduced.accs), subject=1:nrow(full.accs), category=titles
)
mean.diff.accs <- apply(diff.accs, c(1,3), mean)
mean.diff.accs

# Save
outfile <- file.path(outdir, "category_accuracies_leave_region_out.rda")
save(full.accs, reduced.accs, diff.accs, file=outfile)



#' ## Individual Region Category Accuracies
#' 
#' Load the previously run classification analyses from each subject.
#' 
#+ load-data
sub.rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")

#' ### Combine Accuracies
#' 
#' We are only looking at the lasso results. Note that we also combine the 
#' pvalues.
#' 
#+ combine-individual-accs

#laply(1:length(sub.rets), function(i) {
#  laply(1:7, function(j) {
#    dim(sub.rets[[1]]$res$alpha1$individual$res[[1]]$res$probs)[3]
#  })
#})
get.best.accs <- function(res) {
  ave.acc  <- as.numeric(as.character(res$best$stats$Accuracy))
  cat.accs <- t(sapply(1:length(res$preds), function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    acc <- category_stats(res$obs, pred, prob)$Accuracy
    as.numeric(as.character(acc))
  }))
  best.accs <- apply(cat.accs, 2, max)
  names(best.accs) <- categories.title
  best.accs <- c(Average=ave.acc, best.accs)
  return(best.accs)
}
indiv.df <- ldply(1:length(sub.rets), function(si) {
  #cat(si, "\n")
  indiv.res <- sub.rets[[si]]$res$alpha1$individual$res
  ldply(1:length(indiv.res), function(ri) {
    #cat(ri,".")
    res  <- indiv.res[[ri]]$res
    accs <- get.best.accs(res)
    data.frame(
      subject=si, 
      roi=sroi.names[ri], 
      category=names(accs), 
      accuracy=accs
    )
  }, .parallel=T)
}, .progress="text")
head(indiv.df)

# Save
outfile <- file.path(outdir, "category_accuracies_individual_region.csv")
save(indiv.df, file=outfile)





#' # Standard Space
#' 
#' I want to look at the region importance using those standard space ROIs
#+ stuff

# TODO: add setup

dat.subs.std <- laply(subjects, function(subject) {
  cat(subject, "\n")
  sublst       <- load.data.std(subject)
  scale(sublst$dat)
})

tmp        <- load.data.std(subjects[1])
grouping   <- tmp$grouping
rois       <- tmp$rois
ys         <- tmp$yfactor
categories <- tmp$categories # this overwrites the old categories variable
roi.df$nvoxs <- table(grouping)

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

#' Combine the ROIs that are bilateral into one ROI
#' 
#+ data-select
# create a new roi grouping vector
old.grouping <- grouping
tab <- table(roi.df$names)
for (x in names(tab[tab>1])) {
  inds <- which(roi.df$names==x)
  grouping[grouping %in% inds] <- inds[1] # set everything to the first one
}
# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
# get the unique numbers
urois <- sort(unique(grouping))


#' ### Full Model
#' 
#+ full-model-analysis
# first fix NAs in 6th subject
dat.subs.std[6,,1074] <- rowMeans(dat.subs.std[6,,-1074])
lst.sfull <- llply(1:20, function(si) {
  sdat    <- dat.subs.std[si,,]
  rois    <- grouping
  yfactor <- ys
  runs    <- runs
  
  #inds    <- rois < 7 # to keep only the four regions
  #rois    <- rois[inds]
  #sdat    <- sdat[,inds]
  
  # Get the cross-folds
  nrepeats <- 1
  uruns    <- sort(unique(runs))
  nruns    <- length(uruns)
  runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
  foldsI   <- lapply(runFolds, function(x) {
    runs %in% x
  })
  
  # Model
  res <- wrap.model(sdat, yfactor, foldsI, alpha=1, nlambda=100, parallel=T)
  res <- save.vars(res, rois)
  
  # Return the results
  return(res)
}, .parallel=F, .progress="text")
names(lst.sfull) <- subjects

#' ### Reduced Model
#' 
#' Remove one region at a time.
#' 
#+ reduced-model-analysis
cur.urois <- unique(grouping)
roi.combns <- llply(1:length(cur.urois), function(i) {
  cur.urois[-i]
})
lst.sreduced <- llply(1:20, function(si) {
  sdat    <- dat.subs.std[si,,]
  rois    <- grouping
  yfactor <- ys
  runs    <- runs
  
  #inds    <- rois < 7 # to keep only the four regions
  #rois    <- rois[inds]
  #sdat    <- sdat[,inds]
  
  # Get the cross-folds
  nrepeats <- 1
  uruns    <- sort(unique(runs))
  nruns    <- length(uruns)
  runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
  foldsI   <- lapply(runFolds, function(x) {
    runs %in% x
  })
  
  # Model
  lst.res <- llply(roi.combns, function(reduced.urois) {
    x <- sdat[, rois %in% reduced.urois] # take subset of data
    res <- wrap.model(x, yfactor, foldsI, alpha=1, nlambda=100, parallel=T)
    res <- save.vars(res, rois[rois %in% reduced.urois])
    res
  }, .parallel=F)
  names(lst.res) <- sub("[RL]\ ", "", sroi.names[1:4])
  
  # Return the results
  return(lst.res)
}, .parallel=F, .progress="text")
names(lst.sreduced) <- subjects

#' ### Category Accuracies
#' 
#' We find the best accuracies in each category.
#' 
#+ category-accuracies2
# TODO: make sure average accuarcy works and stuff
get.best.accs <- function(res) {
  ave.acc  <- as.numeric(as.character(res$best$stats$Accuracy))
  cat.accs <- t(sapply(1:100, function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    acc <- category_stats(res$obs, pred, prob)$Accuracy
    as.numeric(as.character(acc))
  }))
  best.accs <- apply(cat.accs, 2, max)
  names(best.accs) <- categories.title
  best.accs <- c(average=ave.acc, best.accs)
  return(best.accs)
}
full.accs    <- laply(lst.sfull, function(x) get.best.accs(x$res), .parallel=T)
reduced.accs <- laply(1:7, function(ri) {
  laply(lst.sreduced, function(x) get.best.accs(x[[ri]]$res), .parallel=T)
}, .progress="text")
dimnames(reduced.accs)[[1]] <- sub("[RL]\ ", "", sroi.names) # Remove R and L 

# subtract each of the fulls from the reduced
titles <- c("Average", categories.title)
diff.accs <- aaply(reduced.accs, .(1), function(x) full.accs*100 - x*100)
dimnames(diff.accs) <- list(
  roi=sroi.names, subject=1:nrow(full.accs), category=titles
)
mean.diff.accs <- apply(diff.accs, c(1,3), mean)
mean.diff.accs

sd.diff.accs <- apply(diff.accs, c(1,3), sd)

tval.diff.accs  <- apply(diff.accs, c(1,3), function(x) t.test(x)$statistic)
pval.diff.accs  <- apply(diff.accs, c(1,3), function(x) t.test(x)$p.value)
fdr.diff.accs   <- matrix(p.adjust(pval.diff.accs, "fdr"), 
                          nrow(pval.diff.accs), ncol(pval.diff.accs), 
                          dimnames = dimnames(pval.diff.accs))
wpval.diff.accs  <- apply(diff.accs, c(1,3), function(x) wilcox.test(x)$p.value)
wfdr.diff.accs   <- matrix(p.adjust(wpval.diff.accs, "fdr"), 
                          nrow(pval.diff.accs), ncol(pval.diff.accs), 
                          dimnames = dimnames(pval.diff.accs))

res <- lst.sfull[[2]]$res
get.best.accs <- function(res) {
  res$stats[100,]
  ave.acc  <- as.numeric(as.character(res$best$stats$Accuracy))
  cat.accs <- t(sapply(1:100, function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    acc <- category_stats(res$obs, pred, prob)$Accuracy
    as.numeric(as.character(acc))
  }))
  mean.accs <- rowMeans(cat.accs)
  best.accs <- cat.accs[which.max(mean.accs),]
  names(best.accs) <- categories.title
  best.accs <- c(average=ave.acc, best.accs)
  return(best.accs)
}
full.accs    <- laply(lst.sfull, function(x) get.best.accs(x$res), .parallel=T)
reduced.accs <- laply(1:4, function(ri) {
  laply(lst.sreduced, function(x) get.best.accs(x[[ri]]$res), .parallel=T)
}, .progress="text")
dimnames(reduced.accs)[[1]] <- sub("[RL]\ ", "", sroi.names[1:4]) # Remove R and L 

# subtract each of the fulls from the reduced
titles <- c("Average", categories.title)
diff.accs <- aaply(reduced.accs, .(1), function(x) round(full.accs*100 - x*100, 1))
dimnames(diff.accs) <- list(
  roi=names(reduced.accs), subject=1:nrow(full.accs), category=titles
)
mean.diff.accs2 <- apply(diff.accs, c(1,3), mean)
mean.diff.accs2



# let's try to get the maximum per category
res <- lst.sfull[[2]]$res
get.b
retcats <- t(sapply(1:100, function(li) { # nlambdas
  prob <- res$probs[,,li]
  pred <- res$preds[[li]]
  colnames(prob) <- levels(pred)
  category_stats(res$obs, pred, prob)$Accuracy
}))
colnames(retcats) <- categories.title
plot.ts(retcats)
#plot.ts(retcats, plot.type = "single", col=1:4)
# Second we find the best lambda and get its betas for each category
max.retcats <- apply(retcats, 2, which.max)
