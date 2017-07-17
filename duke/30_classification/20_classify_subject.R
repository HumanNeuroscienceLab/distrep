#!/usr/bin/env Rscript --vanilla

#' The last thing that I can do is manually include or exclude groups of the
#' regions that I have. I could even do every combination of the groups!

# This script will run the classification analysis of the voxelwise data 
# looping through each anatomical region for one subject. It uses glmnet.

# System Arguments
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("continuing with second subject as default")
  subject <- NULL
  nthreads <- 20
} else if (length(args) == 2) {
  subject <- args[1]
  nthreads <- as.numeric(args[2])
  cat(subject, "-", nthreads, "\n")
} else {
  stop("error in reading command-line arguments")
}

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


#' ## Load Data
#' 
#' We get the brain activity (beta-series) for each category.
#' Note for now this is for one subject only.
#' 
#+ load-data
cat("Load\n")

# Let's first get the brain activity for one subject
if (is.null(subject)) subject <- subjects[2]

s <- sprintf
# Reads in the neurosynth ROIs
# There will be a total of 7 ROIs
read.rois <- function(subject) {
  # Read in all ROIs
  roifn    <- file.path(roidir, subject, "neurosynth_visual_funcspace.nii.gz")
  rois     <- read.mask(roifn, NULL)
  rois
}
roi    <- read.rois(subject)
mask   <- roi>10 # exclude L FFA
only.roi <- roi[mask]
vox.grouping <- only.roi

# Combine the ROIs that are bilateral into one ROI
old.grouping <- vox.grouping 
tab <- table(roi.df$names)
for (x in names(tab[tab>1])) {
  inds <- which(roi.df$names==x)
  vals <- roi.df$vals[inds]
  vox.grouping[vox.grouping %in% vals] <- vals[1] # set everything to the first one
}
# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
# get the unique numbers
urois <- sort(unique(vox.grouping))
## redo roi.df
roi.df <- roi.df[roi.df$vals %in% urois,]
roi.df$fullnames <- sroi.names


# Read in the beta-series data
dat <- ldply(categories, function(category) {
  datfn  <- file.path(taskdir, subject, "Categories", s("beta_series_spmg1.reml/beta_series_%s.nii.gz", category))
  dat    <- read.big.nifti(datfn)
  dat    <- deepcopy(dat, mask)
  mdat   <- as.matrix(dat)
  rm(dat); gc()
  mdat
}, .progress="text")
dat <- as.matrix(dat) # observations x

#' We get the associated visual category labels for each trial
#+ setup-labels
# Get the categories
ntrials.per.category <- laply(categories, function(category) {
  datfn  <- file.path(taskdir, subject, "Categories", s("beta_series_spmg1.reml/beta_series_%s.nii.gz", category))
  hdr    <- read.nifti.header(datfn)
  hdr$dim[4]
})
classes <- factor(rep(categories.title, ntrials.per.category))
y    <- as.numeric(classes)

# let's compile the runs associated with each trial
# this way we can do a k-fold cross-validation over the runs
timing <- read.csv("/data1/ffg05/notes/timing.csv")
tmp    <- ddply(timing[,-1], .(category), function(x) data.frame(run=x$run))
tmp    <- tmp[tmp$category!="Circles",]
tmp    <- tmp[tmp$category!="Objects",]
tmp$category <- factor(tmp$category)
if(any(tmp$category != classes)) {
  print(which(tmp$category != classes))
  stop("timing file categories != classes")
}
runs   <- tmp$run


#' ## Functions
#' 
#+ functions
get_stats <- function(prob, pred, obs) {
  require(caret)
  suppressMessages(require(Metrics))
  
  if (!is.factor(pred)) stop("pred must be a factor")
  if (!is.factor(obs)) stop("obs must be a factor")
  
  colnames(prob) <- levels(obs)
  data <- data.frame(
    pred=pred, 
    obs=obs, 
    prob
  )
  
  prob_stats <- lapply(levels(data[, "pred"]), function(class) {
    #Grab one-vs-all data for the class
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs  <- ifelse(data[,  "obs"] == class, 1, 0)
    prob <- data[,class]
    
    #Calculate one-vs-all AUC and logLoss and return
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('ROC', 'logLoss')
    
    return(prob_stats) 
  })
  prob_stats
  
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  #Aggregate and average class-wise stats
  #Todo: add weights
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats)
  #Aggregate overall stats
  overall_stats <- c(CM$overall)
  #Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevalence', 'Detection Prevalence')]
  #Clean names and return
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  
  stats
}

# some redundancy with above
category_stats <- function(obs, pred, prob) {
  suppressMessages(require(Metrics))
  
  ldply(levels(pred), function(class) {
    #Grab one-vs-all data for the class
    pred <- ifelse(pred == class, 1, 0)
    obs  <- ifelse(obs == class, 1, 0)
    prob0 <- prob[,class]
    
    #Calculate accuracy, one-vs-all AUC, and logLoss
    accuracy <- function(obs, pred) sum(pred & obs)/sum(obs)
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(accuracy(obs, pred), auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('Accuracy', 'ROC', 'logLoss')
    
    cbind(class=class, t(prob_stats))
  })
}

# Model Fitting Functions
wrap.model <- function(x, obs, foldsI, alpha=1, nlambda=100) {
  # 1. Get a range of Lambdas
  prefit <- glmnet(x, as.numeric(obs), family="multinomial", alpha=alpha, 
                   nlambda=nlambda)
  lambdas <- prefit$lambda
  
  # 2. Fit the model using a leave-one-run out procedure
  fits <- llply(1:length(foldsI), function(k) {
    indsT   <- foldsI[[k]]
    classesT<- obs[indsT]
    xT      <- x[indsT,]
    fit     <- glmnet(xT, as.numeric(classesT), family="multinomial", 
                      alpha=alpha, lambda=lambdas)
    fit
  }, .parallel=T)
  
  # 3. Get the probabilities and predictions for the outcome (categories)
  probs <- array(NA, c(nrow(x), length(levels(obs)), length(lambdas)))
  preds <- matrix(NA, nrow(x), length(lambdas))
  for (k in 1:length(foldsI)) {
    indsT   <- foldsI[[k]]
    # note: if the fit didn't converge than some lambda fits weren't directly
    # calculated. and so glmnet will use linear interpolations to make the 
    # predictions
    probs0 <- predict(fits[[k]], newx=x[!indsT,], type="response", s=lambdas)
    preds0 <- predict(fits[[k]], newx=x[!indsT,], type="class", s=lambdas)
    #       # to check for which lambda's converged
    #       lambdas0 <- fits[[k]]$lambda
    #       linds    <- lambdas %in% lambdas0 # lambda inds
    oinds          <- which(!foldsI[[k]]) # obs inds
    probs[oinds,,] <- probs0
    preds[oinds,]  <- preds0
    rm(probs0, preds0)
  }
  # relabel the prediction outputs
  preds <- alply(preds, 2, function(xx) {
    factor(as.numeric(xx), levels=unique(as.numeric(obs)), 
           labels=levels(obs))
  })
  
  # 4. Model statistics like accuracy to help determine best lambda etc
  #    (based on caret)...compiles everything into one data frame
  all.stats <- ldply(1:length(preds), function(i) {
    oinds <- !is.na(preds[[i]]) # remove any observations without predictions
    row   <- get_stats(probs[oinds,,i], preds[[i]][oinds], obs[oinds])
    data.frame(alpha=alpha, lambda=lambdas[i], t(row))
  }, .parallel=T)
  
  # 5. Refit the model using the full data
  #    this way we can get the betas and feature information
  final.fits <- glmnet(x, as.numeric(obs), family="multinomial", alpha=alpha, 
                       lambda=lambdas)
  # Number of non-zero features in any class for each lambda
  betas <- coef(final.fits, s=lambdas)
  betas <- laply(betas, as.matrix) # make into array
  nfeats <- aaply(betas, .(3), function(x) sum(colSums(x[,-1]!=0)!=0))
  # add this to the statistics
  all.stats <- cbind(all.stats, nfeats)
  all.stats$lind <- 1:nrow(all.stats)
  #     # If you wanted to do the same thing as above with each fold fit
  #     tmp <- sapply(1:10, function(k) {
  #       betas <- laply(1:5, function(i) as.matrix(fits[[k]]$beta[[i]]))
  #       feats <- apply(betas, c(2,3), function(xx) sum(xx!=0)>0)
  #       nfeats <- colSums(feats)
  #       nfeats
  #     })
  #     rowMeans(tmp)
  
  # 6. Gather information for the best model fit
  # this gives me the average accuracy, logLoss, ROC, and nfeats
  lind <- which.max(all.stats$Accuracy)
  best.stat <- all.stats[lind,]
  # get the betas
  betas   <- coef(final.fits, s=lambdas[lind])
  betas   <- laply(betas, as.matrix) # make into array
  obetas  <- betas      # keep intercept
  betas   <- betas[,-1] # remove the intercept
  
  # 7. Gather information for the best model fit per category
  # get the accuracy, logloss, and ROC
  pred <- preds[[lind]]
  prob <- probs[,,lind]
  colnames(prob) <- levels(pred)
  cat_stats <- category_stats(obs, pred, prob)
  # get the # of non-zero features for each category
  nfeats <- aaply(betas, .(1), function(x) sum(x!=0))
  cat_stats <- cbind(cat_stats, nfeats)
  best.category.stat <- data.frame(alpha=alpha, lambda=lambdas[lind], cat_stats)
  ## We also want to save the overlap between categories
  feats   <- t((betas!=0)*1)
  overlap <- crossprod(feats)
  dimnames(overlap) <- list(levels(obs), levels(obs))
  perc.overlap <- overlap/nrow(feats)
  
  # Return
  list(alpha=alpha, lambdas=lambdas, fold.fits=fits, final.fits=final.fits, 
       stats=all.stats, obs=obs, preds=preds, probs=probs, 
       best=list(stats=best.stat, category.stats=best.category.stat, 
                 betas=betas, intercept=obetas[,1], 
                 overlap=overlap, perc.overlap=perc.overlap))
}

wrap.alphas <- function(x, obs, alphas=c(0,0.5,1), nlambda=100) {
  cat("Fit Models\n")
  res <- llply(alphas, function(alpha) {
    cat("alpha", alpha, "\n")
    wrap.model(x, obs, alpha, nlambda)
  })
  names(res) <- sprintf("alpha%s", as.character(alphas))
  res
}

#' I wonder if it would be best to determine the best fit model, and then figure
#' out what using recursive feature elimination, the particular number of regions
#' that are needed. Let's assume here that we pick one alpha...
#' 
#' For the recursive feature elimination, we first want to run the full model
#' and save those results, then we would want to run 
#' 
#+ rcfe
wrap.rcfe <- function(rois, sdat, obs, foldsI, alpha=1, nlambda=100) {
  # Do everything
  cat("Full Model\n")
  x       <- sdat
  res     <- wrap.model(x, obs, foldsI, alpha, nlambda)
  # save
  full.res  <- res
  full.acc  <- res$best$stats$Accuracy
  full.pval <- res$best$stats$AccuracyPValue
  
  # Do a leave one item out classification
  ref.acc  <- res$best$stats$Accuracy
  ref.urois<- roi.df$vals
  partial.res <- partial.acc <- partial.pval <- partial.rois <- NULL
  
  cat("Partial Model\n")
  for (iter in 1:nrow(roi.df)) {
    cat("+++\nIteration:", iter, "\n")
    cat("fitting\n")
    roi.combns <- llply(1:length(ref.urois), function(i) ref.urois[-i])
    lores <- llply(roi.combns, function(urois) {
      x <- sdat[,rois %in% urois]
      res <- wrap.model(x, obs, foldsI, alpha, nlambda)
      res
    }, .progress="text")
    ## remove those combinations where removal of ROI leads to no change or 
    ## improvement in the accuracy
    cat("determining best\n")
    accs <- laply(lores, function(xx) xx$best$stats$Accuracy)
    keep.rois <- (accs - ref.acc)<0
    if (sum(keep.rois)==0) keep.rois[-which.min(accs)] <- T
    if (sum(!keep.rois)==0) {
      good.rois <- as.character(roi.df$fullnames[roi.df$vals%in%ref.urois[keep.rois]])
      good.rois <- sub(" ", "-", good.rois)
      cat("Keeping ROIs:", paste(good.rois, collapse=", "), "\n")
      
      cat("fitting final\n")
      ref.urois <- ref.urois[keep.rois]
      x <- sdat[,rois %in% ref.urois]
      res <- wrap.model(x, obs, foldsI, alpha, nlambda)
      partial.res <- res
      partial.acc <- res$best$stats$Accuracy
      partial.pval<- res$best$stats$AccuracyPValue
      partial.rois<- ref.urois
      
      cat("nothing to remove, stopping!!!\n")
      break
    } else {
      bad.rois <- as.character(roi.df$fullnames[roi.df$vals%in%ref.urois[!keep.rois]])
      bad.rois <- sub(" ", "-", bad.rois)
      cat("Removing ROIs:", paste(bad.rois, collapse=", "), "\n")
      
      cat("fitting final\n")
      ref.urois <- ref.urois[keep.rois]
      x <- sdat[,rois %in% ref.urois]
      res <- wrap.model(x, obs, foldsI, alpha, nlambda)
      partial.res <- res
      partial.acc <- res$best$stats$Accuracy
      partial.pval<- res$best$stats$AccuracyPValue
      partial.rois<- ref.urois
      
      if (res$best$stats$Accuracy < ref.acc) cat("warning, accuracy lower\n")
      ref.acc <- res$best$stats$Accuracy
      
      if (length(ref.urois) == 1) {
        good.rois <- as.character(roi.df$fullnames[roi.df$vals%in%ref.urois[keep.rois]])
        good.rois <- sub(" ", "-", good.rois)
        cat("Keeping ROIs:", paste(good.rois, collapse=", "), "\n")
        
        cat("left with one roi, stopping!!!\n")
        break
      }
    }
  }
  
  cat("Individual ROI Models\n")
  imods <- llply(roi.df$vals, function(uroi) {
    x <- sdat[,rois==uroi,drop=F]
    res <- wrap.model(x, obs, foldsI, alpha, nlambda)
    list(res=res, 
         acc=res$best$stats$Accuracy, pval=res$best$stats$AccuracyPValue)
  }, .progress="text")
  
  list(
    full=list(res=full.res, acc=full.acc, pval=full.pval), 
    partial=list(res=partial.res, acc=partial.acc, pval=partial.pval, rois=partial.rois), 
    individual=list(res=imods, acc=laply(imods, function(x) x$acc), pval=laply(imods, function(x) x$pval))
  )
}



cat("Set Training Parameters\n")
nrepeats <- 1
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
folds    <- lapply(runFolds, function(x) {
  which(runs %in% x)
})
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})


cat("Scale Data\n")
sdat <- scale(dat)


cat("Run\n")
# Loop through each alpha
alphas  <- seq(0,1,by=0.5)
nlambda <- 100
obs     <- factor(classes)
rois <- vox.grouping

suppressMessages(library(doMC))
registerDoMC(nthreads)

res.rfces <- llply(alphas, function(alpha) {
  cat("\nAlpha", as.character(alpha), "\n")
  wrap.rcfe(rois, sdat, obs, foldsI, alpha=alpha, nlambda=nlambda)
})
names(res.rfces) <- sprintf("alpha%s", as.character(alphas))


cat("Save\n")
subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
save(res.rfces, classes, roi.names, vox.grouping, roi.df, 
     file=file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))

