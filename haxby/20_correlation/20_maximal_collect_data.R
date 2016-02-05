# This script will compile all the haxby data together for easier reading
# It will make use of the neurosynth rois

###
# SETUP
###

# Load the packages etc
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(12)
## paths
base       <- "/data1/ffg05"
base2      <- "/mnt/nfs/psych4/haxby01"
roidir     <- file.path(base2, "rois/mni2func")
scriptdir  <- file.path(base, "scripts/distrep/haxby")
taskdir    <- file.path(base2, "analysis/task_activity")
timedir    <- file.path(base, "notes")
corrdir    <- file.path(base, "analysis/haxby/correlation") # output directory
## subjects
subjects   <- as.character(read.table(file.path(scriptdir, "sublist.txt"))[,1])
## categories
categories      <- c("face", "house", "cat", "bottle", "scissors", "shoe", "chair", "scrambled")
categories.title<- c("Face", "House", "Cat", "Bottle", "Scissors", "Shoe", "Chair", "Scrambled")
## runs
runs       <- c("even", "odd")
runs.title <- c("Even", "Odd")
## # of voxels
lst.nvoxs  <- c(10, 20, 40, 80, 160)


###
# FUNCTIONS
###

s <- sprintf

# Reads in the ventral temporal region
read.roi <- function(subject) {
  roifn    <- file.path(roidir, subject, "ventral_temporal.nii.gz")
  rois     <- read.mask(roifn)
  rois
}



###
# CHECK ROIs
###

# Let's first check the sizes of all the ROIs
roidats <- laply(subjects, function(subject) {
  cat(subject, "\n")
  roi    <- read.roi(subject)
  table(roi*1)
})
print(roidats) # around 1k voxels per subject


###
# DATA SETUP
###

get_maxinds <- function(subject, run, category, vtc, topk=30, use.zstats=TRUE) {
  statdir  <- file.path(taskdir, subject, s("%s_runs.reml", run), "stats")
  beta_fn  <- file.path(statdir, s("coef_%s_gt_all.nii.gz", category))
  
  betas  <- read.nifti.image(beta_fn)
  betas  <- betas[vtc]
  
  if (use.zstats) {
    zstat_fn <- file.path(statdir, s("zstat_%s_gt_all.nii.gz", category))
    
    zstats <- read.nifti.image(zstat_fn)
    zstats <- zstats[vtc]
    
    top_inds <- order(zstats, decreasing=T)[1:topk]
  } else {
    top_inds <- order(betas, decreasing=T)[1:topk]
  }
  
  return(top_inds)
}

get_betas <- function(subject, run, category, mask) {

  
  return(betas)
}

# Betas for every subject
sub.betas <- llply(subjects, function(subject) {
  vtc <- read.roi(subject) # Read in the VTC
  dat <- laply(runs, function(run) {
    laply(categories, function(category) {
      statdir  <- file.path(taskdir, subject, s("%s_runs.reml", run), "stats")
      beta_fn  <- file.path(statdir, s("coef_%s.nii.gz", category))
      betas    <- read.nifti.image(beta_fn)
      betas    <- betas[vtc]
      betas
    })
  })
  dimnames(dat) <- list(run=runs.title, category=categories.title, voxel=NULL)
  dat
}, .progress="text")
names(sub.betas) <- subjects

# Zstats for every subject
sub.zstats <- llply(subjects, function(subject) {
  vtc <- read.roi(subject) # Read in the VTC
  dat <- laply(runs, function(run) {
    laply(categories, function(category) {
      statdir  <- file.path(taskdir, subject, s("%s_runs.reml", run), "stats")
      zstat_fn <- file.path(statdir, s("zstat_%s_gt_all.nii.gz", category))
      zstats   <- read.nifti.image(zstat_fn)
      zstats   <- zstats[vtc]
      zstats
    })
  })
  dimnames(dat) <- list(run=runs.title, category=categories.title, voxel=NULL)
  dat
}, .progress="text")
names(sub.zstats) <- subjects

# Loop through and get the maximal regions of activity
topvox.betas <- llply(lst.nvoxs, function(topk) { # top category-selective voxels to pick
  dat <- laply(subjects, function(subject) {
    betas  <- sub.betas[[subject]]
    zstats <- sub.zstats[[subject]]
    laply(1:length(runs), function(ri) {
      laply(1:length(categories), function(ci) {
        # category-selective activity pattern from ri run
        activity_map <- zstats[ri,ci,]
        # find top voxels (maximally selective)
        top_inds <- order(activity_map, decreasing=T)[1:topk]
        # select those voxels in the beta map of both even and odd runs
        max.betas.even <- betas[runs=="even",ci,top_inds]
        max.betas.odd  <- betas[runs=="odd",ci,top_inds]
        
        cbind(even=max.betas.even, odd=max.betas.odd)
      })
    })
  }, .parallel=TRUE)
  dimnames(dat) <- list(subject=subjects, run.mask=runs.title, 
                        category=categories.title, top.voxel=NULL, 
                        run.data=runs.title)
  dat <- aperm(dat, c(1,2,3,5,4))
  dat
}, .progress="text")
names(topvox.betas) <- lst.nvoxs


# Now we concatenate the data across subjects for later group analysis
library(abind)
dim(topvox.betas[[1]]) # subjs, run.mask, categ, vox, run.data
# scale each subject's / roi's data
grp.topvox <- aaply(topvox.betas[[1]], .(2,3,4), function(x) {
  sx <- apply(x, 2, scale)
  as.vector(sx)
})
names(dimnames(grp.topvox))[[4]] <- "top.voxel"


###
# SAVE
###

save(subjects, roi.names, roi.names.title, categories, categories.title, lst.nvoxs, 
     sub.betas, sub.zstats, topvox.betas, grp.topvox, 
     file=file.path(corrdir, "maximal_activity_patterns.rda"))


