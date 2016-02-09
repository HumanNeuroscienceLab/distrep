# This script will compile all the duke data together for easier reading
# It will make use of the maximal responses as the rois

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
roidir     <- file.path(base, "analysis/rois")
scriptdir  <- file.path(base, "scripts/distrep/duke")
taskdir    <- file.path(base, "analysis/task_activity")
timedir    <- file.path(base, "notes")
corrdir    <- file.path(base, "analysis/duke/correlation") # output directory
## subjects
subjects   <- as.character(read.table(file.path(scriptdir, "sublist.txt"))[,1])
## categories
categories       <- c("faces", "fruits", "letters", "vehicles")
categories.title <- c("Faces", "Fruits", "Letters", "Vehicles")
## runs
runs       <- c("even", "odd")
runs.title <- c("Even", "Odd")
## # of voxels
lst.nvoxs  <- c(10, 20, 30, 40, 80, 160, 320, 640)


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

# Betas for every subject
sub.betas <- llply(subjects, function(subject) {
  vtc <- read.roi(subject) # Read in the VTC
  dat <- laply(runs, function(run) {
    laply(categories, function(category) {
      statdir  <- file.path(taskdir, subject, "Categories", 
                            s("%sruns_task_activity_spmg1.reml", run), 
                            "stats")
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
      statdir  <- file.path(taskdir, subject, "Categories", 
                            s("%sruns_task_activity_spmg1.reml", run), 
                            "stats")
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
library(abind)
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
        max.betas.even <- betas[runs=="even",,top_inds]
        max.betas.odd  <- betas[runs=="odd",,top_inds]
        
        abind(even=max.betas.even, odd=max.betas.odd, along=-1)
      })
    })
  }, .parallel=TRUE)
  dimnames(dat) <- list(subject=subjects, 
                        run.mask=runs.title, category.mask=categories.title, 
                        run.data=runs.title, category.data=categories.title, 
                        top.voxel=NULL)
  dat
}, .progress="text")
names(topvox.betas) <- lst.nvoxs


# Now we concatenate the data across subjects for later group analysis
dim(topvox.betas[[1]]) # subjs, run.mask, categ, vox, run.data
# scale each subject's / roi's data
grp.topvox.betas <- llply(topvox.betas, function(xx) {
  dat <- aaply(xx, .(2,3,4,5), function(x) {
    sx <- apply(x, 2, scale)
    as.vector(sx)
  }, .parallel=TRUE)
  names(dimnames(dat))[5] <- "voxel"
  dat
})



###
# SAVE
###

save(subjects, runs, runs.title, categories, categories.title, lst.nvoxs, 
     sub.betas, sub.zstats, topvox.betas, grp.topvox.betas, 
     file=file.path(corrdir, "maximal_activity_patterns.rda"))
