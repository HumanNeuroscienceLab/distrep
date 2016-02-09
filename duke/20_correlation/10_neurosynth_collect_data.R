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
roidir     <- file.path(base, "analysis/rois")
scriptdir  <- file.path(base, "scripts/distrep/duke")
taskdir    <- file.path(base, "analysis/task_activity")
timedir    <- file.path(base, "notes")
corrdir    <- file.path(base, "analysis/duke/correlation") # output directory
## subjects
subjects   <- as.character(read.table(file.path(scriptdir, "sublist.txt"))[,1])
## rois
roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", "R")
roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", "M1", "M1", "A1", "A1")
roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
                            fullnames=paste(roi.hemis, roi.names))
roi.df        <- roi.df[-1,]
roi.names.title <- as.character(roi.df$fullnames)
roi.names     <- sub(" ", "_", tolower(roi.names.title))
## categories
categories <- c("faces", "fruits", "letters", "vehicles")
categories.title <- c("Faces", "Fruits", "Letters", "Vehicles")
## runs
runs       <- c("even", "odd")
runs.title <- c("Even", "Odd")


###
# FUNCTIONS
###

s <- sprintf

# Reads in the neurosynth ROIs
read.rois <- function(subject) {
  roifn    <- file.path(roidir, subject, "neurosynth_visual_funcspace.nii.gz")
  rois     <- read.mask(roifn, NULL)
  rois
}


###
# CHECK ROIs
###

# Let's first check the sizes of all the ROIs
roidats <- laply(subjects, function(subject) {
  cat(subject, "\n")
  roi    <- read.rois(subject)
  table(roi)
})
roidats <- roidats[,-c(1:2)] # 1: 0s and 2: L FFA
colnames(roidats) <- roi.df$fullnames
print(roidats)
# note: the L FFA voxel size can be down to 5 voxels...so it's been taken out.


###
# DATA SETUP
###

# We get the brain activity for each roi, run, subject, and category
# and collect it together
lldat <- llply(subjects, function(subject) {
  cat(subject, "\n")
  
  # Read in all the ROIs
  roi    <- read.rois(subject)
  
  ldat <- llply(1:length(roi.names), function(ri) {
    # Get ROI
    # roi.name <- roi.names[ri]
    one.roi <- roi == roi.vals[ri]
    
    dat <- laply(runs, function(run) {
      laply(categories, function(category) {
        datfn  <- file.path(taskdir, subject, "Categories", 
                            s("%sruns_task_activity_spmg1.reml", run), 
                            "stats", s("coef_%s.nii.gz", category))
        dat    <- read.nifti.image(datfn)
        dat    <- dat[one.roi]
        dat
      })
    }, .parallel=T)
    dimnames(dat) <- list(run=runs.title, category=categories.title, voxel=NULL)
    dat <- aperm(dat, c(3,2,1))
    dat
  }, .progress="text")
  names(ldat) <- roi.names
  ldat
})
names(lldat) <- subjects


# Now we concatenate the data across subjects
# for later group level analyses
library(abind)
dim(lldat[[1]][[1]]) # note: voxels x categories x even/odd
# scale each subject's / roi's data
cldat <- llply(1:length(roi.names), function(ri) {
  ldat <- llply(1:length(subjects), function(si) {
    dat <- lldat[[si]][[ri]]
    dat <- aaply(dat, 3, scale)
    dat <- aperm(dat, c(2,3,1))
    dat
  }, .parallel = T)
  do.call(abind, c(ldat, along=1))
})


###
# SAVE
###

save(subjects, roi.names, roi.names.title, roi.df, categories, categories.title, 
     lldat, cldat, 
     file=file.path(corrdir, "neurosynth_activity_patterns.rda"))

