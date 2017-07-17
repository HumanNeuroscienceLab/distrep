#!/usr/bin/env R

#subject <- "44172"
subjects <- as.character(read.table("/data1/ffg05/scripts/distrep/duke/sublist.txt")[,1])

# ROI vals
roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", "R")
roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", "M1", "M1", "A1", "A1")
roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
                            fullnames=paste(roi.hemis, roi.names))
roi.df        <- roi.df[-1,] # skip L FFA

# Basic paths
basedir <- "/data1/ffg05/analysis"
roidir  <- file.path(basedir, "rois")
taskdir <- file.path(basedir, "task_activity")
predir  <- file.path(basedir, "preprocessed")

for (subject in subjects) {
  cat(subject, "\n")
  
  # Get the input files
  roifile <- file.path(roidir, subject, "neurosynth_visual_funcspace.nii.gz")
  maskfile<- file.path(predir, subject, "Categories/mask.nii.gz")
  funcfile<- file.path(predir, subject, "Categories/filtered_func_data.nii.gz")

  # Get the output files
  outdir <- file.path(taskdir, subject, "Categories/task_activity.rois")
  if (!file.exists(outdir)) dir.create(outdir)

  # Read in the stuff
  suppressMessages(library(niftir))
  mask <- read.mask(maskfile) # not really needed
  rois <- as.vector(read.nifti.image(roifile))
  func <- read.big.nifti(funcfile)

  # Loop through each unique name
  library(plyr)
  ts.rois <- t(daply(roi.df, .(names), function(x) {
    rowMeans(func[,rois %in% x$vals])
  }))
  
  # Remove the column mean
  #ts.rois <- scale(ts.rois, scale=F, center=T)

  # Reorder
  ts.rois <- ts.rois[,c("FFA","PPA","LOC","VWF","V1","M1","A1")]

  # Save
  outfile <- file.path(outdir, "neurosynth_rois.1D")
  write.table(ts.rois, file=outfile, row.names=F, col.names=F)
}
