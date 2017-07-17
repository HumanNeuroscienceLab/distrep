#!/usr/bin/env Rscript

# Let's read in the univariate and multivariate results

### Load Packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
library(RColorBrewer)
library(ggplot2)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)


#' ### Variables
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
categories <- c("faces", "fruits", "letters", "vehicles")
categories.title <- c("Faces", "Fruits", "Letters", "Vehicles")
# cols
dat.cols <- rev(brewer.pal(11, "RdYlBu"))
cat.cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # Faces, Letters, Fruits, & Vehicles
#roi.cols <- brewer.pal(8, "Set3")[-c(1,7)]
roi.cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)]
# subjects
subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])
## output for saving
outdir <- "/data1/ffg05/analysis/combined/Rout"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/combined/plots"
dir.create(plotdir, showWarnings = FALSE)

### Load Classification Data
grouping <- NULL
infile <- file.path("/data1/ffg05/analysis/duke/classification", "group_fullmodels.rda")
load(infile)
urois <- sort(unique(grouping)) # is this needed?


### Load GLM Data

## ROIs
roipath <- "/data1/ffg05/analysis/rois/group/neurosynth_visual_rois_3mm.nii.gz"
old.rois<- read.mask(roipath, NULL)
mask    <- old.rois!=0
## replace old labels with newer scheme
new.vals <- c(0, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7)
rois    <- vector("numeric", length(old.rois))
for (i in 1:length(roi.vals)) {
  rois[old.rois==roi.vals[i]] <- new.vals[i]
}
sroi.names <- c("FFA", "PPA", "LOC", "VWFA", "V1", "M1/2", "A1/2")
## get other grouping
grouping2 <- grouping*0
for (i in 1:length(urois)) {
  grouping2[grouping==urois[i]] <- i
}

## Task Activity
base <- "/data1/ffg05/analysis/combined/task_std_3mm"
subjects <- as.character(read.table("/data1/ffg05/scripts/distrep/duke/sublist.txt")[,1])

roi.dats <- laply(1:length(subjects), function(si) {
  subfnames <- sprintf("%s/%s_zstat_%s_gt_all.nii.gz", base, subjects[si], categories)
  
  dats <- lapply(subfnames, function(sf) read.nifti.image(sf)[mask])
  
  roi.dats <- laply(1:length(dats), function(di) {
    unlist(llply(1:length(sroi.names), function(ri) {
      dats[[di]][rois[mask]==ri]
    }))
  })
  
  roi.dats
}, .parallel=T)
dimnames(roi.dats) <- list(subject=subjects, category=categories.title, vox=NULL)

grp.roi.dats <- apply(roi.dats, c(2,3), mean)

dim(grp.roi.dats)
tmp1 <- grp.res$alpha1$full$res$best$betas[,grouping2==1]
tmp2 <- grp.res$alpha1$individual[[1]]$res$best$betas


lapply(1:4, function(i) {
  table(abs(grp.roi.dats[i,])>1, grp.res$alpha1$full$res$best$betas[i,]!=0)
})


# TODO: for each subject, want to see the overlap between glm with individual
# and with the group
# RESULT: expect that more overlap of glm with combined than individual.
