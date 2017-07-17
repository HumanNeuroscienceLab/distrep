
#' ## Setup
#' 
#' ### Packages
#' 
#+ packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)
nthreads <- 24

library(plyr)
library(ggplot2)
library(RColorBrewer)
library(caret)
library(glmnet)
library(Matrix)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(nthreads)

source("/data1/ffg05/scripts/distrep/duke/30_classification/z_barplot_funs.R")


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
outdir <- "/data1/ffg05/analysis/task_activity/group/Rout"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/task_activity/group/plots"
dir.create(plotdir, showWarnings = FALSE)


#' ### ROI
#' 
#+ roi
roipath <- "/data1/ffg05/analysis/rois/group/neurosynth_visual_rois.nii.gz"
old.rois<- read.mask(roipath, NULL)
## replace old labels with newer scheme
new.vals <- c(0, 1, 2, 2, 3, 3, 4, 5, 5, 6, 6, 7, 7)
rois    <- vector("numeric", length(old.rois))
for (i in 1:length(roi.vals)) {
  rois[old.rois==roi.vals[i]] <- new.vals[i]
}
sroi.names <- c("FFA", "PPA", "LOC", "VWFA", "V1", "M1/2", "A1/2")


#' ## Task
#' 
#' ### Load
#' 
#+ task-load
base <- "/data1/ffg05/analysis/task_activity/group/categories_task_smoother.mema/subjects"
subdirs <- list.files(base, pattern="_stats", full.names=T)
fnames  <- sprintf("zstat_%s_gt_all.nii.gz", categories)

roi.dats <- laply(1:length(subjects), function(si) {
  subfnames <- file.path(subdirs[si], fnames)
  subject   <- sub("_stats", "", basename(subdirs[si]))
  maskfile  <- file.path(dirname(subdirs[si]), sprintf("%s_mask.nii.gz", subject))
  
  mask <- read.mask(maskfile)
  dats <- lapply(subfnames, function(sf) read.nifti.image(sf)[mask])
  
  roi.dats <- laply(1:length(dats), function(di) {
    laply(1:length(sroi.names), function(ri) {
      mean(dats[[di]][rois[mask]==ri])
    })
  })
  
  roi.dats
}, .parallel=T)
dimnames(roi.dats) <- list(subject=subjects, category=categories.title, roi=sroi.names)

#' ### Plot
#' 
#+ task-plot
dfw_long <- reshape2::melt(roi.dats, value.name="zstats")
dfwc <- summarySEwithin(dfw_long, measurevar="zstats", withinvars=c("category", "roi"),
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc$category <- factor(dfwc$category, levels=categories.title[c(1,3,2,4)], 
                        labels=categories.title[c(1,3,2,4)])
dfwc
## plot
p <- ggplot(dfwc, aes(x=roi, y=zstats, fill=category)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_linerange(aes(ymin=zstats-se, ymax=zstats+se), 
                position=position_dodge(.9)) + 
#  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylab("Zstats") + 
  bar_theme() + 
  scale_fill_manual(values=cat.cols) + 
  theme(
    axis.title.x = element_blank()
  )
plot(p)
# save
outfile <- file.path(plotdir, "univariate_roi_barplot.pdf")
ggsave(outfile, p, width=8, height=6)

