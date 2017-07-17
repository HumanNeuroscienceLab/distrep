#!/usr/bin/env Rscript

if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

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
plotdir <- "/data1/ffg05/analysis/task_activity/group/plots"


library(plyr)
library(abind)
dat <- laply(subjects, function(subject) {
  # Get the input files
  roifile <- file.path(roidir, subject, "neurosynth_visual_funcspace.nii.gz")
  maskfile<- file.path(predir, subject, "Categories/mask.nii.gz")
  funcfile<- file.path(predir, subject, "Categories/filtered_func_data.nii.gz")

  # Get the input files
  indir   <- file.path(taskdir, subject, "Categories/task_activity.rois")
  bucfile <- file.path(indir, "stats_bucket.1D")
  xmatfile<- file.path(indir, "xmat.1D")

  # system(sprintf("cat %s", bucfile))
  # system(sprintf("head -n 20 %s", xmatfile))

  # Read in the bucket
  ## not sure why there are 8 rows
  ## split the fstat, tstat, and betas
  bucket <- read.table(bucfile)
  fstat  <- bucket[,1]
  betas  <- bucket[,seq(2,ncol(bucket),by=2)]
  tstats <- bucket[,seq(3,ncol(bucket),by=2)]

  # Assign the right names
  cnames <- c("circles", "faces", "fruits", "letters", "objects", "vehicles", "circles_gt_all", "faces_gt_all", "fruits_gt_all", "letters_gt_all", "objects_gt_all", "vehicles_gt_all")
  colnames(betas)  <- cnames
  colnames(tstats) <- cnames
  
  # Return
  abind(betas, tstats, along=3)
}, .progress="text")

dimnames(dat) <- list(
  subject=subjects, 
  roi=c("FFA","PPA","LOC","PPA","VWF","V1","M1","A1"), 
  contrast=dimnames(dat)[[3]], 
  stat=c("beta", "tstat")
)

# Let's select the contrasts and stat that we want
## note: i duplicated the ppa column so remove one of them
sdat <- dat[,-4,c(8,10,9,12),]

# Ok so let's see the significant effects
apply(sdat[,,,1], c(2,3), function(x) t.test(x)$statistic)
qt(0.05/2, 19, lower.tail=F) # signif threshold
apply(sdat[,,,1], c(2,3), function(x) qt(binom.test(sum(x>0), length(subjects), alternative="greater")$p.value, 19, lower.tail=F))

# So we'll end up using the t-statistic approach to get the significance
grp.tstats <- apply(sdat[,,,1], c(2,3), function(x) t.test(x)$statistic)
thr <- qt(0.05/2, 19, lower.tail=F)
grp.sig <- grp.tstats > thr

# And then we'll need to convert each subject's tstat to a zstat and average
# note the df is 1958
sign.tdat <- sign(sdat[,,,2])
pval.tdat <- pt(sdat[,,,2], 1958, lower.tail=F)/2
pval.tdat[sign.tdat<0] <- pt(sdat[,,,2][sign.tdat<0], 1958, lower.tail=T)/2
zval.tdat <- qt(pval.tdat, Inf, lower.tail=F) * sign.tdat

# Average the z-values
mean.zvals <- apply(zval.tdat, 2:3, mean)
sd.zvals   <- apply(zval.tdat, 2:3, sd)
sem.zvals  <- sd.zvals/sqrt(20)


#' ### Plot
#' 
#+ task-plot
library(ggplot2)
library(RColorBrewer)
source("/data1/ffg05/scripts/distrep/duke/30_classification/z_barplot_funs.R")
dimnames(zval.tdat)$contrast <- c("Faces", "Letters", "Fruits", "Vehicles")
dimnames(zval.tdat)$roi[c(4,6,7)] <- c("VWFA", "M1/2", "A1/2")
dfw_long <- reshape2::melt(zval.tdat, value.name="zstats")
dfwc <- summarySEwithin(dfw_long, measurevar="zstats", withinvars=c("roi", "contrast"),
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
#dfwc$contrast <- factor(dfwc$contrast, levels=categories.title[c(1,3,2,4)], 
#                        labels=categories.title[c(1,3,2,4)])
dfwc$ctype <- factor(dfwc$roi, levels=levels(dfwc$roi), labels=rep(c("category-selective", "control"),c(4,3)))
cat.cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # Faces, Letters, Fruits, & Vehicles
## plot and save category-selective
p <- ggplot(subset(dfwc, ctype=="category-selective"), aes(x=roi, y=zstats, fill=contrast)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_linerange(aes(ymin=zstats-se, ymax=zstats+se), 
                 position=position_dodge(.9)) + 
  #  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylim(min(dfwc$zstats-dfwc$se), max(dfwc$zstats+dfwc$se)) + 
  ylab("Zstats") + 
  bar_theme() + 
  scale_fill_manual(values=cat.cols) + 
  theme(
    axis.title.x = element_blank()
  )
plot(p)
outfile <- file.path(plotdir, "univariate_roi_barplot_category-selective.pdf")
ggsave(outfile, p, width=8, height=6)
## plot and save control
p <- ggplot(subset(dfwc, ctype=="control"), aes(x=roi, y=zstats, fill=contrast)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_linerange(aes(ymin=zstats-se, ymax=zstats+se), 
                 position=position_dodge(.9)) + 
  #  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylim(min(dfwc$zstats-dfwc$se), max(dfwc$zstats+dfwc$se)) + 
  ylab("Zstats") + 
  bar_theme() + 
  scale_fill_manual(values=cat.cols) + 
  theme(
    axis.title.x = element_blank()
  )
plot(p)
outfile <- file.path(plotdir, "univariate_roi_barplot_control.pdf")
ggsave(outfile, p, width=8, height=6)
