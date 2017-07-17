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
library(ggthemes)
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
titles <- c("Average", categories.title)
# new roi names
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
# colors
roi.cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)]
cat.cols <- brewer.pal(8, "Set2")[c(8,3,2,5,1)] # Total, Faces, Letters,
# outputs
outdir <- "/data1/ffg05/analysis/duke/classification"
plotdir <- file.path(outdir, "plots")

cat("Source\n")
summarySEwithin <- NULL; bar_theme <- NULL; mytheme <- NULL
source("/data1/ffg05/scripts/distrep/duke/30_classification/z_barplot_funs.R")



#' ## Load Data
#' 
#+ load
full.accs <- NULL; reduced.accs <- NULL; diff.accs <- NULL
infile <- file.path(outdir, "category_accuracies_leave_region_out.rda")
load(infile)
infile <- file.path(outdir, "category_accuracies_individual_region.csv")
indiv.df <- read.csv(infile)

#' ## Stats
#' 
#+ stats
mean.diff.accs <- apply(diff.accs, c(1,3), mean)
sd.diff.accs <- apply(diff.accs, c(1,3), sd)

round(mean.diff.accs, 2)
round(sd.diff.accs, 2)

tval.diff.accs  <- apply(diff.accs, c(1,3), function(x) t.test(x, alternative="greater")$statistic)
pval.diff.accs  <- apply(diff.accs, c(1,3), function(x) t.test(x, alternative="greater")$p.value)
fdr.diff.accs   <- matrix(p.adjust(pval.diff.accs, "fdr"), 
                          nrow(pval.diff.accs), ncol(pval.diff.accs), 
                          dimnames = dimnames(pval.diff.accs))
wpval.diff.accs  <- apply(diff.accs, c(1,3), function(x) wilcox.test(x, alternative="greater")$p.value)
wfdr.diff.accs   <- matrix(p.adjust(wpval.diff.accs, "fdr"), 
                           nrow(pval.diff.accs), ncol(pval.diff.accs), 
                           dimnames = dimnames(pval.diff.accs))


#' ## Plots
#' 
#' ### Bar Plot of Average Data
#' 
#+ bar-plot
dfw_long <- reshape2::melt(t(diff.accs[,,1]), value.name="change.accuracy")
names(dfw_long)[1] <- "subject"
names(dfw_long)[2] <- "roi"
dfwc <- summarySEwithin(dfw_long, measurevar="change.accuracy", withinvars="roi",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc

max.ylim <- (max(dfwc$change.accuracy)+max(dfwc$se))*1.02
# plot
cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
p <- ggplot(dfwc, aes(x=roi, y=change.accuracy, fill=roi)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=change.accuracy-se, ymax=change.accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim = c(-2, max.ylim)) + 
  scale_fill_manual(values=cols) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank(), 
    legend.position="none"
  )
# plot
plot(p)
# save
outfile <- file.path(plotdir, "region_importance_leave_one_out.pdf")
ggsave(outfile, p, width=5, height=3)


#' ### Other Plot for Each Category
#'
#'  Same as the correlation plots.
#'
#+ other-plot
library(scales)
mytheme <- theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title = element_blank())
idf        <- reshape2::melt(mean.diff.accs[,-1], value.name="accuracy")
names(idf)[1:2] <- c("roi", "category")
idf$pvalue <- reshape2::melt(pval.diff.accs[,-1])$value
idf$fdr.pvalue <- p.adjust(idf$pvalue, method="fdr")
idf$label  <- ifelse(idf$fdr.pvalue<0.05, '*', '')
idf$category    <- factor(idf$category, levels=rev(levels(idf$category)))
idf$accuracy <- idf$accuracy * -1

p <- ggplot(idf, aes(roi, category, fill = accuracy)) +
  geom_tile(colour="white", size=1) +
  geom_text(aes(label=label), size=15, colour="white") + 
  scale_fill_gradient2(low = brewer.pal(9,"Blues")[9], high=brewer.pal(9,"Reds")[9]) + 
  #scale_fill_gradientn(name="Accuracy", 
  #                     colours=brewer.pal(9,"YlGn"), 
  #                     na.value="grey92") + #, 
  #                     #breaks=c(0.5,0.4,0.3)) + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") +
  mytheme 
plot(p)
## save
ggsave(filename=file.path(plotdir, "region_importance_leave_one_out_per_category.pdf"), 
       plot=p, width=12, height=6)



#' ## Region Importance
#' 
#+ plot-region-importance
dfw_long <- reshape2::melt(diff.accs, value.name="change.accuracy")
dfwc <- summarySEwithin(dfw_long, measurevar="change.accuracy", 
                        withinvars=c("roi", "category"), idvar="subject", 
                        na.rm=FALSE, conf.interval=.95)
dfwc
# plot
p <- ggplot(dfwc, aes(x=roi, y=change.accuracy, color=category, shape=category)) +
  geom_pointrange(aes(ymin=change.accuracy-se, ymax=change.accuracy+se), 
                  position=position_dodge()) + 
  scale_shape_manual(values=c(3,15,16,17,18)) + 
  scale_color_manual(values=cat.cols) + 
  #coord_cartesian(ylim = c(-2, 18)) + 
  scale_y_continuous(breaks=seq(0,16,by=4)) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank()
  )
plot(p)
# save
ggsave(filename=file.path(plotdir, "region_imp_by_cat.pdf"), 
       plot=p, width=6, height=4)

#' ## Individual Region Category Accuracy
#' 
#+ plot-individual-accuracy
dfwc2 <- summarySEwithin(indiv.df, measurevar="accuracy", 
                          withinvars=c("roi", "category"), idvar="subject", 
                          na.rm=FALSE, conf.interval=.95)
dfwc2
# plot
# plot
p <- ggplot(dfwc2, aes(x=roi, y=accuracy, color=category, shape=category)) +
  geom_pointrange(aes(ymin=accuracy-se, ymax=accuracy+se)) + 
  scale_shape_manual(values=c(3,15,16,17,18)) + 
  scale_color_manual(values=cat.cols) + 
  #coord_cartesian(ylim = c(-2, 18)) + 
  #scale_y_continuous(breaks=seq(0,16,by=4)) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank()
  )
plot(p)
# save
outfile <- file.path(plotdir, "individual_accuracy_by_cat.pdf")
ggsave(outfile, p, width=6, height=4)






