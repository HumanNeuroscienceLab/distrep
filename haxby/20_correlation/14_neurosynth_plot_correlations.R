# Our goal for this script is to generate plots comparing the patterns of 
# activity between pairs of categories for every neurosynth based ROI.

# We generate outputs for the 'group-average' thresholded for significant 
# results as well as the percent correct category detection measure.
#
# Note that this builds on `12_neurosynth_partial_correlations.R`.


###
# SETUP
###

# using dev version of ggplot??
# install_github("hadley/scales")
# install_github("hadley/ggplot2")

# If no libraries found...load below:
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
library(ggplot2)
library(RColorBrewer)
require(grid)
require(gridExtra)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)

# Set some variables
base       <- "/data1/ffg05"
corrdir    <- file.path(base, "analysis/haxby/correlation") # output directory

# Load the data
ridge.grp <- NULL; ridge.subs <- NULL; ridge.detect <- NULL
load(file.path(corrdir, "neurosynth_partial_correlations.rda"))
roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", "R")
roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", "M1", "M1", "A1", "A1")
roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
                            fullnames=paste(roi.hemis, roi.names))
roi.df        <- roi.df[-1,] # remove L FFA since it's small


###
# Functions
###

# generic plot theme
mytheme <- theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title = element_blank())


###
# WITHIN-CATEGORY CORRELATIONS
###

# We first focus on plotting information found in a given category (same 
# category).
# 
# I start here by getting the relevant data, the diagonal of each region's 
# partial correlation matrix.
zthr <- 1.96
same.cat.grp.thr <- aaply(ridge.grp, .(1), function(x) {
  diag(x[,,1] * (x[,,2]>zthr))
})

# We convert the data array to a data frame. And apply the following:
# 
# - Any 0 correlations are set to NA so they are colored grey when plotting.
# - Create column with the ROI hemisphere
# - Create column with the ROI lobe
# - Create column with ROI classification (e.g., visual or motor)
# - Create column with index of each ROI
# - Create column to use to order the ROIs based on the hemi, lobes, and 
# indices.
# - Finally, we reorder the categories for a nicer display.
ncats    <- dim(same.cat.grp.thr)[2]
dat      <- reshape2::melt(same.cat.grp.thr)
colnames(dat)[2:3] <- c("category", "correlation")
dat$correlation[dat$correlation==0] <- NA  ## threshold
dat$category <- factor(as.character(dat$category), ## reorder for better display
                       levels=rev(categories.title))
# if want to see ordering
head(dat, 2)

# Now we can plot all these new results. Again this is a plot of within 
# category information across our regions.
lims <- range(ridge.grp[,,,1]*(ridge.grp[,,,2]>zthr), na.rm=T)
if (lims[1]<0) stop("min is less than 0")
lims[1] <- 0
p <- ggplot(dat, aes(roi, category, fill = correlation)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradientn(name="Correlation", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92", 
                       limits=lims) + 
  ggtitle("Comparisons Between the Same Category") +
  xlab("Brain Regions") +
  ylab("") + 
  mytheme #+ 
#theme(axis.text.x = element_blank())
plot(p)
ggsave(filename=file.path(corrdir, "neurosynth_plot_correlations_within_category.pdf"), 
       plot=p, width=14.9, height=6.69)



###
# PERCENT ACCURACY FOR CATEGORY DETECTION
###

# Format the data
ave.ridge.detect <- aaply(ridge.detect, .(2,3), mean)*100
ncats    <- dim(ave.ridge.detect)[2]
dat      <- reshape2::melt(ave.ridge.detect)
colnames(dat)[2:3] <- c("category", "accuracy")
dat$accuracy[dat$accuracy==0] <- NA  ## threshold
dat$category <- factor(as.character(dat$category), ## reorder for better display
                       levels=rev(categories.title))

# Plot
lims <- range(ave.ridge.detect, na.rm=T)
if (lims[1]<0) stop("min is less than 0")
lims[1] <- 0
p <- ggplot(dat, aes(roi, category, fill = accuracy)) +
  geom_tile(colour="white", size=0.5) + 
  scale_fill_gradientn(name="Percent Accuracy", 
                       colours=rev(brewer.pal(11, "RdYlBu")), 
                       na.value="grey92", 
                       limits=lims) + 
  ggtitle("Percent Correct Category Detection") +
  xlab("Brain Regions") +
  ylab("") + 
  mytheme
plot(p)
ggsave(filename=file.path(corrdir, "neurosynth_plot_category_detection.pdf"), 
       plot=p, width=14.9, height=6.69)

