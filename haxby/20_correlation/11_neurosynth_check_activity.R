# Let's check if the activity levels are looking okay
# That is for FFA is faces the highest activity level

# Load packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)

# load the ridge.pcor function based on the parcor package
ridge.pcor <- NULL
source("/data1/ffg05/scripts/repsim/R/pcor.R")

# Set some variables
base       <- "/data1/ffg05"
corrdir    <- file.path(base, "analysis/haxby/correlation") # output directory

# Now load the data from the previous script
# we set all the variables to be loaded to NULL for Rstudio code checker
lldat <- NULL; cldat <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL; subjects <- NULL
load(file.path(corrdir, "neurosynth_activity_patterns.rda"))


# We get the activity levels averaged across subjects and runs
# Note that these are from the betas, which aren't face > everything-else
activity.levels <- laply(roi.names, function(ri) {
  rowMeans(sapply(lldat, function(x) apply(x[[ri]], c(2), mean)))
})
rownames(activity.levels) <- roi.names.title
#library(pander)
#pandoc.table(activity.levels, round=1)
#library(htmlTable)
#htmlTable(round(activity.levels, 1))

# Collapse L/R rois together
# Get the average of the activity levels across some of the ROIs
new.rois <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
tmp <- activity.levels
activity.levels <- matrix(0, length(new.rois), ncol(activity.levels))
activity.levels[1,] <- tmp[1,]
activity.levels[2,] <- colMeans(tmp[2:3,])
activity.levels[3,] <- colMeans(tmp[4:5,])
activity.levels[4,] <- tmp[6,]
activity.levels[5,] <- colMeans(tmp[7:8,])
activity.levels[6,] <- colMeans(tmp[9:10,])
activity.levels[7,] <- colMeans(tmp[11:12,])
dimnames(activity.levels) <- list(roi=new.rois, category=colnames(tmp))

# Plot Table
library(pander)
pandoc.table(activity.levels, round=1)
library(htmlTable)
htmlTable(round(activity.levels, 1))

# Convert to data-frame
adf <- reshape2::melt(activity.levels, value.name="activity")
adf$category <- factor(as.character(adf$category), 
                       levels=rev(levels(adf$category)))

# Plot
library(ggplot2)
library(RColorBrewer)
mytheme <- theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title = element_blank())

p <- ggplot(adf, aes(roi, category, fill = activity)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradientn(name="Beta", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92") + 
  ggtitle("Brain Activity") +
  xlab("Brain Regions") +
  ylab("") + 
  mytheme
plot(p)
dir.create("/data1/ffg05/analysis/haxby/correlation/plots", showWarnings=F)
ggsave(filename=file.path(corrdir, "plots", "activity_levels.pdf"), 
       plot=p, width=9, height=6)
