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
corrdir    <- file.path(base, "analysis/duke/correlation") # output directory

# Now load the data from the previous script
# we set all the variables to be loaded to NULL for Rstudio code checker
lldat <- NULL; cldat <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL; subjects <- NULL
load(file.path(corrdir, "neurosynth_activity_patterns.rda"))


# We get the activity levels averaged across subjects and runs
# Note that these are from the betas, which aren't face > everything-else
activity.levels <- laply(rois, function(ri) {
  rowMeans(sapply(lldat, function(x) apply(x[[ri]], c(2), mean)))
})
rownames(activity.levels) <- roi.names.title
library(pander)
pandoc.table(activity.levels, round=1)
library(htmlTable)
htmlTable(round(activity.levels, 1))
