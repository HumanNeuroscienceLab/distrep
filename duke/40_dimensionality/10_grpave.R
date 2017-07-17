#' This page will deal with getting the output from the group-level analyses
#' 
#' ## Setup
#' 
#' First, we load all the packages
#' 
#+ load-packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)
library(plyr)
library(RColorBrewer)
library(corrplot)
library(ggplot2)
library(caret)
library(glmnet)
library(Matrix)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)

#' Let's load some needed functions including
#' 
#' - `load.data`: Does what it says, loads the data for one subject
#' - `get_stats`: Used to get classification performance
#' - `category_stats`: Used to get classification performance for each category
#' 
#+ load-funs
load.data.std <- NULL
source("/data1/ffg05/scripts/repsim/scratch/stats_lab_funs.R")

#' Now we gather useful variables like the list of categories for the 
#' different stimuli and the subject
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
outdir <- "/data1/ffg05/analysis/duke/dimensionality"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/duke/dimensionality/plots"
dir.create(plotdir, showWarnings = FALSE)


#' Let's get the rois, grouping, etc information
#' 
#+ vars2 
tmp        <- load.data.std(subjects[1])
grouping   <- tmp$grouping
rois       <- tmp$rois
ys         <- tmp$yfactor
categories <- tmp$categories # this overwrites the old categories variable
roi.df$nvoxs <- table(grouping)

#' ## Data
#'  
#' Load in all the subject data
#' 
#+ data-load
dat.subs <- laply(subjects, function(subject) {
  cat(subject, "\n")
  sublst       <- load.data.std(subject)
  scale(sublst$dat)
})

#' Combine the ROIs that are bilateral into one ROI
#' 
#+ data-select
# create a new roi grouping vector
old.grouping <- grouping 
tab <- table(roi.df$names)
for (x in names(tab[tab>1])) {
  inds <- which(roi.df$names==x)
  grouping[grouping %in% inds] <- inds[1] # set everything to the first one
}
# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
# get the unique numbers
urois <- sort(unique(grouping))

#' Average the data across subjects
#' 
#+ data-average
grpdat <- apply(dat.subs, c(2,3), mean, na.rm=T)
# weird NA for 6th subject at voxel 1074...exclude from average

#' Average across subjects but shuffle the trial ordering to get random
#' data
#' 
#+ randomize-average
rand.dat.subs <- dat.subs
for (i in 1:dim(dat.subs)[1]) {
  #rand.dat.subs[i,,] <- rand.dat.subs[i,,sample(dim(dat.subs)[3])]
  rand.dat.subs[i,,] <- rand.dat.subs[i,sample(dim(dat.subs)[2]),sample(dim(dat.subs)[3])]
}
rand.grpdat <- apply(rand.dat.subs, c(2,3), mean, na.rm=T)


#' cluster the trials (rows) in each individual region to get a nicer 
#' visualization.
#' 
#' For now we only save the new ordering so we can flexibly use it later.
#' 
#+ clust-dat
clustdat <- matrix(0, nrow(grpdat), ncol(grpdat))
roi.orders <- lapply(urois, function(roi) {
  # Get subset of data associated with the ROI
  inds <- which(grouping == roi)
  tmpdat <- grpdat[,inds]
  # Compute the hierarchical clustering and get ordering of data
  d   <- dist(t(tmpdat))
  hcc <- hclust(d, method="ward.D2")
  ddc <- as.dendrogram(hcc)
  colInd <- order.dendrogram(ddc)
  # Return the ordering
  inds[colInd]
})
roi.orders <- unlist(roi.orders)
categ.orders <- lapply(categories, function(categ) {
  # Gcategbset of data associated with the ROI
  inds <- which(as.character(ys) == categ)
  tmpdat <- grpdat[inds,roi.orders]
  # Compute the hierarchical clustering and get ordering of data
  d   <- dist(tmpdat)
  hcc <- hclust(d, method="ward.D2")
  ddc <- as.dendrogram(hcc)
  rowInd <- order.dendrogram(ddc)
  # Return re-ordered data
  inds[rowInd]
})
categ.orders <- unlist(categ.orders)

#' Needed Functions
#' 
#+ functions
mytheme <- theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title = element_blank())


#' ## Correlations
#' 
#' Now we want to compute the correlations between the different voxels
#' and then average across the regions. 
#' 
#' And then we plot it!
#' 
#+ plot-corrs

# Voxel Correlations
cmat <- cor(grpdat)
# Region Correlations (averages)
n <- length(urois)
region.cmat <- matrix(0, n, n)
for (i in 1:n) {
  cat(i, "/", n, "\n")
  for (j in 1:n) {
    subset.cmat <- cmat[grouping==urois[i],grouping==urois[j]]
    region.cmat[i,j] <- mean(subset.cmat)
  }
}
names <- sub("^[RL]\ ", "", sroi.names)
dimnames(region.cmat) <- list(region1=names, region2=names)
# Convert to data frame for plotting
dfrc <- reshape2::melt(region.cmat, value.name="correlation")
dfrc$region2 <- factor(dfrc$region2, 
                       levels=rev(levels(dfrc$region2)), 
                       labels=rev(levels(dfrc$region2)))
# Plot
(
p <- ggplot(dfrc, aes(region1, region2)) + 
  geom_tile(aes(fill=correlation), colour="white", size=1) + 
  
  scale_fill_gradientn(name="Correlation", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92", 
                       limits = c(0,0.601), 
                       breaks = c(0, 0.2, 0.4, 0.6)) + 
  mytheme
)
# Save
ggsave(filename=file.path(plotdir, "region_ave_cors.pdf"), 
       plot=p, width=7, height=6)


#' We can also plot the correlation of the random data
#' 
#+ plot-rand-data

# Correlate voxels
rand.cmat <- cor(rand.grpdat)
# Collapse into region correlations
n <- length(urois)
region.rand.cmat <- matrix(0, n, n)
for (i in 1:n) {
  cat(i, "/", n, "\n")
  for (j in 1:n) {
    subset.cmat <- rand.cmat[grouping==urois[i],grouping==urois[j]]
    region.rand.cmat[i,j] <- mean(subset.cmat)
  }
}
names <- sub("^[RL]\ ", "", sroi.names)
dimnames(region.rand.cmat) <- list(region1=names, region2=names)
# Convert to data frame for plotting
dfrc <- reshape2::melt(region.rand.cmat, value.name="correlation")
dfrc$region2 <- factor(dfrc$region2, 
                       levels=rev(levels(dfrc$region2)), 
                       labels=rev(levels(dfrc$region2)))
# Plot
(
p <- ggplot(dfrc, aes(region1, region2)) + 
  geom_tile(aes(fill=correlation), colour="white", size=1) + 
  
  scale_fill_gradientn(name="Correlation", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92", 
                       limits = c(0,0.601), 
                       breaks = c(0, 0.2, 0.4, 0.6)) + 
  mytheme
)
# Save
ggsave(filename=file.path(plotdir, "random_region_ave_cors.pdf"), 
       plot=p, width=7, height=6)
