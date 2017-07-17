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
outdir <- "/data1/ffg05/analysis/duke/classification"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/duke/classification/plots"
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


## TEST
#library(randomForest)
#ret <- randomForest(grpdat, y=ys, ntree=125, importance=T, do.trace=F)
#plot.ts(ret$err.rate[,2])
#head(ret$importance)
#?randomForest
#
#rets <- llply(urois, function(ri) {
#  randomForest(grpdat[,grouping!=ri], y=ys, ntree=125, importance=T, do.trace=F)
#}, .parallel=T)
#
#tmp <- sapply(rets, function(x) x$err.rate[125,1])
#names(tmp) <- sroi.names
#round((ret$err.rate[125,1]-tmp)*100*-1, 1)
#
#round(laply(urois, function(ri) colMeans(ret$importance[grouping==ri,1:4]))*1000, 2)


#' ## Classification
#' 
#' ### Cross-Validation Setup
#' 
#' We will be doing a leave-one run out cross-validation. Since we have 10 runs
#' of data for each subject, this could also be considered a 10-fold cross-
#' validation. Presentation of each stimulus from each category is more or less
#' evenly distributed across the runs. Here we will get the indices for each 
#' cross-validation fold.
#' 
#+ cross-validation-setup
nrepeats <- 1
runs     <- tmp$runs # since all subjects the same can use runs from 1 subject
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
folds    <- lapply(runFolds, function(x) {
  which(runs %in% x)
})
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})


#' ## Clustered Data
#' 
#' We will cluster the voxels (columns) in each individual region as well as 
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


#' ### Models
#' 
#' Now we run the ridge, elastic-net, and lasso on the group model. Here we run
#' the full, partial, and individual region model.
#' 
#+ run-models
alphas <- c(0,0.5,1)
nlambda <- 100
grp.res <- llply(alphas, function(alpha) {
  cat("\nAlpha", as.character(alpha), "\n")
  wrap.rcfe(grpdat, ys, foldsI, grouping, alpha=alpha, nlambda=nlambda)
})
names(grp.res) <- sprintf("alpha%s", as.character(alphas))


#' #### Full Model
#' 
#' Let's examine the results of the full model right now.
#' 
#+ full-models
accs <- sapply(grp.res, function(x) x$full$acc)
print(accs)
# we can see from the accs that the elastic-net performs the best

#' Let's save these models for future use (might want to move the plot below as well)
#' 
#+ save-full-models
outfile <- file.path(outdir, "group_fullmodels.rda")
save(grp.res, grouping, roi.orders, categ.orders, categories.title, 
     sroi.names, roi.df, ys, file=outfile)


#' Now we can plot the non-zero betas for the best performing model.
#' To make it easier to see, we will be showing only if the betas are
#' positive or negative.
#' 
#+ plot-full-models
for (alpha in names(grp.res)) {
  cat(alpha, "\n")
  betas <- grp.res[[alpha]]$ful$res$best$betas
  rownames(betas) <- categories.title[c(1,3,2,4)]
  
  ## remove the 0s from the plot
  betas[betas==0] <- NA
  ## reorder the categories (faces, letters, fruits, vehicles)
  ## reorder the columns to clump the rois and related voxels
  betas      <- betas[4:1,roi.orders]
  ## plot with only positive/negative
  x <- (betas>0) - (betas<0)
  ## get the colors
  #dat.cols <- c("blue", "white", "red")
  dat.cols <- rev(c(brewer.pal(3, "Set1")[1], "white", brewer.pal(3, "Set1")[2]))
  #tmp <- c("#943156", "#AA7539", "#679933", "#26596A", brewer.pal(10, "Set3")[c(6,7,10)])
  tmp <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
  cc <- as.character(factor(grouping[roi.orders], levels=urois, labels=tmp))
  ## plot
  heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
          col=dat.cols)
  ## save
  outfile <- file.path(plotdir, sprintf("group_fullmodel_%s.pdf", alpha))
  pdf(outfile, width=4, height=4, paper='special')
  heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
          col=dat.cols)
  dev.off()
}

#' We can also plot the individual results for just the FFA as a comparison
#' 
#+ res
print(grp.res$alpha1$individual[[1]]$res$best$category.stats)
print(grp.res$alpha1$individual[[1]]$res$best$stats)

betas <- grp.res$alpha1$individual[[1]]$res$best$betas
rownames(betas) <- categories.title[c(1,3,2,4)]

## remove the 0s from the plot
betas[betas==0] <- NA
## reorder the categories (faces, letters, fruits, vehicles)
## reorder the columns to clump the rois and related voxels
betas      <- betas[4:1,roi.orders[grouping==1]]
## plot with only positive/negative
x <- (betas>0) - (betas<0)
## get the colors
dat.cols <- rev(c(brewer.pal(3, "Set1")[1], "white", brewer.pal(3, "Set1")[2]))
tmp <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
cc <- rep(tmp[1], ncol(betas))
## plot
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols)
## save
outfile <- file.path(plotdir, "ffa_fullmodel_alpha1.pdf")
pdf(outfile, width=4, height=4, paper='special')
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols)
dev.off()

#' We can now do the same for each of the individual plots
#' 
#+ res2
names(grp.res$alpha1$individual) <- sroi.names
betas <- llply(grp.res$alpha1$individual, function(res) {
  betas <- res$res$best$betas
  betas
})
betas <- do.call(cbind, betas)
rownames(betas) <- categories.title#[c(1,3,2,4)]
## remove the 0s from the plot
betas[betas==0] <- NA
## reorder the categories (faces, letters, fruits, vehicles)
## reorder the columns to clump the rois and related voxels
betas      <- betas[4:1,roi.orders]
## plot with only positive/negative
x <- (betas>0) - (betas<0)
## get the colors
dat.cols <- rev(c(brewer.pal(3, "Set1")[1], "white", brewer.pal(3, "Set1")[2]))
tmp <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
cc <- as.character(factor(grouping[roi.orders], levels=urois, labels=tmp))
## plot
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols)
## save
outfile <- file.path(plotdir, "group_fullmodel_individaul_a1.pdf")
pdf(outfile, width=6, height=6, paper='special')
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols)
dev.off()
