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

#' # Setup
#' 
#' ## Variables
#' 
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
roi.cols <- brewer.pal(8, "Set3")[-c(1,7)]

subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])

#' Let's get the rois, grouping, etc information
#+ vars2 
tmp <- load.data.std(subjects[1])
grouping   <- tmp$grouping
rois       <- tmp$rois
ys         <- tmp$yfactor
categories <- tmp$categories # this overwrites the old categories variable
roi.df$nvoxs <- table(grouping)


#' ## Data
#'  
#' Load in all the subject data
#+ load-subs
dat.subs <- laply(subjects, function(subject) {
  cat(subject, "\n")
  sublst       <- load.data.std(subject)
  scale(sublst$dat)
})

#' Select the right hemisphere stuff + left VWF
#' 
#+ select-data
srois <- c("R FFA", "R PPA", "R LOC", "L VWF", "R V1", "R M1", "R A1")
inds <- grouping %in% which(roi.names.title %in% srois)
dat.subs <- dat.subs[,,inds]
grouping <- grouping[inds]

#' Average the data across subjects
#' 
#+ ave-subs
grpdat <- apply(dat.subs, c(2,3), mean, na.rm=T)
# weird NA for 6th subject at voxel 1074...exclude from average
grpdat[,1074] <- colMeans(dat.subs[-6,,1074])


#' # Classification
#' 
#' ## Cross-Validation Setup
#' 
#' We will be doing a leave-one run out cross-validation. Since we have 10 runs
#' of data for each subject, this could also be considered a 10-fold cross-
#' validation. Presentation of each stimulus from each category is more or less
#' evenly distributed across the runs. Here we will get the indices for each 
#' cross-validation fold.
#' 
#+ cross-validation-setup
nrepeats <- 1
runs     <- tmp$runs
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
folds    <- lapply(runFolds, function(x) {
  which(runs %in% x)
})
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})

#' ## Full Model Classification
#' 
#' I'll start by applying the classification to the full dataset to show case
#' how our implementation of the elastic-net works. Our goal is to select 
#' voxels that are predictive of category membership. Using the elastic-net, we
#' will be selecting sets of correlated voxels.
#' 
#' Also to start, we use the lasso and show the stats for the best model
#' as well as the stats for the best model for each category.
#' 
#+ full-classify
tmp <- aperm(dat.subs, c(2,1,3))
dim(tmp) <- c(prod(dim(tmp)[1:2]), dim(tmp)[3])
all.equal(as.numeric(tmp[,1]), as.numeric(unlist(llply(1:20, function(i) dat.subs[i,,1]))))
tmp[is.na(tmp)] <- 0

ret <- wrap.model(tmp, rep(ys, 20), foldsI, alpha=1, nlambda=100)

ret <- wrap.model(grpdat, ys, foldsI, alpha=1, nlambda=100)
print(subset(ret$best$stats, select=c("alpha", "lambda", "Accuracy", 
                                      "AccuracyPValue", "ROC", "logLoss", 
                                      "nfeats")))
print(ret$best$category.stats)
print(ret$best$overlap)

# Plot it as a heatmap
ret <- res.rfces$alpha0.5$full$res
## betas
betas <- ret$best$betas
betas[betas==0] <- NA
rownames(betas) <- categories.title
## reorder so regions are grouped
oinds      <- order(grouping)
betas      <- betas[,oinds]
rgrouping  <- grouping
## reorder the categories (faces, letters, fruits, vehicles)
betas      <- betas[c(4,2,3,1),]
## get the colols
dat.cols <- rev(brewer.pal(11, "RdYlBu"))
cc <- as.character(factor(rgrouping, levels=sort(unique(rgrouping)), 
                          labels=brewer.pal(7, "Dark2")))
## plot
heatmap(betas, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols, zlim=c(-1,1))
## plot with only positive/negative
x <- (betas>0) - (betas<0)
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=c("blue", "white","red"))






alphas <- c(0,0.5,1)
nlambda <- 100
res.rfces <- llply(alphas, function(alpha) {
  cat("\nAlpha", as.character(alpha), "\n")
  wrap.rcfe(grpdat, ys, foldsI, grouping, alpha=alpha, nlambda=nlambda)
})
names(res.rfces) <- sprintf("alpha%s", as.character(alphas))




# Try wtih full data
res2 <- wrap.rcfe(grpdat, ys, foldsI, grouping, alpha=0.5, nlambda=nlambda)

# Plot it as a heatmap
ret <- res2$full$res
## betas
betas <- ret$best$betas
betas[betas==0] <- NA
rownames(betas) <- categories.title
## reorder so regions are grouped
#rgrouping  <- grouping[grouping %in% res2$reduced$urois]
#oinds      <- order(rgrouping)
oinds      <- order(grouping)
betas      <- betas[,oinds]
rgrouping  <- grouping[oinds]

## reorder the categories (faces, letters, fruits, vehicles)
betas      <- betas[c(4,2,3,1),]
## get the colols
dat.cols <- rev(brewer.pal(11, "RdYlBu"))
cc <- as.character(factor(rgrouping, levels=sort(unique(rgrouping)), 
                          labels=brewer.pal(12, "Set3")))
## plot
heatmap(betas, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols, zlim=c(-1,1))
## plot with only positive/negative
x <- (betas>0) - (betas<0)
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=c("blue", "white","red"))




#' The results above appear to suggest even at the group level localization for
#' category information but a localization that is distributed?

#' We can also visualize the percent of voxels that are used in our lasso model.
#' 
#+ full-classify-viz
# note that the betas here doesn't have the intercept
# we first those voxels that have non-zero features for any category
nonzeros <- apply(ret$best$betas!=0, 2, any)
# then we can make a table with the percent of features needed in each region
tab <- tapply(nonzeros, sublst$grouping, mean)
names(tab) <- roi.names.title
barplot(tab)

ret <- res.rfces$alpha0.5$reduced$res
dat.cols <- rev(brewer.pal(11, "RdYlBu"))
cat.cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # Faces, Letters, Fruits, & Vehicles
roi.cols <- brewer.pal(8, "Set3")[-c(1,7)]
rc   <- as.character(factor(ys, levels=categories, labels=cat.cols))
cc   <- as.character(factor(grouping, levels=rois, labels=roi.cols))
x <- t(ret$best$betas>0) - t(ret$best$betas<0)
heatmap(t(x), 
        Rowv=NA, Colv=NA, labRow=NA, labCol=NA, 
        col = c("blue", "white","red"), 
        RowSideColors = unique(rc), ColSideColors = cc)

# todo: try to use ggplot to make this type of plot with boxes and all
# 
# we can also now save this data to a file in order to visualize this on the
# brain
outdir <- "/data1/ffg05/analysis/classification/group"
infile <- "/data1/ffg05/analysis/rois/group/neurosynth_visual_rois_3mm.nii.gz"
hdr    <- read.nifti.header(infile)
mask   <- read.mask(infile, NULL) > 10
tmpr   <- read.mask("/data1/ffg05/analysis/rois/group/neurosynth_visual_rois_3mm.nii.gz", NULL)
tmpr   <- tmpr[mask]
for (i in 1:length(categories)) {
  categ <- categories[i]
  # we need to reorder the data back to what it was
  os       <- order(tmpr)
  rdat     <- vector("numeric", ncol(ret$best$betas))
  rdat[os] <- ret$best$betas[i,]
  # save
  cat("Saving", categ, "\n")
  write.nifti(rdat, hdr, mask, 
              outfile=file.path(outdir, sprintf("avebetas_lasso_%s.nii.gz", categ)), 
              overwrite=T)
}


tmp <- laply(1:4, function(i) tapply(ret$best$betas[i,], grouping, function(x) sum(x!=0)))
tmp <- as.matrix(tmp)
dimnames(tmp) <- list(categ=categories, roi=roi.names.title)
tmp

print(ret$best$category.stats)

# Apply the recursive feature thingy
rcfe.ret <- wrap.rcfe(grpdat, ys, foldsI, grouping, alpha=1, nlambda=100)

accs  <- c(rcfe.ret$full$acc, rcfe.ret$reduced$acc, sapply(rcfe.ret$individual, function(x) x$acc))
names(accs) <- c("Full", "Reduced", roi.names)
barplot(accs, 
        main="Classification Accuracy", 
        sub=sprintf("Reduced: %s", paste(rcfe.ret$reduced$urois, collapse=", ")))

