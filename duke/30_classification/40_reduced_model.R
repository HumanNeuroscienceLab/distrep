#' This page will deal with getting the output from the group-level analyses
#' in particular with a reduced 4 or 5 region model.
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
library(ggthemes)
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
load.data.std <- NULL; wrap.model <- NULL
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
categories <- c("faces", "letters", 'fruits', "vehicles")
categories.title <- c("Faces", "Letters", "Fruits", "Vehicles")
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
mask       <- tmp$orig.rois
grouping   <- tmp$grouping
rois       <- tmp$rois
ys         <- tmp$yfactor
categories <- tmp$categories # this overwrites the old categories variable
roi.df$nvoxs <- table(grouping)


#' ## Load Data
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
old.mask     <- mask
mask[mask==10] <- 0; mask[mask==12] <- 1; mask[mask==40] <- 6
tab <- table(roi.df$names)
for (x in names(tab[tab>1])) {
  inds <- which(roi.df$names==x)
  grouping[grouping %in% inds] <- inds[1] # set everything to the first one
  vals <- roi.df$vals[inds]
  mask[mask %in% vals] <- inds[1]
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


#' ### Run the Four Models At The Group Level
#' 
#' Now we run the ridge, elastic-net, and lasso on the group model. Here we run
#' the full, partial, and individual region models...just cuz.
#' 
#+ run-models
sgrpdat         <- scale(grpdat)
subset.urois    <- urois[1:4]
subset.grouping <- grouping[grouping %in% subset.urois]
subset.grpdat   <- sgrpdat[,grouping %in% subset.urois]
subset.roi.orders <- roi.orders[grouping %in% subset.urois]
alpha <- 1

#rgrp.res4 <- wrap.rcfe(subset.grpdat, ys, foldsI, subset.grouping, 
#                       alpha=alpha, nlambda=100)
grp.res4 <- wrap.model(subset.grpdat, ys, foldsI, alpha=alpha, nlambda=100)

#' Save the model output
#' 
#+ save-reduced
outfile <- file.path(outdir, "group_4region_model_alpha1.rda")
#grp.res4 <- rgrp.res4$full$res
save(grp.res4, file=outfile)
#rm(grp.res4)
#load(outfile)


#' #### Plot
#' 
#' Now get the full-model and plot
#' 
#+ plot-reduced-model
res <- grp.res4
## TODO: add the accuracy for each category...but need this per subject
betas <- res$best$betas
rownames(betas) <- categories.title
## remove the 0s from the plot
betas[betas==0] <- NA
## reorder the categories (faces, letters, fruits, vehicles)
## reorder the columns to clump the rois and related voxels
betas      <- betas[4:1, subset.roi.orders]
## plot with only positive/negative
x <- (betas>0) - (betas<0)
## get the colors
dat.cols <- rev(c(brewer.pal(3, "Set1")[1], "white", brewer.pal(3, "Set1")[2]))
tmp <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)][urois %in% subset.urois]
cc <- as.character(factor(grouping[subset.roi.orders], levels=subset.urois, labels=tmp))
## plot
outfile <- file.path(plotdir, sprintf("group_4region_fullmodel_alpha%s.pdf", as.character(alpha)))
pdf(outfile, width=4, height=4, paper='special')
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols)
dev.off()
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols)



#' ### Best Model Per Category
#' 
#' We get an individual lambda for each category that gets the best classification
#' performance.
#' 
#+ model-per-category
# First, we get the category statistics at each lambda
res <- grp.res4
retcats <- t(sapply(1:100, function(li) { # nlambdas
  prob <- res$probs[,,li]
  pred <- res$preds[[li]]
  colnames(prob) <- levels(pred)
  category_stats(res$obs, pred, prob)$Accuracy
}))
colnames(retcats) <- categories.title
plot.ts(retcats)
#plot.ts(retcats, plot.type = "single", col=1:4)
# Second we find the best lambda and get its betas for each category
max.retcats <- apply(retcats, 2, which.max)
betas <- t(sapply(1:length(max.retcats), function(i) {
  res$final.fits$beta[[i]][,max.retcats[i]]
}))
rownames(betas) <- categories.title
## remove the 0s from the plot
betas[betas==0] <- NA
## reorder the categories (faces, letters, fruits, vehicles)
## reorder the columns to clump the rois and related voxels
betas      <- betas[4:1, subset.roi.orders]
## plot with only positive/negative
x <- (betas>0) - (betas<0)
## get the colors
#dat.cols <- c("blue", "white", "red")
dat.cols <- rev(c(brewer.pal(3, "Set1")[1], "white", brewer.pal(3, "Set1")[2]))
#tmp <- c("#943156", "#AA7539", "#679933", "#26596A", brewer.pal(10, "Set3")[c(6,7,10)])
tmp <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
cc <- as.character(factor(grouping[subset.roi.orders], levels=urois, labels=tmp))
## fix for vehicles (remember to remove if it shows up in the plot)
x[1,1] <- -1
## plot + save
outfile <- file.path(plotdir, sprintf("group_4region_best_alpha%s.pdf", as.character(alpha)))
pdf(outfile, width=4, height=4, paper='special')
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, col=dat.cols)
dev.off()
## just plot
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, col=dat.cols)


#' ### Accuracies by Lambda
#' 
#' We can get each categories accuracy by lambda and plot it.
#' 
#+ plot-accuracies-by-lambda
res <- grp.res4
cat.accs <- t(sapply(1:100, function(li) { # nlambdas
  prob <- res$probs[,,li]
  pred <- res$preds[[li]]
  colnames(prob) <- levels(pred)
  acc <- category_stats(res$obs, pred, prob)$Accuracy
  acc <- as.numeric(as.character(acc))
  acc
}))
dimnames(cat.accs) <- list(lamdba=res$lambdas, category=categories.title)
## for pasting in keynote
write.table(cat.accs, row.names=F, col.names=F, quote=F)
write.table(round(res$lambdas, 4), row.names=F, col.names = F, quote=F)

#' We can also plot our results by the df instead of lambda.
#' 
#' For this, we need to average across any times there are multiple 
#' values of df for a given category.
#' 
#+ plot-by-df
## get the dfs
dfmat <- cbind(average=res$final.fits$df/4, t(res$final.fits$dfmat))
cat.accs2 <- cbind(total=res$stats$Accuracy, cat.accs)
titles <- c("Average", categories.title)
## collapse duplicate dfs together by averaging
df.df <- ldply(1:5, function(ci) {
  cat.df <- data.frame(df=dfmat[,ci], acc=cat.accs2[,ci])
  cat.df <- ddply(cat.df, .(df), function(x) {
    c(accuracy=mean(x$acc))
  })
  data.frame(category=titles[ci], cat.df)
})
## get the maximal dfs + accuracies for each category
max.linds <- apply(cat.accs2, 2, which.max)
max.df <- ldply(1:5, function(i) {
  data.frame(category=titles[i], 
    df=as.numeric(dfmat[max.linds[i],i]), 
    accuracy=as.numeric(cat.accs2[max.linds[i],i]))
})
## plot
p <- ggplot(df.df, aes(x=df, y=accuracy, color=category)) + 
  #geom_vline(xintercept=max.df$df[1], linetype="dashed") + 
  geom_line(size=1.5) + 
  geom_point(data=max.df, size=4) + 
  scale_color_manual(values=c("grey",cat.cols)) + 
  theme_gdocs()
plot(p)
## save
outfile <- file.path(plotdir, "category_accuracies_by_df.pdf")
ggsave(outfile, p, width=8, height=6)


#' Number of Features in Each Region
#' 
#+ nfeats-regions
fac <- factor(subset.grouping, levels=unique(subset.grouping), labels=sroi.names[1:4])
nfeats.all <- laply(1:4, function(fi) {
  betas <- res$final.fits$beta[[fi]]
  nfeats.pos <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x>0)))
  nfeats.neg <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x<0)))
  abind::abind(pos=nfeats.pos, neg=nfeats.neg, along = 3)
})
nfeats.all <- abind::abind(apply(nfeats.all, 2:4, mean), nfeats.all, along=1)
rev.cats <- c("Average",categories.title)
dimnames(nfeats.all) <- list(category=rev.cats, region=sroi.names[1:4], 
                             df=res$final.fits$df, direction=c("pos","neg"))
#dimnames(nfeats.all) <- list(category=rev.cats, region=sroi.names[1:4], 
#                             df=1/res$lambdas, direction=c("pos","neg"))
df2 <- reshape2::melt(nfeats.all, value.name="nfeats")
## collapse duplicate dfs together by averaging
df.df0 <- ldply(1:5, function(ci) {
  cat.df <- subset(df2, category==rev.cats[ci])
  cat.df <- ddply(cat.df, .(category, region, direction, df), function(x) {
    c(nfeats=mean(x$nfeats))
  })
  cat.df
})
df.df <- subset(df.df0, category!="Average")
## this is if you want to collapse across df for each category
#df.df0 <- ldply(1:5, function(ci) {
#  cat.df <- subset(df2, category==rev.cats[ci])
#  cat.df <- ddply(cat.df, .(region, direction), function(x) {
#    x$df <- dfmat[,ci]
#    x
#  })
#  cat.df <- ddply(cat.df, .(category, region, direction, df), function(x) {
#    c(nfeats=mean(x$nfeats))
#  })
#  cat.df
#})
#df.df <- subset(df.df0, category!="Average")
## other approach
#dimnames(nfeats.all) <- list(category=rev.cats, region=sroi.names[1:4], 
#                             df=res$lambdas, direction=c("pos","neg"))
#df.df0 <- reshape2::melt(nfeats.all, value.name="nfeats")
#df.df <- subset(df.df0, category!="Average")
## smooth
#df.df <- ddply(df.df, .(category, region, direction), function(dat) {
#  newx <- smooth.spline(as.numeric(dat$df), as.numeric(dat$nfeats), df=25)
#  data.frame(df=dat$df, nfeats=dat$nfeats, smooth.nfeats=rev(newx$y))
#})
cat.cols <- brewer.pal(8, "Set2")[c(8,3,2,5,1)] # Total, Faces, Letters, Fruits, & Vehicles
dfs <- res$final.fits$df
sel.df <- dfs[which.max(rowMeans(cat.accs))]
p <- ggplot(df.df, aes(x=df, y=nfeats, color=category)) + 
  geom_vline(xintercept=sel.df, linetype="dashed") + 
  geom_line(size=1.5) + 
  #geom_point(data=max.df, size=4) + 
  scale_color_manual(values=cat.cols[-1]) + 
  theme_gdocs() + 
  #scale_x_reverse() + 
  facet_grid(direction~region)
plot(p)
## save
outfile <- file.path(plotdir, "category_nfeats_by_df_01.pdf")
ggsave(outfile, p, width=8, height=6)


sel.df2 <- dfs[apply(cat.accs, 2, which.max)]
sel.df2 <- data.frame(category=categories.title, df=sel.df2)
p <- ggplot(df.df, aes(x=df, y=nfeats, color=region)) + 
  geom_vline(aes(xintercept=df), linetype="dotted", data=sel.df2) + 
  geom_vline(xintercept=sel.df, linetype="dashed") + 
  geom_line(size=1.5) + 
  #geom_point(data=max.df, size=4) + 
  scale_color_manual(values=roi.cols) + 
  theme_gdocs() + 
  #scale_x_reverse() + 
  facet_grid(direction~category)
plot(p)
## save
outfile <- file.path(plotdir, "category_nfeats_by_df_02.pdf")
ggsave(outfile, p, width=8, height=6)


#' ## Probabilities
#' 
#' Let's get the trial category probabilities plotted
res <- grp.res4
lind <- res$best$stats$lind
probs <- t(res$probs[,,lind])
rownames(probs) <- categories.title
probs <- probs[4:1,]
#cols      <- rev(brewer.pal(11, "RdBu"))
cols      <- rev(brewer.pal(11, "Spectral"))
## column colors
cc <- as.character(factor(res$obs, levels=levels(res$obs), labels=cat.cols))
## just plot
heatmap((probs>0.5)*probs, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, col=cols)
heatmap(probs, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, col=cols, zlim=c(0,1), scale="none")
tmp <- probs*(probs>0.5)
#tmp[tmp==0] <- NA
heatmap(tmp, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, col=cols, zlim=c(0,1), scale="none")
## plot + save
outfile <- file.path(plotdir, sprintf("group_4region_best_alpha%s.pdf", as.character(alpha)))
pdf(outfile, width=4, height=4, paper='special')
heatmap(probs, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, col=cols)
dev.off()





#' ### Other Test
#' 
#' I'm going to now try with repeats on the fold
#' 
#+ run-models2
## get the folds again
nrepeats <- 10
runFolds <- createMultiFolds(uruns, k=5, times = nrepeats)
folds    <- lapply(runFolds, function(x) {
  which(runs %in% x)
})
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})
## run the new model
grp.res4.5s <- llply(seq(1,length(foldsI),by=5), function(i) {
  wrap.model(subset.grpdat, ys, foldsI[i:(i+4)], alpha=alpha, nlambda=100, parallel=T)
}, .parallel=F, .progress="text")
## see how it compares
tmp2 <- laply(grp.res4.5s, function(res) {
  retcats <- t(sapply(1:100, function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    category_stats(res$obs, pred, prob)$Accuracy
  }))
  colnames(retcats) <- categories.title
  print(apply(retcats, 2, which.max))
  apply(retcats, 2, as.numeric)
})


#' Okay that didn't work
#' 
#' so we also try running a separate model for each category
#' 
#+ run-models3
nrepeats <- 1
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})

ret <- llply(levels(ys), function(fac) {
  ys2 <- factor((ys==fac)*1, levels=c(0,1), labels=c(paste("Not", fac, sep=""), fac))
  tmp <- wrap.model(subset.grpdat, ys2, foldsI, alpha=alpha, nlambda=100, family="multinomial")
  tmp
}, .progress="text")
names(ret) <- levels(ys)

ldply(levels(ys), function(fac) {
  tmp <- ret[[fac]]
  li <- which.max(tmp$stats$Kappa)
  cbind(fac=fac, tmp$stats[li,])
})

ldply(levels(ys), function(fac) {
  tmp <- ret[[fac]]
  li <- which.max(tmp$stats$ROC)
  cbind(fac=fac, tmp$stats[li,])
})

ldply(levels(ys), function(fac) {
  tmp <- ret[[fac]]
  cbind(fac=fac, tmp$best$stats)
})




grp.res42 <- wrap.model(subset.grpdat, ys, foldsI, alpha=alpha, nlambda=50)
res <- grp.res42

retcats <- t(sapply(1:length(res$lambdas), function(li) { # nlambdas
  prob <- res$probs[,,li]
  pred <- res$preds[[li]]
  colnames(prob) <- levels(pred)
  category_stats(res$obs, pred, prob)$Accuracy
}))
colnames(retcats) <- categories.title
retcats <- apply(retcats, 2, as.numeric)
print(apply(retcats, 2, which.max))


cat.accs <- t(sapply(1:length(res$lambdas), function(li) { # nlambdas
  prob <- res$probs[,,li]
  pred <- res$preds[[li]]
  colnames(prob) <- levels(pred)
  acc <- category_stats(res$obs, pred, prob)$Accuracy
  acc <- as.numeric(as.character(acc))
  acc
}))
dimnames(cat.accs) <- list(lamdba=res$lambdas, category=categories.title)

tmpfun1 <- function(res) {
  cat.accs <- t(sapply(1:length(res$lambdas), function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    acc <- category_stats(res$obs, pred, prob)$Accuracy
    acc <- as.numeric(as.character(acc))
    acc
  }))
  dimnames(cat.accs) <- list(lamdba=res$lambdas, category=categories.title)
  dfmat <- cbind(average=res$final.fits$df/4, t(res$final.fits$dfmat))
  cat.accs2 <- cbind(total=res$stats$Accuracy, cat.accs)
  titles <- c("Average", categories.title)
  ## collapse duplicate dfs together by averaging
  df.df <- ldply(1:5, function(ci) {
    cat.df <- data.frame(df=dfmat[,ci], acc=cat.accs2[,ci])
    cat.df <- ddply(cat.df, .(df), function(x) {
      c(accuracy=mean(x$acc))
    })
    data.frame(category=titles[ci], cat.df)
  })
  ## get the maximal dfs + accuracies for each category
  max.linds <- apply(cat.accs2, 2, which.max)
  max.df <- ldply(1:5, function(i) {
    data.frame(category=titles[i], 
               df=as.numeric(dfmat[max.linds[i],i]), 
               accuracy=as.numeric(cat.accs2[max.linds[i],i]))
  })
  
  max.df
}


# Lastly try some type of permutation
alphas <- seq(0,1,by=0.25)
wrp.grp.res4 <- llply(alphas, function(alpha) {
  wrap.model(subset.grpdat, ys, foldsI, alpha=alpha, nlambda=100)
}, .progress="text")
names(wrp.grp.res4) <- alphas

llply(alphas, function(alpha) {
  rdf <- tmpfun1(wrp.grp.res4[[as.character(alpha)]])
  cbind(alpha=alpha, rdf)
})

grp.res4.p0 <- wrap.model(subset.grpdat, ys, foldsI, alpha=alpha, nlambda=100)
grp.res4.z0 <- wrap.model(subset.grpdat, ys, foldsI, alpha=0.5, nlambda=100)

grp.res4.p1 <- wrap.model(subset.grpdat, sample(ys), foldsI, alpha=alpha, nlambda=50)
grp.res4.p2 <- wrap.model(subset.grpdat, sample(ys), foldsI, alpha=alpha, nlambda=50)
grp.res4.p3 <- wrap.model(subset.grpdat, sample(ys), foldsI, alpha=alpha, nlambda=50)
grp.res4.p4 <- wrap.model(subset.grpdat, sample(ys), foldsI, alpha=alpha, nlambda=50)

tmpfun1(grp.res4.p0)
tmpfun1(grp.res4.p1)
tmpfun1(grp.res4.p2)
tmpfun1(grp.res4.p3)
tmpfun1(grp.res4.p4)

tmpfun1(grp.res4.z0)


tmpfun0 <- function(res) {
  cat.accs <- t(sapply(1:length(res$lambdas), function(li) { # nlambdas
    prob <- res$probs[,,li]
    pred <- res$preds[[li]]
    colnames(prob) <- levels(pred)
    acc <- category_stats(res$obs, pred, prob)$Accuracy
    acc <- as.numeric(as.character(acc))
    acc
  }))
  dimnames(cat.accs) <- list(lamdba=res$lambdas, category=categories.title)
  apply(cat.accs, 2, which.max)
}
