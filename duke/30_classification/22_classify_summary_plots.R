

#' ## Setup
#' 
#+ setup
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
library(RColorBrewer)
library(ggplot2)

subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])
outdir <- "/data1/ffg05/analysis/duke/classification"
plotdir <- file.path(outdir, "plots")

roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", "R")
roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", "M1", "M1", "A1", "A1")
roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
                            fullnames=paste(roi.hemis, roi.names))
roi.df        <- roi.df[-1,]
sroi.names    <- roi.names[roi.hemis=="L"]


#' Needed functions
#' 
#+ funs
source("/data1/ffg05/scripts/distrep/duke/30_classification/z_barplot_funs.R")


#' ## Load Data
#' 
#' Load the previously run classification analyses from each subject.
#' 
#+ load-data
sub.rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")


#' ## Regions Selected
#' 
#' Let me list all the rois that were selected here

lst.rois <- llply(sub.rets, function(x) {
  ## get the roi number
  selected.rois <- x$res$alpha1$partial$rois
  ## convert to the roi name
  selected.rois <- as.character(roi.df$names[roi.df$vals%in%selected.rois])
  
  selected.rois
}) 
set.seed(42*42)
inds <- sort(sample(20, 10))
print(inds)
x <- lst.rois[inds]
names(x) <- inds
x


#' ## Accuracies
#' 
#' Collect the accuracies for the full model and the reduced models.
#' Summarize this collection by taking the mean and also the number of subjects
#' that choose X model as their best model.
#' 
#' This will get results for ridge, enet, and lasso.
#' 
#' ### Full
#' 
#+ accuracies-full
accs.full <- laply(sub.rets, function(x) {
  sapply(x$res, function(xx) xx$full$acc)
})
## summarize (lasso is always the best)
colMeans(accs.full)
apply(accs.full, 2, sd)
#' ### Reduced
#' 
#+ accuracies-reduced
accs.reduced <- laply(sub.rets, function(x) {
  sapply(x$res, function(xx) xx$partial$acc)
})
colMeans(accs.reduced)
apply(accs.reduced, 2, sd)
#' ### Individual
#' 
#+ accuracies-individual
accs.individual <- laply(sub.rets, function(x) { # subj, regions, alpha
  ret <- sapply(x$res, function(xx) xx$individual$acc)
  rownames(ret) <- sroi.names
  ret
})
apply(accs.individual, 2:3, mean)
apply(accs.individual, 2:3, sd)
mean(accs.individual[,,3]) # mean accuracy for lasso
sd(accs.individual[,,3])
apply(accs.individual, 2:3, function(x) t.test(x)$p.value)

#' Run some t-tests comparing the model accuracies
#' 
#+ accuracies-test
## full-model
wilcox.test(accs.full[,1], accs.full[,3], paired=T) # should be significant
t.test(accs.full[,1], accs.full[,3], paired=T) # should be significant
t.test(accs.full[,1], accs.full[,2], paired=T) # also significant
t.test(accs.full[,3], accs.full[,2], paired=T) # not significant
## reduced-model
wilcox.test(accs.reduced[,1], accs.reduced[,3], paired=T) # should be significant
t.test(accs.reduced[,1], accs.reduced[,3], paired=T) # should be significant
t.test(accs.reduced[,1], accs.reduced[,2], paired=T) # also significant
t.test(accs.reduced[,3], accs.reduced[,2], paired=T) # not significant


#' Convert data to data-frame for later plotting
#' 
#+ accuracies-to-df
methods <- c("Ridge", "Elastic-Net", "Lasso")
methods <- factor(methods, levels=methods, labels=methods)
## full-model
dimnames(accs.full) <- list(subject=1:nrow(accs.full), method=methods) # this is so reshape uses the correct labels
dfw_long <- reshape2::melt(accs.full*100, value.name="accuracy")
dfwc <- summarySEwithin(dfw_long, measurevar="accuracy", withinvars="method",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc
## reduced-model
dimnames(accs.reduced) <- list(subject=1:nrow(accs.reduced), method=methods) # this is so reshape uses the correct labels
dfw_long <- reshape2::melt(accs.reduced*100, value.name="accuracy")
dfwc2 <- summarySEwithin(dfw_long, measurevar="accuracy", withinvars="method",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc2
## individual-models
dimnames(accs.individual) <- list(subject=1:nrow(accs.reduced), region=sroi.names, method=methods) # this is so reshape uses the correct labels
dfw_long <- reshape2::melt(accs.individual*100, value.name="accuracy")
dfwc3 <- summarySEwithin(dfw_long, measurevar="accuracy", withinvars=c("method", "region"),
                         idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc3

#' ### Plot Accuracies
#' 
#' This will be me making those pretty barplots for publication
#' 
#' Let's start with the accuracies for the full-model.
#' 
#+ plot-acc-full-model
# y-limit
max.ylim <- 50
if (max(dfwc$accuracy)+max(dfwc$sd) > max.ylim) stop("accuracy higher than ylim max")
# plot
p <- ggplot(dfwc, aes(x=method, y=accuracy, fill=method)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank()
  )
# with the asterix for significance relative to the ridge
label.df <- data.frame(method = c("Elastic-Net", "Lasso"),
                       accuracy = c(48, 48))
p2 <- p + geom_text(data = label.df, label = "***")
# save
outfile <- file.path(plotdir, "accuracies_fullmodel_barplot.pdf")
ggsave(outfile, p2, width=8, height=6)
# plot
plot(p2)

#' Now we can do the reduced model
#' 
#+ plot-acc-reduced-model
max.ylim <- 52
if (max(dfwc2$accuracy)+max(dfwc2$sd) > max.ylim) stop("accuracy higher than ylim max")
# plot
p <- ggplot(dfwc2, aes(x=method, y=accuracy, fill=method)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank()
  )
# with the asterix for significance relative to the ridge
# oh but nothing is significant...although lasso is marginal
label.df <- data.frame(method = c("Lasso"),
                       accuracy = c(51))
p2 <- p + geom_text(data = label.df, label = "~")
# save
outfile <- file.path(plotdir, "accuracies_reducedmodel_barplot.pdf")
ggsave(outfile, p2, width=8, height=6)
# plot
plot(p2)


#' # Percent of Non-Zero Voxels
#' 
#' For the paper, we compile the percent of non-zero features for each subject
#' in the full model and the individual model
#' 
#+ percent-feats
## combined regions
perc.feats <- sapply(sub.rets, function(x) {
  nfeats <- x$res$alpha1$full$res$best$stats$nfeats
  tfeats <- sub.rets[[1]]$res$alpha1$full$res$final.fits$dim[1]
  (nfeats/tfeats)*100
})
mean(perc.feats)
sd(perc.feats)
## individual regions
perc.feats2 <- laply(sub.rets, function(x) {
  sapply(x$res$alpha1$individual$res, function(res) {
    nfeats <- res$res$best$stats$nfeats
    tfeats <- res$res$final.fits$dim[1]
    (nfeats/tfeats)*100
  })
})
mean(apply(perc.feats2, 1, mean))
sd(apply(perc.feats2, 1, mean))
## compare
t.test(perc.feats, apply(perc.feats2, 1, mean), paired=T)

##' Finally let's get the individual region model
##' 
##' Only show the results for the alpha 1 or lasso.
##' First, we get the accuracy for each individual region across the four
##' categories.
##' 
##+ plot-acc-individual
#sdf <- subset(dfwc3, method=="Lasso")
#max.ylim <- 40
#p <- ggplot(sdf, aes(x=region, y=accuracy, fill=region)) +
#  geom_bar(position="dodge", stat="identity") + 
#  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),
#                width=.1,
#                position=position_dodge(.9)) + 
#  coord_cartesian(ylim = c(25, max.ylim)) + 
#  ylab("Classification Accuracy") + 
#  bar_theme() + 
#  theme(
#    axis.title.x = element_blank()
#  )
#print(p)


#' ## Number of Features
#' 
#' Collect the number of features for the full model and the reduced models.
#' 
#+ nfeats
tmp <- laply(sub.rets, function(x) {
  sapply(x$res, function(xx) xx$full$res$best$stats$nfeats)
})
mean(tmp[,3]/tmp[,1]) # can just mention this number: 27% of the voxels ()
sd(tmp[,3]/tmp[,1])/sqrt(nrow(tmp))
mean(tmp[,3]) # 140
mean(tmp[,1]) # 510

mean(tmp[,2]/tmp[,1]) # 39%
