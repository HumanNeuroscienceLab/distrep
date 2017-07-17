# This script gets accuracies for each category across the
# - full model
# - reduced model
# - reduced model with the best lambda for each category

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
#' - `bar_theme`
#' 
#+ load-funs
load.data.std <- NULL; bar_theme <- NULL; summarySEwithin <- NULL
source("/data1/ffg05/scripts/repsim/scratch/stats_lab_funs.R")
source("/data1/ffg05/scripts/distrep/duke/30_classification/z_barplot_funs.R")

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
cat.cols <- brewer.pal(8, "Set2")[c(8,3,2,5,1)] # Total, Faces, Letters, Fruits, & Vehicles
#roi.cols <- brewer.pal(8, "Set3")[-c(1,7)]
roi.cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)]
# subjects
subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])
## output for saving
outdir <- "/data1/ffg05/analysis/duke/classification"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/duke/classification/plots"
dir.create(plotdir, showWarnings = FALSE)


#' ## Full Data
#' 
#' This is for the full model.
#' 
#' ### Load Data
#' 
#' Load the previously run classification analyses from each subject.
#' 
#+ full-load-data
sub.rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")

#' ### Extract Accuracy
#' 
#' We collect the accuracy for each category 
#' 
#+ full-extract-accuracy
# TODO: add the total accuracy
cat.accs <- t(sapply(sub.rets, function(x) {
  res <- x$res$alpha1$full$res
  
  accs <- res$best$category.stats$Accuracy
  accs <- as.numeric(as.character(accs))
  names(accs) <- as.character(res$best$category.stats$class)
  
  acc <- as.numeric(as.character(res$best$stats$Accuracy))
  accs <- c(Total=acc, accs)
  
  accs
}))
cat.accs <- cat.accs[,c(1,2,4,3,5)]
dimnames(cat.accs) <- list(subject=1:20, category=c("Total", categories.title))

apply(cat.accs, 2, mean)
apply(cat.accs, 2, sd)

#' ### Tests
#' 
#' FIX 
#' 
#+ full-tests
t.test(cat.accs[,1+1], cat.accs[,2+1], paired=T) # sig
t.test(cat.accs[,1+1], cat.accs[,3+1], paired=T) # sig
t.test(cat.accs[,1+1], cat.accs[,4+1], paired=T) # sig

t.test(cat.accs[,3+1], cat.accs[,2+1], paired=T) # sig
t.test(cat.accs[,3+1], cat.accs[,4+1], paired=T) # sig

t.test(cat.accs[,4+1], cat.accs[,2+1], paired=T) # not sig

#' ### Plot
#' 
#' Get that bar plot action
#' 
#+ full-plot
# convert the data to a data-frame
dfw_long <- reshape2::melt(cat.accs*100, value.name="accuracy")
dfwc <- summarySEwithin(dfw_long, measurevar="accuracy", withinvars="category",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc$sd[dfwc$category=="Total"] <- sd(cat.accs[,1]*100)
dfwc$se[dfwc$category=="Total"] <- sd(cat.accs[,1]*100)/sqrt(nrow(cat.accs))
dfwc
# plot
max.ylim <- 64
if (max(dfwc$accuracy)+max(dfwc$se) > max.ylim) stop("accuracy higher than ylim max")
# plot
p <- ggplot(dfwc, aes(x=category, y=accuracy, fill=category)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  scale_fill_manual(values=cat.cols) + 
  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank(),
    legend.position="none"
  )
# plot
plot(p)
# save
outfile <- file.path(plotdir, "category_accuracies_fullmodel_barplot.pdf")
ggsave(outfile, p, width=8, height=6)



#' ## Reduced Data
#' 
#' This is for the reduced model per subject...so the regions will vary
#' from subject to subject.
#' 
#' ### Load Data
#' 
#' We can use the already loaded data.
#' 
#' ### Extract Accuracy
#' 
#+ reduced-extract-accuracy
cat.accs <- t(sapply(sub.rets, function(x) {
  res <- x$res$alpha1$partial$res
  
  accs <- res$best$category.stats$Accuracy
  accs <- as.numeric(as.character(accs))
  names(accs) <- as.character(res$best$category.stats$class)
  
  acc <- as.numeric(as.character(res$best$stats$Accuracy))
  accs <- c(Total=acc, accs)
  
  accs
}))
cat.accs <- cat.accs[,c(1,2,4,3,5)]
dimnames(cat.accs) <- list(subject=1:20, category=c("Total", categories.title))

#' ### Tests
#' 
#+ full-tests
t.test(cat.accs[,1+1], cat.accs[,2+1], paired=T) # sig
t.test(cat.accs[,1+1], cat.accs[,3+1], paired=T) # sig
t.test(cat.accs[,1+1], cat.accs[,4+1], paired=T) # sig

t.test(cat.accs[,3+1], cat.accs[,2+1], paired=T) # sig-barely
t.test(cat.accs[,3+1], cat.accs[,4+1], paired=T) # not sig

t.test(cat.accs[,4+1], cat.accs[,2+1], paired=T) # not sig

#' ### Plot
#' 
#' Get that bar plot action
#' 
#+ full-plot
# convert the data to a data-frame
dfw_long <- reshape2::melt(cat.accs*100, value.name="accuracy")
dfwc2 <- summarySEwithin(dfw_long, measurevar="accuracy", withinvars="category",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc2$sd[dfwc2$category=="Total"] <- sd(cat.accs[,1]*100)
dfwc2$se[dfwc2$category=="Total"] <- sd(cat.accs[,1]*100)/sqrt(nrow(cat.accs))
dfwc2
# combine with the full-model
dfwc0 <- rbind(
  cbind(model="full", dfwc), 
  cbind(model="reduced", dfwc2)
)

# plot
max.ylim <- 64
if (max(dfwc0$accuracy)+max(dfwc0$se) > max.ylim) stop("accuracy higher than ylim max")
# plot
p <- ggplot(dfwc0, aes(x=category, y=accuracy, fill=category, alpha=model)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  scale_alpha_manual(values=c(1,0.7)) + 
  scale_fill_manual(values=cat.cols) + 
  coord_cartesian(ylim = c(25, max.ylim)) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank(), 
    legend.position="none"
  )
# plot
plot(p)
# save
outfile <- file.path(plotdir, "category_accuracies_full+reduced_barplot.pdf")
ggsave(outfile, p, width=8, height=6)





cat.accs <- t(sapply(1:100, function(li) { # nlambdas
  prob <- res$probs[,,li]
  pred <- res$preds[[li]]
  colnames(prob) <- levels(pred)
  acc <- category_stats(res$obs, pred, prob)$Accuracy
  as.numeric(as.character(acc))
}))
colnames(cat.accs) <- categories.title

colMeans(cat.accs)


