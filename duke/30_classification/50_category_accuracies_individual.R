
#' ## Setup
#' 
#+ setup
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
library(RColorBrewer)
library(ggplot2)

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



# simplify the roi names list as well
sroi.names <- c("FFA", "PPA", "LOC", "VWF", "V1", "M1", "A1")



#' ### Extract Accuracy
#' 
#' We collect the accuracy for each category at each individual ROI
#' 
#+ full-extract-accuracy
# TODO: add the total accuracy
cat.accs <- laply(sub.rets, function(x) {
  ares <- x$res$alpha1$individual$res
  names(ares) <- sroi.names
  
  laply(ares, function(x) {
    res <- x$res
    
    accs <- res$best$category.stats$Accuracy
    accs <- as.numeric(as.character(accs))
    names(accs) <- as.character(res$best$category.stats$class)
    
    acc <- as.numeric(as.character(res$best$stats$Accuracy))
    accs <- c(Total=acc, accs)
    
    accs
  })
})
cat.accs <- cat.accs[,,c(1,2,4,3,5)]
dimnames(cat.accs) <- list(subject=1:20, roi=sroi.names, category=c("Total", categories.title))

cat.means <- apply(cat.accs, 2:3, mean)
cat.sds   <- apply(cat.accs, 2:3, sd)
cat.sems  <- cat.sds/sqrt(20)
cat.pvals <- apply(cat.accs, 2:3, function(x) t.test(x, mu=0.25, alternative="greater")$p.value)
cat.tvals <- apply(cat.accs, 2:3, function(x) t.test(x, mu=0.25, alternative="greater")$statistic)
cat.fdrs  <- matrix(p.adjust(cat.pvals, method="fdr"), nrow(cat.pvals), ncol(cat.pvals))


colMeans(cat.means)
rowMeans(cat.means)

tmp <- (cat.fdrs < 0.05)*1
dimnames(tmp) <- dimnames(cat.tvals)
tmp # note: the average classification accuracy is significant in every region; and each category is also virtually significant. These results are very similar to the correlation based findings


# Run the ANOVA for roi x category
cat.df <- reshape2::melt(cat.accs, value.name="accuracy")
head(cat.df)
fit <- aov(accuracy ~ roi*category + Error(subject/(roi*category)), data=cat.df)
summary(fit)



#' ## Plot
#' 
#' Plot the individual category results for the individual regions.
#' 
#+ plot-individual
library(scales)
mytheme <- theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title = element_blank())
idf        <- reshape2::melt(cat.means[,-1], value.name="accuracy")
idf$pvalue <- reshape2::melt(cat.pvals[,-1])$value
idf$fdr.pvalue <- p.adjust(idf$pvalue, method="fdr")
idf$label  <- ifelse(idf$fdr.pvalue<0.05, '*', '')
idf$category    <- factor(idf$category, levels=rev(levels(idf$category)))

p <- ggplot(idf, aes(roi, category, fill = accuracy)) +
  geom_tile(colour="white", size=1) +
  geom_text(aes(label=label), size=15, colour="white") + 
  #scale_fill_gradient2(low=muted("green"), mid=muted("yellow"), high=muted("red"), midpoint=0.35) +
  #scale_fill_gradient2(low=cols[3], mid=cols[2], high=cols[1], midpoint=0.3) + 
  scale_fill_gradientn(name="Accuracy", 
                       colours=brewer.pal(9,"YlGn"), 
                       na.value="grey92", 
                       breaks=c(0.5,0.4,0.3)) + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") +
  mytheme 
plot(p)
ggsave(filename=file.path(plotdir, "individual_region+category_accuracies.pdf"), 
       plot=p, width=12, height=6)


idf$category    <- factor(idf$category, levels=rev(levels(idf$category)))
reshape2::acast(idf, roi ~ category, value.var="accuracy")

