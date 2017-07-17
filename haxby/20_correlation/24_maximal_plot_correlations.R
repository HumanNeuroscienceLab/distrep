#' ---
#' title: "Category Information with Maximal Responses in the Haxby Dataset"
#' author: "Zarrar Shehzad"
#' date: "February 15th, 2016"
#' output: 
#'   html_document:
#'     fig_caption: true
#'     toc: true
#'     number_sections: true
#' ---

#' Our goal for this script is to generate plots comparing the patterns of 
#' activity between pairs of categories for every maximal response. That is 
#' we get the maximal response in the even or odd dataset for faces and each of
#' the other categories, and then apply each of those masks to the data for 
#' faces and the other categories.
#'
#' We generate outputs for the 'group-average' thresholded for significant 
#' results as well as the percent correct category detection measure.
#'
#' Note that this builds on `12_maximal_partial_correlations.R`.
#'
#' ## Load
#' 
#' ### Packages
#' 
#' using dev version of ggplot??
#' install_github("hadley/scales")
#' install_github("hadley/ggplot2")
#' 
#+ setup
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
library(ggplot2)
library(RColorBrewer)
suppressMessages(require(grid))
suppressMessages(require(gridExtra))
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)
library(pander) # for tables


#' ### Data and Variables
#' 
#' We load the data and some relevant variables related to the ROIs
#' 
#+ load
# Set some variables first
base       <- "/data1/ffg05"
corrdir    <- file.path(base, "analysis/haxby/correlation") # output directory
dir.create(file.path(corrdir, "plots"), showWarnings=F)
# Load
# note that the NULL setting is for the rstudio debugger
subjects <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL
ridge.grp <- NULL; ridge.subs <- NULL; ridge.detect <- NULL
load(file.path(corrdir, "maximal_partial_correlations.rda"))


#' ## Functions
#' 
#' Any functions that we need. Here only the generic plot theme.
#' 
#+ functions
mytheme <- theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.title = element_blank())


#' ## Select Data
#' 
#' Select only the 30 voxel data.
#' 
#+ select-data
# Get only 30 voxels
# Collapse across even/odd run masks
## grp
subset.ridge.grp <- ridge.grp[3,,,,,]
subset.ridge.grp <- apply(subset.ridge.grp, 2:5, mean)
## subs
subset.ridge.subs <- ridge.subs[3,,,,,,]
subset.ridge.subs <- apply(subset.ridge.subs, c(1,3:6), mean)
## detect
subset.ridge.detect <- ridge.detect[3,,,,]
subset.ridge.detect <- apply(subset.ridge.detect, c(1,3:4), mean)


#' ## Format Data
#' 
#' ### Group Average Results
#' 
#' We want to format our partial correlations into a dataframe for easier 
#' analysis. In the process, we also need to correct for multiple comparisons 
#' and add two new variables: 
#' - category.pair (same or different)
#' - region.type (visual or control)
#' 
#+ format-data-group
# Convert to a dataframe
gdf <- reshape2::melt(subset.ridge.grp[,,,1], value.name="cor")
gdf$zval <- reshape2::melt(subset.ridge.grp[,,,2])$value
colnames(gdf) <- c("category.mask", "even", "odd", "cor", "zval")
# correct for multiple comparisons
gdf$fdr.pval <- p.adjust(pt(gdf$zval, Inf, lower.tail=F), "fdr")
gdf$thr.cor <- gdf$cor * (gdf$fdr.pval < 0.05)
gdf$thr.cor[gdf$thr.cor==0] <- NA # so the 0 box is blank when displaying plots
# reorder for better display (face first)
gdf$category.mask <- factor(as.character(gdf$category.mask),
                            levels=rev(levels(gdf$category.mask)))
gdf$even <- factor(as.character(gdf$even),
                   levels=rev(levels(gdf$even)))
#gdf$odd <- factor(as.character(gdf$odd),
#                  levels=levels(gdf$odd))
# add the category pair column
gdf$category.pair <- "different"
gdf$category.pair[gdf$even==gdf$odd] <- "same"
gdf$category.pair <- factor(gdf$category.pair)
## reorder for better display
gdf$category.pair <- factor(as.character(gdf$category.pair),
                            levels=rev(levels(gdf$category.pair)))
# rearrange columns
gdf <- gdf[,c("category.mask", "category.pair", "even", "odd", 
              "cor", "zval", "fdr.pval", "thr.cor")]


#' ### Subject Results
#' 
#+ format-data-subjects
sdf <- reshape2::melt(subset.ridge.subs[,,,,1], value.name="cor")     # Convert to a dataframe
colnames(sdf)[3:4] <- c("even", "odd")
# Fischer-z transform
sdf$fisherz <- atanh(sdf$cor)
sdf$fisherz[is.infinite(sdf$fisherz)] <- atanh(0.999)
# Z-values
sdf$zval <- reshape2::melt(subset.ridge.subs[,,,,2])$value
# add the category pair column
sdf$category.pair <- "different"
sdf$category.pair[sdf$even==sdf$odd] <- "same"
sdf$category.pair <- factor(sdf$category.pair)
## reorder for better display
sdf$category.mask <- factor(as.character(sdf$category.mask),
                            levels=rev(levels(sdf$category.mask)))
sdf$even <- factor(as.character(sdf$even),
                   levels=rev(levels(sdf$even)))
#sdf$odd <- factor(as.character(sdf$odd),
#                  levels=levels(sdf$odd))
sdf$category.pair <- factor(as.character(sdf$category.pair),
                            levels=rev(levels(sdf$category.pair)))
# reorder columns
sdf <- sdf[,c("subject", "category.mask", "category.pair", "even", "odd", 
              "cor", "fisherz", "zval")]
head(sdf)


#' ### Category Detection Results
#' 
#+ format-data-detect
ddf <- reshape2::melt(subset.ridge.detect, value.name="accuracy")
# get the average and standard-deviation of the accuracy
# as well as the p-value using a binomial test
ddf <- ddply(ddf, .(category.mask, category.data), function(x) {
  xmean <- mean(x$accuracy)
  xsd   <- sd(x$accuracy)
  xpval <- binom.test(round(mean(x$accuracy)*20), 20, p=0.25, alternative="greater")$p.value
  data.frame(mean=xmean, sd=xsd, pvalue=xpval)
})
# correct for multiple comparisons using fdr
ddf$fdr.pvalue <- p.adjust(ddf$pvalue, method="fdr")
sum(ddf$pvalue < 0.05); sum(ddf$fdr.pvalue < 0.05)
# reorder for better visual
#ddf$category.mask <- factor(as.character(ddf$category.mask), 
#                            levels=levels(ddf$category.mask))
#ddf$category.data <- factor(as.character(ddf$category.data), 
#                            levels=levels(ddf$category.data))
# ok so let's see what we've got
head(ddf)


#' ## Unique Information
#' 
#' First we analyze and plot the amount of unique information (correlation) for
#' each category.
#' 
#' ### Visual Areas for Same Categories
#' 
#' #### ANOVA
#' 
#' Our first analysis is an ANOVA to see if the correlations vary
#' by region. If the information is evenly distributed, then it shouldn't 
#' matter which region we sample.
#' 
#+ unique-information-anova
# Setup the dataframe
dat <- subset(sdf, category.pair=="same")    # Select only same category pairs
dat$category.data <- dat$even        # Category column (even/odd same)
# Select columns
dat <- subset(dat, select=c("subject", "category.mask", "category.data", "cor", "fisherz", "zval"))

#' #### Plots
#' 
#' We can now plot the information for each same category pair. Note that here
#' we are using the group-average correlation results that I have calculated. 
#' Also note that we will plot both the visual and the control regions here.
#' 
#' First, we plot the within-category information for each region and category.
#' 
#+ unique-information-plot
# select only same category pairs
dat <- subset(gdf, category.pair=="same")
names(dat)[3] <- "category.data"
dat <- dat[,-4]
# collapse 
dat$category.data <- factor(as.character(dat$category.data), 
                            levels=rev(levels(dat$category.data)))
# finally plot
p <- ggplot(dat, aes(category.data, category.mask, fill = thr.cor)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradientn(name="Correlation", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92") + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") + #facet_grid(.~region.type, scales="free_x") + 
  mytheme #+ 
#theme(axis.text.x = element_blank())
plot(p)
ggsave(filename=file.path(corrdir, "plots", "maximal_unique_info_same_category.pdf"), 
       plot=p, width=12, height=10)

# wondering if this table might do a better job?
mat <- reshape2::acast(category.mask ~ category.data, data=dat, value.var='thr.cor')[8:1,]
print(round(mat, 2))

# also wondering about making a plot that takes the averages from each row/column
mat[is.na(mat)] <- 0
mean.catdata <- colMeans(mat, na.rm=T)
mean.catmask <- rowMeans(mat, na.rm=T)
round(mean.catdata, 2)
round(mean.catmask, 2)


#' ## All Categories
#' 
#' TODO: see what they do in order to the small objects mask!
#' 
#' We want to make a nice plot with all the categories
#' 
#' ### Plot
#' 
#+ all-plot
dat <- gdf
dat$category.mask <- factor(as.character(dat$category.mask), 
                            levels=rev(levels(dat$category.mask)))
p <- ggplot(dat, aes(category.mask, even, fill = thr.cor)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradientn(name="Correlation", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92") + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") + facet_grid(odd~.) + 
  mytheme #+ 
#theme(axis.text.x = element_blank())
plot(p)
ggsave(filename=file.path(corrdir, "plots", "maximal_unique_info_all_pairs.pdf"), 
       plot=p, width=12, height=12)


#' ### Summary Bar Plot
#' 
#' In here, we simply want to visualize the number of significant patterns
#' for those pair of categories across regions.
#' 
#+ summary-all-plot
dat <- ddply(gdf, .(even, odd), function(x) c(num.sig=sum(x$fdr.pval<0.05)))
dat$even <- factor(as.character(dat$even),
                   levels=rev(levels(dat$even)))
#cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # face, letters, fruits, car
p <- ggplot(data=dat, aes(x=even, y=num.sig, fill=odd)) +
  geom_bar(stat="identity", position=position_dodge(), colour="grey10") + 
  xlab("") + ylab("Number of Significant Regions") + 
  #scale_fill_manual(values=cols) + 
  theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title = element_blank())
plot(p)
#dat$odd <- factor(as.character(dat$odd),
#                   levels=rev(levels(dat$odd)))
#p <- ggplot(dat, aes(even, odd, fill = factor(num.sig, levels=8:0))) +
#  geom_tile(colour="white", size=0.5) +
#  scale_fill_manual(name="# of Sig", 
#                       values=rev(c("grey92", brewer.pal(8,"YlOrRd")))) + 
#  ggtitle("Comparisons Between Categories") +
#  xlab("Brain Regions") +
#  ylab("") + 
#  mytheme
#plot(p)
ggsave(filename=file.path(corrdir, "plots", "maximal_unique_info_summary_all_pairs.pdf"), 
       plot=p, width=8, height=6)



#' ## Categorization Accuracy
#' 
#' ### ANOVA
#' 
#' Actually, let's ignore this section and focus on the plot.
#'  
#+ cat-acc-anova
dat <- reshape2::melt(subset.ridge.detect, value.name="accuracy")
#dat$category.mask <- factor(as.character(dat$category.mask),
#                            levels=levels(dat$category.mask))
#dat$category.data <- factor(as.character(dat$category.data),
#                            levels=levels(dat$category.data))

# Run the anova but only within the visual areas
#fit <- glm(accuracy ~ roi*category + Error(subjects), family=binomial(), data=dat)
#fit <- aov(accuracy ~ roi*category + Error(subjects), data=subset(dat, region.type=="visual"))
#res <- summary(fit)
#res
# for plotting in latex
#xtable::xtable(res, digits=c(0,0,2,2,1,3))


#' ### Table
#' 
#' This should be the table showing all that we need.
#' 
#+ cat-acc-table
mean.arr <- reshape2::acast(ddf, category.mask ~ category.data, value.var="mean")
sd.arr <- reshape2::acast(ddf, category.mask ~ category.data, value.var="sd")
pval.arr <- reshape2::acast(ddf, category.mask ~ category.data, value.var="fdr.pvalue")
mat <- matrix(sprintf("%i%s", round(mean.arr*100), gsub(" ", "", add.significance.stars(pval.arr))), 
              nrow(mean.arr), ncol(mean.arr), dimnames=dimnames(mean.arr))
xtable::xtable(mat)

#+ cat-acc-table2, results='asis'
pandoc.table(mat, 
             caption = "Shows the average categorization accuracy across subjects.")

#' ### Plot
#' 
#+ cat-acc-plot
ddf$label <- ifelse(ddf$fdr.pvalue<0.05, '*', '')
ddf$category.mask <- factor(ddf$category.mask, 
                            levels=rev(levels(ddf$category.mask)), 
                            labels=rev(levels(ddf$category.mask)))
(p <- ggplot(ddf, aes(category.data, category.mask, fill = mean)) +
  geom_tile(colour="white", size=1) + 
  geom_text(aes(label=label), size=12, colour="white") + 
  scale_fill_gradientn(name="Classification Accuracy", 
                       colours=brewer.pal(9,"YlGn"), 
                       na.value="grey92") + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") + 
  mytheme)
ggsave(filename=file.path(corrdir, "plots", "maximal_categorization_accuracy.pdf"), 
       plot=p, width=6, height=5)




#' ## All voxels
#' 
#' We will plot the pattern similarity for each of the different voxel size set.
#' 
#+ all-voxels
dat <- apply(ridge.grp, c(1,3:6), mean)
dat <- dat[,,,,1] # select the cor values
dat <- reshape2::melt(dat, value.name="cor")
dat <- dat[dat$even==dat$odd,] # only keep same
dat <- dat[,-4]
colnames(dat) <- c("topvox", "category.mask", "category.data", "cor")
dat <- subset(dat, topvox!=30)
#cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # face, letters, fruits, car
p <- ggplot(data=dat, aes(x=topvox, y=cor, color=category.data)) +
  geom_line() + 
  facet_wrap(~category.mask, ncol=2) + 
  xlab("Most Selective Voxels") + ylab("Partial Correlation") + 
  theme_minimal(14) + 
  #scale_color_manual(values=cols) + 
  scale_x_log10(breaks=lst.nvoxs[-3]) + 
  theme(axis.ticks = element_blank(), 
        axis.title = element_blank())
plot(p)
ggsave(filename=file.path(corrdir, "plots", "maximal_corr_vox_lineplots.pdf"), 
       plot=p, width=8, height=6)

dat <- apply(ridge.detect, c(1,4:5), mean)
dat <- reshape2::melt(dat, value.name="cor")
dat <- subset(dat, topvox!=30)
p <- ggplot(data=dat, aes(x=topvox, y=cor, color=category.data)) +
  geom_line() + 
  facet_wrap(~category.mask, ncol=2) + 
  xlab("Most Selective Voxels") + ylab("Partial Correlation") + 
  theme_minimal(14) + 
  #scale_color_manual(values=cols) + 
  scale_x_log10(breaks=lst.nvoxs[-3]) + 
  theme(axis.ticks = element_blank(), 
        axis.title = element_blank())
plot(p)
ggsave(filename=file.path(corrdir, "plots", "maximal_categorization_vox_lineplots.pdf"), 
       plot=p, width=8, height=6)





