#' ---
#' title: "Category Information with Neurosynth ROIs in the Duke Dataset"
#' author: "Zarrar Shehzad"
#' date: "February 8th, 2016"
#' output: 
#'   html_document:
#'     fig_caption: true
#'     toc: true
#'     number_sections: true
#' ---

#' Our goal for this script is to generate plots comparing the patterns of 
#' activity between pairs of categories for every neurosynth based ROI.
#'
#' We generate outputs for the 'group-average' thresholded for significant 
#' results as well as the percent correct category detection measure.
#'
#' Note that this builds on `12_neurosynth_partial_correlations.R`.
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
corrdir    <- file.path(base, "analysis/duke/correlation") # output directory
# Load
# note that the NULL setting is for the rstudio debugger
subjects <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL
ridge.grp <- NULL; ridge.subs <- NULL; ridge.detect <- NULL
load(file.path(corrdir, "neurosynth_partial_correlations.rda"))
# Setup other relevant variables
#roi.vals      <- c(10, 12, 20, 22, 30, 32, 40, 50, 52, 60, 62, 70, 72)
#roi.hemis     <- c("L", "R", "L", "R", "L", "R", "L", "L", "R", "L", "R", "L", #"R")
#roi.names     <- c("FFA", "FFA", "PPA", "PPA", "LOC", "LOC", "VWF", "V1", "V1", #"M1", "M1", "A1", "A1")
#roi.df        <- data.frame(vals=roi.vals, hemi=roi.hemis, names=roi.names, 
#                            fullnames=paste(roi.hemis, roi.names))
#roi.df        <- roi.df[-1,] # remove L FFA since it's small


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
#' To make the analysis more interpretable, we will be selecting the regions of
#' interest. Specifically, we will take all the right hemisphere ROIs excluding
#' the VWF.
#' 
#+ select-data
sinds <- roi.names.title %in% c("R FFA", "R PPA", "R LOC", "L VWF", "L V1", "L M1", "L A1")
# subject
subset.ridge.subs <- ridge.subs[,sinds,,,]
dimnames(subset.ridge.subs)$roi <- sub("[LR] ", "", roi.names.title[sinds])
# group
subset.ridge.grp <- ridge.grp[sinds,,,]
dimnames(subset.ridge.grp)$roi <- sub("[LR] ", "", roi.names.title[sinds])
# detect
subset.ridge.detect <- ridge.detect[,sinds,]
dimnames(subset.ridge.detect)$roi <- sub("[LR] ", "", roi.names.title[sinds])


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
gdf <- subset(gdf, select=c("roi", "even", "odd", "cor", "zval"))
# correct for multiple comparisons
gdf$fdr.pval <- p.adjust(pt(gdf$zval, Inf, lower.tail=F), "fdr")
gdf$thr.cor <- gdf$cor * (gdf$fdr.pval < 0.05)
gdf$thr.cor[gdf$thr.cor==0] <- NA # so the 0 box is blank when displaying plots
# reorder for better display (face first)
gdf$even <- factor(as.character(gdf$even),
                  levels=rev(levels(gdf$even)))
# add the category pair column
gdf$category.pair <- "different"
gdf$category.pair[gdf$even==gdf$odd] <- "same"
gdf$category.pair <- factor(gdf$category.pair)
## reorder for better display
gdf$category.pair <- factor(as.character(gdf$category.pair),
                            levels=rev(levels(gdf$category.pair)))
# add the region type column
gdf$region.type <- "visual"
gdf$region.type[gdf$roi %in% c("M1", "A1")] <- "control"
gdf$region.type <- factor(gdf$region.type)
## reorder for better display
gdf$region.type <- factor(as.character(gdf$region.type),
                          levels=rev(levels(gdf$region.type)))
# rearrange columns
gdf <- gdf[,c("region.type", "roi", "category.pair", "even", "odd", 
            "cor", "zval", "fdr.pval", "thr.cor")]

#' ### Subject Results
#' 
#+ format-data-subjects
sdf <- reshape2::melt(subset.ridge.subs[,,,,1], value.name="cor")     # Convert to a dataframe
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
sdf$category.pair <- factor(as.character(sdf$category.pair),
                            levels=rev(levels(sdf$category.pair)))
# add the region type column
sdf$region.type <- "visual"
sdf$region.type[sdf$roi %in% c("M1", "A1")] <- "control"
sdf$region.type <- factor(sdf$region.type)
## reorder for better display
sdf$region.type <- factor(as.character(sdf$region.type),
                          levels=rev(levels(sdf$region.type)))
# reorder columns
sdf <- sdf[,c("subjects", "region.type", "roi", "category.pair", "even", "odd", 
              "cor", "fisherz", "zval")]
head(sdf)

#' ### Category Detection Results
#' 
#+ format-data-detect
ddf <- reshape2::melt(subset.ridge.detect, value.name="accuracy")
# get the average and standard-deviation of the accuracy
# as well as the p-value using a binomial test
ddf <- ddply(ddf, .(roi, category), function(x) {
  xmean <- mean(x$accuracy)
  xsd   <- sd(x$accuracy)
  xpval <- binom.test(sum(x$accuracy), 20, p=0.25, alternative="greater")$p.value
  data.frame(mean=xmean, sd=xsd, pvalue=xpval)
})
# correct for multiple comparisons using fdr
ddf$fdr.pvalue <- p.adjust(ddf$pvalue, method="fdr")
sum(ddf$pvalue < 0.05); sum(ddf$fdr.pvalue < 0.05)
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
dat$category <- dat$even        # Category column (even/odd same)
# Focus for now on visual areas
dat <- subset(dat, region.type=="visual")
dat$roi <- factor(dat$roi)
# Select columns
dat <- subset(dat, select=c("subjects", "roi", "category", "cor", "fisherz", "zval"))
# Run the anova but only within the visual areas
fit <- aov(fisherz ~ roi*category + Error(subjects), data=dat)
res <- summary(fit)
# for plotting in latex
xtable::xtable(res, digits=c(0,0,2,2,1,3))

#' We can see that the main effects of ROI and category are significant 
#' but the interaction effect.
#+ unique-information-anova-table, results='asis'
panderOptions("digits", 2)
pander(res)
# maybe for later?
# ddply(sdf, .(roi), colwise(function(x) tanh(mean(x)), .(fisherz)))
# ddply(sdf, .(roi), colwise(function(x) tanh(sd(x)), .(cor)))


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
names(dat)[4] <- "category"
dat <- dat[,-5]
# finally plot
p <- ggplot(dat, aes(roi, category, fill = thr.cor)) +
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradientn(name="Correlation", 
                       colours=brewer.pal(9,"YlOrRd"), 
                       na.value="grey92") + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") + facet_grid(.~region.type, scales="free_x") + 
  mytheme #+ 
#theme(axis.text.x = element_blank())
plot(p)

#' ### Visual vs Control Regions for Same Category
#' 
#' Let's also look at the comparison between the control regions and visual
#' regions.
#' 
#+ unique-info-visual-vs-control
dat <- subset(sdf, category.pair=="same")    # Select only same category pairs
dat$category <- dat$even        # Category column (even/odd same)
# Select columns
dat <- subset(dat, select=c("subjects", "region.type", "roi", "category", "cor", "fisherz", "zval"))
# ANOVA result
fit <- aov(fisherz ~ region.type*category + Error(subjects), data=dat)
res <- summary(fit)
# for plotting in latex
xtable::xtable(res, digits=c(0,0,2,2,1,3))

#+ unique-info-visual-vs-control-table, results='asis'
panderOptions("digits", 2)
pander(res)




#' ## Categorization Accuracy
#' 
#' ### ANOVA
#' 
#' Actually, let's ignore this section and focus on the plot.
#'  
#+ cat-acc-anova
dat <- reshape2::melt(subset.ridge.detect, value.name="accuracy")
dat$region.type <- "visual"
dat$region.type[dat$roi %in% c("M1", "A1")] <- "control"
dat$region.type <- factor(dat$region.type)
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
mean.arr <- reshape2::acast(ddf, roi ~ category, value.var="mean")
pval.arr <- reshape2::acast(ddf, roi ~ category, value.var="fdr.pvalue")
mat <- matrix(sprintf("%i%s", mean.arr*100, gsub(" ", "", add.significance.stars(pval.arr))), 
              nrow(mean.arr), ncol(mean.arr), dimnames=dimnames(mean.arr))

#+ cat-acc-table2, results='asis'
pandoc.table(mat, 
             caption = "Shows the average categorization accuracy across subjects.")

