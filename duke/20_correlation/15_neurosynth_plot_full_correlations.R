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
dir.create(file.path(corrdir, "plots"), showWarnings=F)
# Load
# note that the NULL setting is for the rstudio debugger
subjects <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL
full.grp <- NULL; full.subs <- NULL; full.detect <- NULL
load(file.path(corrdir, "neurosynth_full_correlations.rda"))
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


##' ## Select Data
##' 
##' To make the analysis more interpretable, we will be selecting the regions of
##' interest. Specifically, we will take all the right hemisphere ROIs excluding
##' the VWF.
##' 
##+ select-data
#sinds <- roi.names.title %in% c("R FFA", "R PPA", "R LOC", "L VWF", "L V1", "L #M1", "L A1")
## subject
#subset.full.subs <- full.subs[,sinds,,,]
#dimnames(subset.full.subs)$roi <- sub("[LR] ", "", roi.names.title[sinds])
## group
#subset.full.grp <- full.grp[sinds,,,]
#dimnames(subset.full.grp)$roi <- sub("[LR] ", "", roi.names.title[sinds])
## detect
#subset.full.detect <- full.detect[,sinds,]
#dimnames(subset.full.detect)$roi <- sub("[LR] ", "", roi.names.title[sinds])

#' ## Filter Data
#' 
#' To make the analysis more interpretable, we will be selecting the regions of
#' interest. Before we had only taken one hemisphere. Here we will actually be
#' averaging the results from any region with more than one hemisphere.
#' 
#+ select-data
ave.lst <- list(
  ppa = grep("PPA", roi.names.title), 
  loc = grep("LOC", roi.names.title), 
  v1 = grep("V1", roi.names.title), 
  m1 = grep("M1", roi.names.title), 
  a1 = grep("A1", roi.names.title)
)
# subjects
nrois <- length(unique(gsub("[LR] ", "", roi.names.title)))
dims  <- dim(full.subs); dims[2] <- nrois
dns   <- dimnames(full.subs); dns$roi <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
subset.full.subs <- array(NA, dims, dimnames=dns)
subset.full.subs[,1,,,] <- full.subs[,1,,,]
subset.full.subs[,2,,,] <- apply(full.subs[,ave.lst$ppa,,,], c(1,3,4,5), mean)
subset.full.subs[,3,,,] <- apply(full.subs[,ave.lst$loc,,,], c(1,3,4,5), mean)
subset.full.subs[,4,,,] <- full.subs[,6,,,]
subset.full.subs[,5,,,] <- apply(full.subs[,ave.lst$v1,,,], c(1,3,4,5), mean)
subset.full.subs[,6,,,] <- apply(full.subs[,ave.lst$m1,,,], c(1,3,4,5), mean)
subset.full.subs[,7,,,] <- apply(full.subs[,ave.lst$a1,,,], c(1,3,4,5), mean)
# group
dims  <- dim(full.grp); dims[1] <- nrois
dns   <- dimnames(full.grp); dns$roi <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
subset.full.grp <- array(NA, dims, dimnames=dns)
subset.full.grp[1,,,] <- full.grp[1,,,]
subset.full.grp[2,,,] <- apply(full.grp[ave.lst$ppa,,,], 2:4, mean)
subset.full.grp[3,,,] <- apply(full.grp[ave.lst$loc,,,], 2:4, mean)
subset.full.grp[4,,,] <- full.grp[6,,,]
subset.full.grp[5,,,] <- apply(full.grp[ave.lst$v1,,,], 2:4, mean)
subset.full.grp[6,,,] <- apply(full.grp[ave.lst$m1,,,], 2:4, mean)
subset.full.grp[7,,,] <- apply(full.grp[ave.lst$a1,,,], 2:4, mean)
# detect
dims  <- dim(full.detect); dims[2] <- nrois
dns   <- dimnames(full.detect); dns$roi <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")
subset.full.detect <- array(NA, dims, dimnames=dns)
subset.full.detect[,1,] <- full.detect[,1,]
subset.full.detect[,2,] <- apply(full.detect[,ave.lst$ppa,], c(1,3), mean)
subset.full.detect[,3,] <- apply(full.detect[,ave.lst$loc,], c(1,3), mean)
subset.full.detect[,4,] <- full.detect[,6,]
subset.full.detect[,5,] <- apply(full.detect[,ave.lst$v1,], c(1,3), mean)
subset.full.detect[,6,] <- apply(full.detect[,ave.lst$m1,], c(1,3), mean)
subset.full.detect[,7,] <- apply(full.detect[,ave.lst$a1,], c(1,3), mean)


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
gdf <- reshape2::melt(subset.full.grp[,,,1], value.name="cor")
gdf$zval <- reshape2::melt(subset.full.grp[,,,2])$value
gdf <- subset(gdf, select=c("roi", "even", "odd", "cor", "zval"))
# correct for multiple comparisons
gdf$fdr.pval <- p.adjust(pt(gdf$zval, Inf, lower.tail=F), "fdr")
gdf$thr.cor <- gdf$cor * (gdf$fdr.pval < 0.05)
gdf$thr.cor[gdf$thr.cor==0] <- NA # so the 0 box is blank when displaying plots
# reorder for better display (face first)
gdf$even <- factor(as.character(gdf$even),
                   levels=rev(levels(gdf$even)[c(1,3,2,4)]))
gdf$odd <- factor(as.character(gdf$odd),
                  levels=levels(gdf$odd)[c(1,3,2,4)])
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
sdf <- reshape2::melt(subset.full.subs[,,,,1], value.name="cor")     # Convert to a dataframe
# Fischer-z transform
sdf$fisherz <- atanh(sdf$cor)
sdf$fisherz[is.infinite(sdf$fisherz)] <- atanh(0.999)
# Z-values
sdf$zval <- reshape2::melt(subset.full.subs[,,,,2])$value
# add the category pair column
sdf$category.pair <- "different"
sdf$category.pair[sdf$even==sdf$odd] <- "same"
sdf$category.pair <- factor(sdf$category.pair)
## reorder for better display
sdf$even <- factor(as.character(sdf$even),
                   levels=rev(levels(sdf$even)[c(1,3,2,4)]))
sdf$odd <- factor(as.character(sdf$odd),
                  levels=levels(sdf$odd)[c(1,3,2,4)])
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
ddf <- reshape2::melt(subset.full.detect, value.name="accuracy")
# get the average and standard-deviation of the accuracy
# as well as the p-value using a binomial test
ddf <- ddply(ddf, .(roi, category), function(x) {
  xmean <- mean(x$accuracy)
  xsd   <- sd(x$accuracy)
  xpval <- binom.test(round(mean(x$accuracy)*20), 20, p=0.25, alternative="greater")$p.value
  data.frame(mean=xmean, sd=xsd, pvalue=xpval)
})
# correct for multiple comparisons using fdr
ddf$fdr.pvalue <- p.adjust(ddf$pvalue, method="fdr")
sum(ddf$pvalue < 0.05); sum(ddf$fdr.pvalue < 0.05)
# reorder for better visual
ddf$category <- factor(as.character(ddf$category), 
                       levels=rev(levels(ddf$category)[c(1,3,2,4)]))
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
ggsave(filename=file.path(corrdir, "plots", "fullcorr_info_same_category.pdf"), 
       plot=p, width=12, height=6)


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


#' ## All Categories
#' 
#' We want to make a nice plot with all the categories
#' 
#' ### Plot
#' 
#+ all-plot
p <- ggplot(gdf, aes(roi, even, fill = thr.cor)) +
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
ggsave(filename=file.path(corrdir, "plots", "fullcorr_info_all_pairs.pdf"), 
       plot=p, width=12, height=6)


#' ### Summary Bar Plot
#' 
#' In here, we simply want to visualize the number of significant patterns
#' for those pair of categories across regions.
#' 
#+ summary-all-plot
dat <- ddply(gdf, .(even, odd), function(x) c(num.sig=sum(x$fdr.pval<0.05)))
dat$even <- factor(as.character(dat$even),
                   levels=rev(levels(dat$even)))
cols <- brewer.pal(8, "Set2")[c(3,2,5,1)] # face, letters, fruits, car
p <- ggplot(data=dat, aes(x=even, y=num.sig, fill=odd)) +
  geom_bar(stat="identity", position=position_dodge(), colour="grey10") + 
  xlab("") + ylab("Number of Significant Regions") + 
  scale_fill_manual(values=cols) + 
  theme_minimal(14) + 
  theme(axis.ticks = element_blank(), 
        panel.grid.major.x = element_blank(), 
        axis.title = element_blank())
plot(p)
ggsave(filename=file.path(corrdir, "plots", "fullcorr_info_summary_all_pairs.pdf"), 
       plot=p, width=8, height=6)


#' ## Categorization Accuracy
#' 
#' ### ANOVA
#' 
#' Actually, let's ignore this section and focus on the plot.
#'  
#+ cat-acc-anova
dat <- reshape2::melt(subset.full.detect, value.name="accuracy")
dat$region.type <- "visual"
dat$region.type[dat$roi %in% c("M1", "A1")] <- "control"
dat$region.type <- factor(dat$region.type)
dat$category <- factor(as.character(dat$category),
                       levels=levels(dat$category)[c(1,3,2,4)])

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
xtable::xtable(mat)

#+ cat-acc-table2, results='asis'
pandoc.table(mat, 
             caption = "Shows the average categorization accuracy across subjects.")


#' ### Plot
#' 
#+ cat-acc-plot
ddf$label <- ifelse(ddf$fdr.pvalue<0.05, '*', '')
(p <- ggplot(ddf, aes(roi, category, fill = mean)) +
  geom_tile(colour="white", size=1) + 
  geom_text(aes(label=label), size=15, colour="white") + 
  scale_fill_gradientn(name="Classification Accuracy", 
                       colours=brewer.pal(9,"YlGn"), 
                       na.value="grey92") + 
  ggtitle("Comparisons Between Categories") +
  xlab("Brain Regions") +
  ylab("") + 
  mytheme)
ggsave(filename=file.path(corrdir, "plots", "fullcorr_categorization_accuracy.pdf"), 
       plot=p, width=6, height=4)
