#+ load-packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)
library(plyr)
library(RColorBrewer)
library(ggplot2)
library(niftir)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)

load.data <- NULL; wrap.model <- NULL
source("/data1/ffg05/scripts/repsim/scratch/stats_lab_funs.R")

# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")

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


### LOAD THE DATA
dat.subs <- llply(subjects, function(subject) {
  cat("subject:", subject, "\n")
  
  sublst       <- load.data(subject)
  vox.grouping <- sublst$grouping
  
  # Combine the ROIs that are bilateral into one ROI
  tab <- table(roi.df$names)
  for (x in names(tab[tab>1])) {
    inds <- which(roi.df$names==x)
    vox.grouping[vox.grouping %in% inds] <- inds[1] # set everything to the first one
  }
  # get the unique numbers
  urois <- sort(unique(vox.grouping))
  
  # And now just select the 4 regions
  urois <- urois[1:4]
  dat   <- sublst$dat[,vox.grouping %in% urois]
  vox.grouping <- vox.grouping[vox.grouping %in% urois]
  
  # Add back
  sublst$dat      <- dat
  sublst$grouping <- vox.grouping
  sublst$rois     <- urois
  
  return(sublst)
})

### COMPUTE THE COMPONENTS

#get number of components with pca and parallel-like analysis
library(MASS)
getncomps <- function(x, nsim=100, scale=FALSE, parallel=TRUE) {
  r     <- nrow(x)
  c     <- ncol(x)
  if (scale) x <- scale(x)
  true.evals <- svd(x, nu = 0)$d
  perm.evals <- laply(1:nsim, function(i) {
    y <- mvrnorm(n=r, mu=rep(0,c), Sigma=diag(1,c), empirical=F)
    evals <- svd(y, nu = 0)$d
    evals
  }, .parallel=parallel)
  pvals <- sapply(1:length(true.evals), function(ci) {
    vec <- c(true.evals[ci], perm.evals[,ci]) # combine true with permuted
    sum(vec[1]<=vec)/length(vec)
  })
  ncomps  <- sum(pvals<0.05)
  return(ncomps)
}

# compute all together
subs.ncomps <- laply(dat.subs, function(x) {
  ncomps  <- getncomps(x$dat, scale=T)
}, .progress="text")

# compute individually for each region
subs.ncomps.indiv <- laply(dat.subs, function(x) {
  ncomps  <- sapply(x$rois, function(roi) {
    getncomps(x$dat[,x$grouping==roi], scale=T)
  })
}, .progress="text")

# compute combined but only using top compenents in each region
# hmm the eigen vectors here are actually super small
subs.ncomps2 <- laply(1:length(dat.subs), function(i) {
  x <- dat.subs[[i]]
  
  combine.dat  <- llply(1:length(x$rois), function(j) {
    roi <- x$rois[[j]]
    ncomps <- subs.ncomps.indiv[i,j]
    pca.dat <- princomp(scale(x$dat[,x$grouping==roi]), scores=T)$scores
    pca.dat[,1:ncomps]
  })
  combine.dat <- do.call("cbind", combine.dat)
  
  getncomps(combine.dat, scale=T)
}, .progress="text")

# So I can now make some type of plot here
# maybe i can copy the pretty plot from the QAP poster
# I might want to also include something with how the # of features changed
library(ggplot2)
library(grid)
library(RColorBrewer)

source("duke/40_dimensionality/qa_plot_measures.R")

mat <- cbind(subs.ncomps.indiv, apply(subs.ncomps.indiv, 1, sum),  subs.ncomps)
dimnames(mat) <- list(subject=1:nrow(mat), roi=c("FFA", "PPA", "LOC", "VWFA", "Combined\nPredicted", "Combined\nActual"))
df <- reshape2::melt(mat, value.name="ncomps")

p1 <- ggplot(df, aes(x=roi, y=ncomps, color=roi, fill=roi)) +
  geom_boxplot(lwd=1.5, width=0.5) + 
  stat_summary(geom="crossbar", width=0.3, fatten=2, color="white", fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) + 
  ylab("# of Components") +
  xlab("") + 
  theme_bw() + 
  theme(panel.border=element_blank()) +
  theme(panel.grid.major.x= element_blank()) + 
  theme(panel.grid.minor  = element_blank()) + 
  theme(panel.grid.major.y= element_line(color="grey70")) + 
  theme(axis.title.x      = element_blank()) +  
  theme(axis.title.y      = element_text(family = "Century Gothic", face = "plain", 
                                         size=18, angle=90, vjust=1)) +  
  theme(axis.text.x       = element_text(family = "Century Gothic", face = "plain", 
                                         size=16, 
                                         margin=margin(0.15,0.15,0.15,0.15,"lines"))) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5, 
                                         margin=margin(0.15,0.15,0.15,0.15,"lines"))) + 
  theme(axis.ticks.y = element_line(color="grey")) +
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.x      = element_blank()) +
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines")) + 
  theme(legend.position="none")
plot(p1)
#ggsave("anat_icc_btw.png", width=4, height=3, dpi=100)

