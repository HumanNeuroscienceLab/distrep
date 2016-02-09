# This builds on the `10_neurosynth_collect_data.R` script
# which collected all the activity patterns together. Here we 
# compute the partial correlations using the linear ridge regression.
# The correlations are between the activity patterns of pairs of cateogiries 
# for a given ROI.

###
# SETUP
###

# Load packages
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(plyr)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(24)

# load the ridge.pcor function based on the parcor package
ridge.pcor <- NULL
source("/data1/ffg05/scripts/repsim/R/pcor.R")

# Set some variables
base       <- "/data1/ffg05"
corrdir    <- file.path(base, "analysis/duke/correlation") # output directory

# Now load the data from the previous script
# we set all the variables to be loaded to NULL for Rstudio code checker
lldat <- NULL; cldat <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL; subjects <- NULL
load(file.path(corrdir, "neurosynth_activity_patterns.rda"))


###
# PARTIAL CORRELATIONS VIA RIDGE REGRESSION
###

# We will be correlating each of the activation patterns between the even and
# odd runs for the same category or for different categories.
# 
# We run this first for each subject and then for the entire group.
library(abind)
# subjects
ridge.subs <- laply(lldat, function(ldat) {
  cat('.')
  laply(ldat, function(dat) {
    X <- dat[,,1]
    Y <- dat[,,2]
    # calculate the ridge
    res <- ridge.pcor(X, Y, scale=T, parallel=F)
    # combine results together
    res <- abind(res$cor, res$zvals, rev.along=0, new.names=c("cor", "zvals"))
    res
  }, .parallel=T)
}) # subjects x regions x category x category x measure(cor-or-zvals)
dimnames(ridge.subs) <- list(subjects=subjects, roi=roi.names.title, 
                             even=categories.title, odd=categories.title, 
                             measure=c("cor", "zvals"))
# groups
ridge.grp <- laply(cldat, function(dat) {
  X <- dat[,,1]
  Y <- dat[,,2]
  # calculate the ridge
  res <- ridge.pcor(X, Y, scale=T, parallel=F)
  # combine results together
  res <- abind(res$cor, res$zvals, rev.along=0, new.names=c("cor", "zvals"))
  res
}, .parallel=T)
dimnames(ridge.grp) <- list(roi=roi.names.title, 
                            even=categories.title, odd=categories.title, 
                            measure=c("cor", "zvals"))


###
# CATEGORY DETECTION ACCURACY
###

# For each subject, we see how often an ROIs pattern of activity from the even
# run matches more with the same category pattern in the odd run and not from
# the different category pattern in the odd run. And vice-versa.

ridge.detect <- aaply(ridge.subs[,,,,1], .(1,2), function(cmat) {
  #cmat <- ridge.subs[1,1,,,1]
  # NOTE: even/odd should be same because cmat is a correlation matrix
  even_detect <- laply(1:nrow(cmat), function(ri) {
    all(cmat[ri,ri] > cmat[ri,-ri])*1
  })
  odd_detect <- laply(1:ncol(cmat), function(ci) {
    all(cmat[ci,ci] > cmat[-ci,ci])*1
  })
  detect <- (even_detect + odd_detect)/2
  names(detect) <- categories.title
  detect
}, .parallel=TRUE)
names(dimnames(ridge.detect))[3] <- "category"
# average across subjects (but won't save this)
ave.ridge.detect <- aaply(ridge.detect, .(2,3), mean)
print(round(ave.ridge.detect*100)) 

# this is the category detection with the betas 
# (so we need to recompute the partial correlations)
ridge.detect.betas <- laply(lldat, function(ldat) {
  cat('.')
  laply(ldat, function(dat) {
    X <- dat[,,1]
    Y <- dat[,,2]
    # calculate the ridge
    res <- ridge.pcor(X, Y, scale=T, parallel=F, extra.output=T)
    # take only the betas
    betas <- res$betas
    # get the detection
    even_detect <- laply(1:nrow(betas), function(ri) {
      all(betas[ri,ri] > betas[ri,-ri])*1
    })
    odd_detect <- laply(1:ncol(betas), function(ci) {
      all(betas[ci,ci] > betas[-ci,ci])*1
    })
    detect <- (even_detect + odd_detect)/2
    names(detect) <- categories.title
    detect
  }, .parallel=T)
}) # subjects x regions x category x category x measure(cor-or-zvals)
dimnames(ridge.detect.betas) <- list(subject=subjects, roi=roi.names.title, 
                                     category=categories.title)
# average but again don't save ## interesting that this is lower
ave.ridge.detect <- aaply(ridge.detect.betas, .(2,3), mean)
print(round(ave.ridge.detect*100)) 


###
# SAVE
###

save(subjects, roi.names, roi.names.title, 
     categories, categories.title, 
     ridge.grp, ridge.subs, ridge.detect, ridge.detect.betas, 
     file=file.path(corrdir, "neurosynth_partial_correlations.rda"))
