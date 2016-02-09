# This builds on the `20_maximal_collect_data.R` script
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
library(ggplot2)
library(RColorBrewer)
suppressMessages(library(niftir))
suppressMessages(library(doMC))
registerDoMC(20)

# load the ridge.pcor function based on the parcor package
ridge.pcor <- NULL
source("/data1/ffg05/scripts/repsim/R/pcor.R")

# Set some variables
base       <- "/data1/ffg05"
corrdir    <- file.path(base, "analysis/haxby/correlation") # output directory

# Now load the data from the previous script
# we set all the variables to be loaded to NULL for Rstudio code checker
sub.betas <- NULL; sub.zstats <- NULL; topvox.betas <- NULL; grp.topvox.betas <- NULL
runs <- NULL; runs.title <- NULL; lst.nvoxs <- NULL
categories <- NULL; categories.title <- NULL; subjects <- NULL
load(file.path(corrdir, "maximal_activity_patterns.rda"))


###
# PARTIAL CORRELATIONS VIA RIDGE REGRESSION
###

# We will be correlating each of the activation patterns between the even and
# odd runs for the same category or for different categories.
# 
# We run this first for each subject and then for the entire group.
library(abind)
# subjects
ridge.subs <- laply(topvox.betas, function(dat) {
  cat('.')
  
  res <- aaply(dat, .(1,2,3), function(sdat) {
    X <- t(sdat[1,,])
    Y <- t(sdat[2,,])
    # calculate the ridge
    res <- ridge.pcor(X, Y, scale=T, parallel=F)
    # combine results together
    res <- abind(res$cor, res$zvals, rev.along=0, new.names=c("cor", "zvals"))
    res
  }, .parallel=T)
  
  dn <- dimnames(res)
  names(dn)[4:6] <- c("even.category", "odd.category", "measure")
  dimnames(res) <- dn
  
  res # subjects x regions x category x category x measure(cor-or-zvals)
})
names(dimnames(ridge.subs))[[1]] <- "topvox"
dimnames(ridge.subs)[[1]] <- lst.nvoxs

# group
ridge.grp <- laply(grp.topvox.betas, function(dat) {
  cat('.')
  
  res <- aaply(dat, .(1,2), function(sdat) {
    X <- t(sdat[1,,])
    Y <- t(sdat[2,,])
    # calculate the ridge
    res <- ridge.pcor(X, Y, scale=T, parallel=F)
    # combine results together
    res <- abind(res$cor, res$zvals, rev.along=0, new.names=c("cor", "zvals"))
    res
  }, .parallel=T)
  
  dn <- dimnames(res)
  names(dn)[3:5] <- c("even.category", "odd.category", "measure")
  dimnames(res) <- dn
  
  res # subjects x regions x category x category x measure(cor-or-zvals)
})
names(dimnames(ridge.grp))[[1]] <- "topvox"
dimnames(ridge.grp)[[1]] <- lst.nvoxs


###
# CATEGORY DETECTION ACCURACY
###

# For each subject, we see how often an ROIs pattern of activity from the even
# run matches more with the same category pattern in the odd run and not from
# the different category pattern in the odd run. And vice-versa.

ridge.detect <- aaply(ridge.subs[,,,,,,1], .(1,2,3,4), function(cmat) {
  #cmat <- ridge.subs[1,1,,,1]
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
names(dimnames(ridge.detect))[[5]] <- "category.data"
# average across subjects (but won't save this)
ave.ridge.detect <- aaply(ridge.detect, .(1,4,5), mean)
print(round(aperm(ave.ridge.detect, c(2,3,1))*100)) 


###
# OVERAP IN ACTIVITY
###

overlap.subs <- laply(lst.nvoxs, function(topk) {
  laply(sub.zstats, function(zstats) {
    # find top voxels (maximally selective)
    top.inds.arr <- aaply(zstats, .(1,2), function(x) {
      xnew <- vector("logical", length(x))
      xnew[order(x, decreasing=T)[1:topk]] <- T
      xnew
    })
    laply(1:length(runs), function(i) {
      tcrossprod(top.inds.arr[i,,])/topk # get overlap
    })
  })
}, .parallel=TRUE)
dimnames(overlap.subs)[1:3] <- list(topvox=lst.nvoxs, subject=subjects, run=runs.title)
names(dimnames(overlap.subs))[1:3] <- c("topvox", "subject", "run")

overlap.subs.ave <- aaply(overlap.subs, .(1,4,5), mean)
print(round(aperm(overlap.subs.ave, c(2,3,1)), 3)*100) # % overlap


###
# SAVE
###

save(subjects, runs, runs.title, categories, categories.title, lst.nvoxs, 
     ridge.subs, ridge.grp, ridge.detect, overlap.subs, 
     file=file.path(corrdir, "maximal_partial_correlations.rda"))


