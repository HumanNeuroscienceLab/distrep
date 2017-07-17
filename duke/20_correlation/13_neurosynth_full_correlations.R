# This builds on the `10_neurosynth_collect_data.R` script
# which collected all the activity patterns together. Here we 
# compute the partial correlations using the linear regression.
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

# Set some variables
base       <- "/data1/ffg05"
corrdir    <- file.path(base, "analysis/duke/correlation") # output directory

# Now load the data from the previous script
# we set all the variables to be loaded to NULL for Rstudio code checker
lldat <- NULL; cldat <- NULL; roi.names <- NULL; roi.names.title <- NULL
categories <- NULL; categories.title <- NULL; subjects <- NULL
load(file.path(corrdir, "neurosynth_activity_patterns.rda"))


###
# Full Correlations with z-scores
###

t2z <- function(t, kappa) {
  t <- as.matrix(t)
  z <- matrix(0, nrow(t), ncol(t))
  z[t>0] <- qt(pt(t[t>0], kappa-1, lower.tail=F), Inf, lower.tail=F)
  z[t<0] <- qt(pt(t[t<0], kappa-1, lower.tail=T), Inf, lower.tail=T)
  z
}

lm.fcor <- function(X, Y, parallel=F) {
  # Get the correlations
  full.cors <- cor(X, Y)
  # Get the betas (inefficient i know)
  full.betas <- laply(1:ncol(X), function(i) {
    fit  <- lm(Y ~ X[,i])
    sfit <- summary(fit)
    betas <- laply(sfit, function(x) {
      sfit[[1]]$coefficients[2,1]
    })
    betas
  }, .parallel=parallel)
  # Get the zvals
  zvals <- laply(1:ncol(X), function(i) {
    fit  <- lm(Y ~ X[,i])
    sfit <- summary(fit)
    zvals <- laply(sfit, function(x) {
      t2z(x$coefficients[2,3], x$df[2])
    })
    zvals
  }, .parallel=parallel)
  # average
  full.zvals <- (zvals + t(zvals))/2
  # get in the names
  dimnames(full.zvals) <- dimnames(full.cors)
  dimnames(full.betas) <- dimnames(full.betas)
  
  list(cors=full.cors, zvals=full.zvals, betas=full.betas)
}

# We will be correlating each of the activation patterns between the even and
# odd runs for the same category or for different categories.
# 
# We run this first for each subject and then for the entire group.
library(abind)
# subjects
full.subs <- laply(lldat, function(ldat) {
  cat('.')
  laply(ldat, function(dat) {
    X <- dat[,,1]
    Y <- dat[,,2]
    # calculate the linear regression / correlation
    res <- lm.fcor(X, Y, parallel=F)
    # combine results together
    res <- abind(res$cor, res$zvals, rev.along=0, new.names=c("cor", "zvals"))
    res
  }, .parallel=T)
}) # subjects x regions x category x category x measure(cor-or-zvals)
dimnames(full.subs) <- list(subjects=subjects, roi=roi.names.title, 
                             even=categories.title, odd=categories.title, 
                             measure=c("cor", "zvals"))
# groups
full.grp <- laply(cldat, function(dat) {
  X <- dat[,,1]
  Y <- dat[,,2]
  # calculate the linear regression / correlation
  res <- lm.fcor(X, Y, parallel=F)
  # combine results together
  res <- abind(res$cor, res$zvals, rev.along=0, new.names=c("cor", "zvals"))
  res
}, .parallel=T)
dimnames(full.grp) <- list(roi=roi.names.title, 
                            even=categories.title, odd=categories.title, 
                            measure=c("cor", "zvals"))


###
# CATEGORY DETECTION ACCURACY
###

# For each subject, we see how often an ROIs pattern of activity from the even
# run matches more with the same category pattern in the odd run and not from
# the different category pattern in the odd run. And vice-versa.

full.detect <- aaply(full.subs[,,,,1], .(1,2), function(cmat) {
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
names(dimnames(full.detect))[3] <- "category"
# average across subjects (but won't save this)
ave.full.detect <- aaply(full.detect, .(2,3), mean)
print(round(ave.full.detect*100)) 

# this is the category detection with the betas 
# (so we need to recompute the partial correlations)
full.detect.betas <- laply(lldat, function(ldat) {
  cat('.')
  laply(ldat, function(dat) {
    X <- dat[,,1]
    Y <- dat[,,2]
    # calculate the ridge
    res <- lm.fcor(X, Y, parallel=F)
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
dimnames(full.detect.betas) <- list(subject=subjects, roi=roi.names.title, 
                                     category=categories.title)
# average but again don't save ## interesting that this is lower
ave.full.detect <- aaply(full.detect.betas, .(2,3), mean)
print(round(ave.full.detect*100)) 


###
# SAVE
###

save(subjects, roi.names, roi.names.title, 
     categories, categories.title, 
     full.grp, full.subs, full.detect, full.detect.betas, 
     file=file.path(corrdir, "neurosynth_full_correlations.rda"))
