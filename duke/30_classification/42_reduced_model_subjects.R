#' This script will run the subject level analysis

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
#+ load-funs
load.data <- NULL; wrap.model <- NULL
source("/data1/ffg05/scripts/repsim/scratch/stats_lab_funs.R")


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
outdir <- "/data1/ffg05/analysis/duke/classification"
dir.create(outdir, showWarnings = FALSE)
plotdir <- "/data1/ffg05/analysis/duke/classification/plots"
dir.create(plotdir, showWarnings = FALSE)

# simplify the roi names list as well
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")



#' ## Load Data
#' 
#' Now we load the data
#' 
#+ load-data
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

## redo roi.df
roi.df <- roi.df[dat.subs[[1]]$rois,]
roi.df$fullnames <- sroi.names[1:4]



#' ## Run
#' 
#' Run the reduced model on all the subjects
#' 
#+ run-models
lst.res <- llply(1:length(dat.subs), function(si) {
  # Get runs for this person
  runs <- dat.subs[[si]]$runs
  
  # Get the folds
  nrepeats <- 1
  uruns    <- sort(unique(runs))
  nruns    <- length(uruns)
  runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
  foldsI   <- lapply(runFolds, function(x) {
    runs %in% x
  })
  
  # Scale the data
  sdat     <- scale(dat.subs[[si]]$dat)
  
  # Loop through each alpha
  obs      <- dat.subs[[si]]$yfactor
  rois     <- dat.subs[[si]]$grouping
  
  # Run the model!!!
  res      <- wrap.model(sdat, obs, foldsI, alpha=1, nlambda=100, parallel=F)
  res      <- save.vars(res, rois)
  
  res
}, .parallel=T)

#' Also run individual ones again?

#' Let's save these models for future use
#' 
#+ save-models
outfile <- file.path(outdir, "subjects_4region_models.rda")
save(lst.res, dat.subs, roi.df, categories.title, file=outfile)




#' ## AUC along path
#' 
#' ### Load Data
#+ sub-std-data
dat.subs.std <- laply(subjects, function(subject) {
  cat(subject, "\n")
  sublst       <- load.data.std(subject)
  scale(sublst$dat)
})
sublst       <- load.data.std(subjects[2])
grouping     <- sublst$grouping
nrepeats <- 1
runs     <- sublst$runs


#' ### Classify
#' 
#+ sub-std-classify
# Get the folds
nrepeats <- 1
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})
# Get the classification in standard space
lst.sres <- llply(1:length(subjects), function(si) {
  # Loop through each alpha
  sdat     <- dat.subs.std[si,,grouping %in% c(1,2,4,6)]
  if (si==6) { # one subject has a column with NAs
    w <- apply(sdat, 2, function(x) any(is.na(x)))
    sdat[,w] <- rowMeans(sdat[,!w])
  }
  obs      <- sublst$yfactor
  # Run the model!!!
  res      <- wrap.model(sdat, obs, foldsI, alpha=1, nlambda=100, parallel=T)
  res
}, .parallel=F, .progress="text")

#' ### Num of Feats by Varying Lambda
#' 
#' First, let's get all the nfeats.
#' 
#+ feats-by-lambda0
inds     <- grouping %in% c(1,2,4,6)
fac <- factor(grouping[inds], levels=unique(grouping[inds]), labels=sroi.names[1:4])
subs.nfeats.all <- ldply(1:20, function(si) {
  res <- lst.sres[[si]]
  nfeats.all <- laply(1:4, function(fi) {
    betas <- res$final.fits$beta[[fi]]
    nfeats.pos <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x>0)))
    nfeats.neg <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x<0)))
    abind::abind(pos=nfeats.pos, neg=nfeats.neg, along = 3)
  })
  dimnames(nfeats.all) <- list(category=categories.title, region=sroi.names[1:4], 
                               df=res$final.fits$df, direction=c("pos","neg"))
  df <- reshape2::melt(nfeats.all, value.name="nfeats")
  df2 <- ddply(df, .(category, region, direction, df), function(x) {
    c(nfeats=mean(x$nfeats))
  })
  data.frame(subject=si, df2)
}, .parallel=TRUE)

# subs.nfeats.all <- ldply(1:20, function(si) {
#   grp <- dat.subs[[si]]$grouping
#   fac <- factor(grp, levels=unique(grp), labels=sroi.names[1:4])
#   res <- lst.res[[si]]$res
#   nfeats.all <- laply(1:4, function(fi) {
#     betas <- res$final.fits$beta[[fi]]
#     nfeats.pos <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x>0)))
#     nfeats.neg <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x<0)))
#     abind::abind(pos=nfeats.pos, neg=nfeats.neg, along = 3)
#   })
#   dimnames(nfeats.all) <- list(category=categories.title, region=sroi.names[1:4], 
#                                df=res$final.fits$df, direction=c("pos","neg"))
#   df <- reshape2::melt(nfeats.all, value.name="nfeats")
#   df2 <- ddply(df, .(category, region, direction, df), function(x) {
#     c(nfeats=mean(x$nfeats))
#   })
#   data.frame(subject=si, df2)
# }, .parallel=TRUE)

#' ### AUC
#' 
#' #### Get Data
#' 
#+ auc
subs.nfeats.auc <- ddply(subs.nfeats.all, .(subject, category, region, direction),  function(x) {
  #x <- dlply(subs.nfeats.all, .(subject, category, region, direction))[[1]]
  c(auc=caTools::trapz(x$df, x$nfeats))
})
ave.nfeats.auc <- ddply(subs.nfeats.auc, .(category, region, direction), function(x) {
  c(auc=mean(x$auc))
})

#' #### Plot
#' 
#+ auc-plot
# test plot
library(ggthemes)

p <- ggplot(subset(ave.nfeats.auc, direction=="pos"), aes(x=region, y=auc)) +       geom_bar(stat="identity", position=position_dodge(), aes(fill=category)) + 
  theme_few()
plot(p)
ggsave(file.path(plotdir, "auc_pos.pdf"), p, width=8, height=5)

p <- ggplot(subset(ave.nfeats.auc, direction=="neg"), aes(x=region, y=auc)) +       geom_bar(stat="identity", position=position_dodge(), aes(fill=category)) + 
  theme_few()
plot(p)
ggsave(file.path(plotdir, "auc_neg.pdf"), p, width=8, height=5)


#' #### Stats
#' 
#' auc-stats
#tmp <- subset(subs.nfeats.auc, direction=="pos" & region=="R FFA")
#tmp2 <- pairwise.t.test(tmp$auc, tmp$category, p.adjust.method = "none")
#tmp2
tmp2 <- dlply(subset(subs.nfeats.auc, direction=="pos"), .(category), function(tmp) {
  pairwise.t.test(tmp$auc, tmp$region, p.adjust.method = "none")
})
tmp2
dlply(subset(subs.nfeats.auc, direction=="neg"), .(region), function(tmp) {
  pairwise.t.test(tmp$auc, tmp$category, p.adjust.method = "none")
})






### TESTS BELOW ###

# testing out per subject
load(outfile)

sapply(dat.subs, function(x) length(x$grouping))
llply(lst.res, function(x) tmpfun1(x$res))


ret <- laply(1:20, function(si) {
  sgrouping <- dat.subs[[si]]$grouping
  sfac      <- factor(sgrouping, levels=1:4, labels=c("FFA", "PPA", "LOC", "VWFA"))
  sres      <- lst.res[[si]]$res
  ret <- laply(1:4, function(fi) {
    cbind(pos=table(sfac[sres$best$betas[fi,]>0]), neg=table(sfac[sres$best$betas[fi,]<0]))
  })
  dimnames(ret)[[1]] <- categories.title
  
  ret
})
round(apply(ret, 2:4, mean), 1)


# this one tries to get how many subjects have any weights in that region.
# note: not much comes out with just averaging the weights
ret0 <- laply(1:20, function(si) {
  sgrouping <- dat.subs[[si]]$grouping
  sfac      <- factor(sgrouping, levels=1:4, labels=c("FFA", "PPA", "LOC", "VWFA"))
  sres      <- lst.res[[si]]$res
  ret <- laply(1:4, function(fi) {
    betas <- sres$best$betas[fi,]
    cbind(pos=tapply(betas[betas>0], sfac[betas>0], mean), 
          neg=tapply(betas[betas<0], sfac[betas<0], mean))
  })
  dimnames(ret)[[1]] <- categories.title
  
  ret
})
round(apply(ret0, 2:4, function(x) sum(!is.na(x))), 1)

ret

table(sfac[sres$best$betas[1,]<0])

names(dat.subs[[2]])

# Let's try this with the subsample selection
lasso.firstq.multi <- function (x, y, q, sel.level, ...) 
{
  fit   <- glmnet(x, y, dfmax = q, family="multinomial", ...)
  m     <- predict(fit, type = "nonzero")
  m     <- m[[sel.level]]
  #m     <- lapply(1:length(m[[1]]), function(i) {
  #  unique(unlist(sapply(m, function(mm) mm[[i]])))
  #})
  delta <- q - unlist(lapply(m, length))
  delta[delta < 0] <- Inf
  take  <- which.min(delta)
  m[[take]]
}

lst.hret <- llply(1:4, function(li) {
  cat(li,"\n")
  llply(1:20, function(si) {
    hret <- hdi(as.matrix(dat.subs[[si]]$dat), dat.subs[[si]]$y, 
                method="stability", B=100, fraction=0.8, 
                model.selector=lasso.firstq.multi, 
                args.model.selector=list(nlambda=100, sel.level=li, alpha=1), 
                EV=20, threshold=1, 
                ncores=24, parallel=TRUE, verbose=F)
    hret
  }, .progress="text")
})
  
ret1 <- laply(1:20, function(si) {
  sgrouping <- dat.subs[[si]]$grouping
  sfac      <- factor(sgrouping, levels=1:4, labels=c("FFA", "PPA", "LOC", "VWFA"))
  sres      <- lst.res[[si]]$res
  freqs     <- sapply(lst.hret, function(x) x[[si]]$freq)
  ret <- laply(1:4, function(fi) {
    betas <- sres$best$betas
    cbind(pos=table(sfac[freqs[,fi]>0.5 & betas[fi,]>0]), neg=table(sfac[freqs[,fi]>0.5 & betas[fi,]<0]))
  })
  ret
})
round(apply(ret1, 2:4, mean))





#' ## Test Standard Space
#'  
#' It could be that we are getting 
#' 
#+ test-std
dat.subs.std <- laply(subjects, function(subject) {
  cat(subject, "\n")
  sublst       <- load.data.std(subject)
  scale(sublst$dat)
})

sublst       <- load.data.std(subjects[2])

#'  Load group data?

# Get the folds
nrepeats <- 1
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- createMultiFolds(uruns, k=nruns, times = nrepeats)
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})

# Get the classification in standard space
lst.sres <- llply(1:length(subjects), function(si) {
  #cat(si, "\n")
  
  # Loop through each alpha
  sdat     <- dat.subs.std[si,,]
  if (si==6) { # one subject has a column with NAs
    w <- apply(sdat, 2, function(x) any(is.na(x)))
    sdat[,w] <- rowMeans(sdat[,!w])
  }
  obs      <- sublst$yfactor
  rois     <- grouping
  
  # Run the model!!!
  res      <- wrap.model(sdat, obs, foldsI, alpha=1, nlambda=100, parallel=T)
  #res      <- save.model(res, rois)
  
  res
}, .parallel=F, .progress="text")

# Get the table with features
sret <- laply(1:20, function(si) {
  sres      <- lst.sres[[si]]
  
  inds     <- grouping %in% c(1,2,4,6)
  grp      <- grouping[inds]
  betas    <- sres$best$betas[,inds]
  
  fac      <- factor(grp, levels=c(1,2,4,6), labels=c("FFA", "PPA", "LOC", "VWFA"))
  
  ret <- laply(1:4, function(fi) {
    cbind(pos=table(fac[betas[fi,]>0]), neg=table(fac[betas[fi,]<0]))
  })
  dimnames(ret)[[1]] <- categories.title
  
  ret
})
round(apply(sret, 2:4, mean), 1)


#' Get the number of feats by varying lambda
#+ feats-by-lambda
inds     <- grouping %in% c(1,2,4,6)
fac <- factor(grouping[inds], levels=unique(grouping[inds]), labels=sroi.names[1:4])
subs.nfeats.all <- ldply(1:20, function(si) {
  res <- lst.sres[[si]]
  nfeats.all <- laply(1:4, function(fi) {
    betas <- res$final.fits$beta[[fi]]
    nfeats.pos <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x>0)))
    nfeats.neg <- apply(betas, 2, function(x) tapply(x, fac, function(x) sum(x<0)))
    abind::abind(pos=nfeats.pos, neg=nfeats.neg, along = 3)
  })
  dimnames(nfeats.all) <- list(category=categories.title, region=sroi.names[1:4], 
                               df=res$final.fits$df, direction=c("pos","neg"))
  df <- reshape2::melt(nfeats.all, value.name="nfeats")
  df2 <- ddply(df, .(category, region, direction, df), function(x) {
    c(nfeats=mean(x$nfeats))
  })
  data.frame(subject=si, df2)
}, .parallel=TRUE)


ave.nfeats.all <- ddply(subs.nfeats.all, .(category, region, direction, df), function(x) {
  c(nfeats=mean(x$nfeats))
})

udfs <- sort(unique(ave.nfeats.all$df))
smooth.nfeats.all <- ddply(ave.nfeats.all, .(category, region, direction), function(x) {
  #x <- subset(df2, category=="Faces" & region=="R FFA" & direction=="pos")
  ssfit <- smooth.spline(x$df, x$nfeats)
  pred <- predict(ssfit, udfs)
  data.frame(df=pred$x, nfeats=pred$y)
})

# use subject data to calculate the AUC for each person

ggplot(subset(smooth.nfeats.all, category=="Faces" & region=="R FFA" & direction=="pos"), aes(x=df, y=nfeats)) + geom_line()

ggplot(subset(ave.nfeats.all, category=="Faces" & region=="R FFA" & direction=="pos"), aes(x=df, y=nfeats)) + geom_line()

ggplot(subset(df2, category=="Faces" & region=="R FFA" & direction=="pos"), aes(x=df, y=nfeats)) + geom_line()

p <- ggplot(smooth.nfeats.all, aes(x=df, y=nfeats, color=category)) + 
  #geom_vline(xintercept=-log10(sel.lambda), linetype="dashed") + 
  geom_line(size=1.5) + 
  #geom_point(data=max.df, size=4) + 
  scale_color_manual(values=cat.cols[-1]) + 
  theme_gdocs() + 
  #scale_x_reverse() + 
  facet_grid(direction~region)
plot(p)
