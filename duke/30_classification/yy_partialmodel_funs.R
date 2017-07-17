load.data <- function(subject) {
  base       <- "/data1/ffg05"
  roidir     <- file.path(base, "analysis/rois")
  taskdir    <- file.path(base, "analysis/task_activity")
  categories <- c("Faces", "Letters", "Fruits", "Vehicles")
  
  # Reads in the neurosynth ROIs
  # There will be a total of 12 ROIs
  roifn    <- file.path(roidir, subject, "neurosynth_visual_funcspace.nii.gz")
  rois     <- read.mask(roifn, NULL)
  mask     <- rois>10 # exclude L FFA
  grouping <- rois[mask]
  grouping <- as.numeric(as.factor(grouping)) # redo to be 1:n
  roi.names<- c("R FFA", "L PPA", "R PPA", "L LOC", "R LOC", "L VWF", 
                "L V1", "R V1", "L M1", "R M1", "L A1", "R A1")
  
  # Read in the beta-series data
  dat <- ldply(categories, function(category) {
    datfn  <- file.path(taskdir, subject, "Categories", "beta_series_spmg1.reml", 
                        sprintf("beta_series_%s.nii.gz", tolower(category)))
    dat    <- read.big.nifti(datfn)
    dat    <- deepcopy(dat, mask)
    mdat   <- as.matrix(dat)
    rm(dat); gc()
    mdat
  }, .progress="text")
  dat <- as.matrix(dat) # observations x
  
  # Reorder the columns to have all the ROIs in order
  os       <- order(grouping)
  grouping <- grouping[os]
  dat      <- dat[,os]
  
  # Get the associated categories
  # i.e. the labels associated with each observation in `dat`
  ntrials.per.category <- laply(categories, function(category) {
    datfn  <- file.path(taskdir, subject, "Categories", "beta_series_spmg1.reml", 
                        sprintf("beta_series_%s.nii.gz", tolower(category)))
    hdr    <- read.nifti.header(datfn)
    hdr$dim[4]
  })
  classes <- factor(rep(categories, ntrials.per.category), levels=categories)
  y       <- as.numeric(classes)
  
  # For a leave-one run out cross-validation
  # Get the runs associated with each observation
  timing <- read.csv("/data1/ffg05/notes/timing.csv")
  tmp    <- ddply(timing[,-1], .(category), function(x) data.frame(run=x$run))
  tmp    <- tmp[tmp$category!="Circles",] # no circles
  tmp    <- tmp[tmp$category!="Objects",] # no objects
  # reorder
  inds   <- unlist(lapply(categories, function(x) which(tmp$category==x)))
  tmp$run      <- tmp$run[inds]
  tmp$category <- factor(tmp$category[inds])
  if(any(tmp$category != classes)) {
    print(which(tmp$category != classes))
    stop("timing file categories != classes")
  }
  runs   <- tmp$run
  
  list(dat=dat, grouping=grouping, yfactor=classes, y=y, runs=runs, 
       rois=sort(unique(grouping)), categories=categories, 
       roi.names=roi.names)
}


# Loads the functional data in standard space
load.data.std <- function(subject) {
  base       <- "/data1/ffg05"
  roidir     <- file.path(base, "analysis/rois")
  taskdir    <- file.path(base, "analysis/task_activity")
  categories <- c("Faces", "Letters", "Fruits", "Vehicles")
  
  # Reads in the neurosynth ROIs
  # There will be a total of 12 ROIs
  roifn    <- file.path(roidir, "group", "neurosynth_visual_rois_3mm.nii.gz")
  rois     <- read.mask(roifn, NULL)
  mask     <- rois>10 # exclude L FFA
  grouping <- rois[mask]
  grouping <- as.numeric(as.factor(grouping)) # redo to be 1:n
  roi.names<- c("R FFA", "L PPA", "R PPA", "L LOC", "R LOC", "L VWF", 
                "L V1", "R V1", "L M1", "R M1", "L A1", "R A1")
  
  # Read in the beta-series data
  dat <- ldply(categories, function(category) {
    datfn  <- file.path(taskdir, subject, "Categories", "beta_series_spmg1.reml", 
                        "reg_standard", 
                        sprintf("beta_series_%s.nii.gz", tolower(category)))
    dat    <- read.big.nifti(datfn)
    dat    <- deepcopy(dat, mask)
    mdat   <- as.matrix(dat)
    rm(dat); gc()
    mdat
  }, .progress="text")
  dat <- as.matrix(dat) # observations x
  
  # Reorder the columns to have 12 ROIs in order
  os       <- order(grouping)
  grouping <- grouping[os]
  dat      <- dat[,os]
  
  # Get the associated categories
  # i.e. the labels associated with each observation in `dat`
  ntrials.per.category <- laply(categories, function(category) {
    datfn  <- file.path(taskdir, subject, "Categories", "beta_series_spmg1.reml", 
                        "reg_standard", 
                        sprintf("beta_series_%s.nii.gz", tolower(category)))
    hdr    <- read.nifti.header(datfn)
    hdr$dim[4]
  })
  classes <- factor(rep(categories, ntrials.per.category), levels=categories)
  y       <- as.numeric(classes)
  
  # For a leave-one run out cross-validation
  # Get the runs associated with each observation
  timing <- read.csv("/data1/ffg05/notes/timing.csv")
  tmp    <- ddply(timing[,-1], .(category), function(x) data.frame(run=x$run))
  tmp    <- tmp[tmp$category!="Circles",] # no circles
  tmp    <- tmp[tmp$category!="Objects",] # no objects
  # reorder
  inds   <- unlist(lapply(categories, function(x) which(tmp$category==x)))
  tmp$run      <- tmp$run[inds]
  tmp$category <- factor(tmp$category[inds])
  if(any(tmp$category != classes)) {
    print(which(tmp$category != classes))
    stop("timing file categories != classes")
  }
  runs   <- tmp$run
  
  list(dat=dat, grouping=grouping, yfactor=classes, y=y, runs=runs, 
       rois=sort(unique(grouping)), categories=categories, 
       roi.names=roi.names)
}

get_stats <- function(prob, pred, obs) {
  require(caret)
  suppressMessages(require(Metrics))
  
  if (!is.factor(pred)) stop("pred must be a factor")
  if (!is.factor(obs)) stop("obs must be a factor")
  
  colnames(prob) <- levels(obs)
  data <- data.frame(
    pred=pred, 
    obs=obs, 
    prob
  )
  
  prob_stats <- lapply(levels(data[, "pred"]), function(class) {
    #Grab one-vs-all data for the class
    pred <- ifelse(data[, "pred"] == class, 1, 0)
    obs  <- ifelse(data[,  "obs"] == class, 1, 0)
    prob <- data[,class]
    
    #Calculate one-vs-all AUC and logLoss and return
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('ROC', 'logLoss')
    
    return(prob_stats) 
  })
  prob_stats
  
  prob_stats <- do.call(rbind, prob_stats)
  rownames(prob_stats) <- paste('Class:', levels(data[, "pred"]))
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"])
  #Aggregate and average class-wise stats
  #Todo: add weights
  class_stats <- cbind(CM$byClass, prob_stats)
  class_stats <- colMeans(class_stats)
  #Aggregate overall stats
  overall_stats <- c(CM$overall)
  #Combine overall with class-wise stats and remove some stats we don't want 
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', 
                                       'Prevalence', 'Detection Prevalence')]
  #Clean names and return
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  
  stats
}

# some redundancy with above
category_stats <- function(obs, pred, prob) {
  suppressMessages(require(Metrics))
  
  ldply(levels(pred), function(class) {
    #Grab one-vs-all data for the class
    pred <- ifelse(pred == class, 1, 0)
    obs  <- ifelse(obs == class, 1, 0)
    prob0 <- prob[,class]
    
    #Calculate accuracy, one-vs-all AUC, and logLoss
    accuracy <- function(obs, pred) sum(pred & obs)/sum(obs)
    cap_prob <- pmin(pmax(prob, .000001), .999999)
    prob_stats <- c(accuracy(obs, pred), auc(obs, prob), logLoss(obs, cap_prob))
    names(prob_stats) <- c('Accuracy', 'ROC', 'logLoss')
    
    cbind(class=class, t(prob_stats))
  })
}

# taken from package heatmap.plus
heatmap.plus <- function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
                          distfun = dist, hclustfun = hclust, reorderfun = function(d, 
                                                                                    w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, 
                                                                                                                                               "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
                          margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 
                            1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
                          labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
                          verbose = getOption("verbose"), ...) 
{
  scale <- if (symm && missing(scale)) 
    "none"
  else match.arg(scale)
  if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
    stop("'x' must be a numeric matrix")
  nr <- di[1]
  nc <- di[2]
  if (nr <= 1 || nc <= 1) 
    stop("'x' must have at least 2 rows and 2 columns")
  if (!is.numeric(margins) || length(margins) != 2) 
    stop("'margins' must be a numeric vector of length 2")
  doRdend <- !identical(Rowv, NA)
  doCdend <- !identical(Colv, NA)
  if (is.null(Rowv)) 
    Rowv <- rowMeans(x, na.rm = na.rm)
  if (is.null(Colv)) 
    Colv <- colMeans(x, na.rm = na.rm)
  if (doRdend) {
    if (inherits(Rowv, "dendrogram")) 
      ddr <- Rowv
    else {
      hcr <- hclustfun(distfun(x))
      ddr <- as.dendrogram(hcr)
      if (!is.logical(Rowv) || Rowv) 
        ddr <- reorderfun(ddr, Rowv)
    }
    if (nr != length(rowInd <- order.dendrogram(ddr))) 
      stop("row dendrogram ordering gave index of wrong length")
  }
  else rowInd <- 1:nr
  if (doCdend) {
    if (inherits(Colv, "dendrogram")) 
      ddc <- Colv
    else if (identical(Colv, "Rowv")) {
      if (nr != nc) 
        stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
      ddc <- ddr
    }
    else {
      hcc <- hclustfun(distfun(if (symm) 
        x
        else t(x)))
      ddc <- as.dendrogram(hcc)
      if (!is.logical(Colv) || Colv) 
        ddc <- reorderfun(ddc, Colv)
    }
    if (nc != length(colInd <- order.dendrogram(ddc))) 
      stop("column dendrogram ordering gave index of wrong length")
  }
  else colInd <- 1:nc
  x <- x[rowInd, colInd]
  labRow <- if (is.null(labRow)) 
    if (is.null(rownames(x))) 
      (1:nr)[rowInd]
  else rownames(x)
  else labRow[rowInd]
  labCol <- if (is.null(labCol)) 
    if (is.null(colnames(x))) 
      (1:nc)[colInd]
  else colnames(x)
  else labCol[colInd]
  if (scale == "row") {
    x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
    sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if (scale == "column") {
    x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
    sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }
  lmat <- rbind(c(NA, 3), 2:1)
  lwid <- c(if (doRdend) 1 else 0.05, 4)
  lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
            4)
  if (!missing(ColSideColors)) {
    if (!is.matrix(ColSideColors)) 
      stop("'ColSideColors' must be a matrix")
    if (!is.character(ColSideColors) || dim(ColSideColors)[1] != 
        nc) 
      stop("'ColSideColors' dim()[2] must be of length ncol(x)")
    lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
    lhei <- c(lhei[1], 0.2, lhei[2])
  }
  if (!missing(RowSideColors)) {
    if (!is.matrix(RowSideColors)) 
      stop("'RowSideColors' must be a matrix")
    if (!is.character(RowSideColors) || dim(RowSideColors)[1] != 
        nr) 
      stop("'RowSideColors' must be a character vector of length nrow(x)")
    lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
                                         1), 1), lmat[, 2] + 1)
    lwid <- c(lwid[1], 0.2, lwid[2])
  }
  lmat[is.na(lmat)] <- 0
  if (verbose) {
    cat("layout: widths = ", lwid, ", heights = ", lhei, 
        "; lmat=\n")
    print(lmat)
  }
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
  if (!missing(RowSideColors)) {
    par(mar = c(margins[1], 0, 0, 0.5))
    rsc = RowSideColors[rowInd, , drop=F]
    rsc.colors = matrix()
    rsc.names = names(table(rsc))
    rsc.i = 1
    for (rsc.name in rsc.names) {
      rsc.colors[rsc.i] = rsc.name
      rsc[rsc == rsc.name] = rsc.i
      rsc.i = rsc.i + 1
    }
    rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
    image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
    if (length(colnames(RowSideColors)) > 1) {
      axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), colnames(RowSideColors), 
           las = 2, tick = FALSE)
    }
  }
  if (!missing(ColSideColors)) {
    par(mar = c(0.5, 0, 0, margins[2]))
    csc = ColSideColors[colInd, , drop=F]
    csc.colors = matrix()
    csc.names = names(table(csc))
    csc.i = 1
    for (csc.name in csc.names) {
      csc.colors[csc.i] = csc.name
      csc[csc == csc.name] = csc.i
      csc.i = csc.i + 1
    }
    csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
    image(csc, col = as.vector(csc.colors), axes = FALSE)
    if (length(colnames(ColSideColors)) > 1) {
      axis(2, 0:(dim(csc)[2] - 1)/(dim(csc)[2] - 1), colnames(ColSideColors), 
           las = 2, tick = FALSE)
    }
  }
  par(mar = c(margins[1], 0, 0, margins[2]))
  if (!symm || scale != "none") {
    x <- t(x)
  }
  if (revC) {
    iy <- nr:1
    ddr <- rev(ddr)
    x <- x[, iy]
  }
  else iy <- 1:nr
  image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
          c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
  axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexCol)
  if (!is.null(xlab)) 
    mtext(xlab, side = 1, line = margins[1] - 1.25)
  axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
       cex.axis = cexRow)
  if (!is.null(ylab)) 
    mtext(ylab, side = 4, line = margins[2] - 1.25)
  if (!missing(add.expr)) 
    eval(substitute(add.expr))
  par(mar = c(margins[1], 0, 0, 0))
  if (doRdend) 
    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
  else frame()
  par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
  if (doCdend) 
    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
  else if (!is.null(main)) 
    frame()
  if (!is.null(main)) 
    title(main, cex.main = 1.5 * op[["cex.main"]])
  invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                              doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
}


# y = data; X = model matrix
qlm <- function(y, X) {
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  n <- nrow(y)
  mdf <- ncol(X)
  rdf <- n - mdf
  
  # beta values
  iXtX <- solve(t(X) %*% X)
  b <- iXtX %*% t(X) %*% y
  
  # errors
  e <- y - X %*% b
  
  # sum-square of errors
  SSE <- colSums(e^2)
  
  # standard error
  MSE <- SSE/rdf
  diXtX <- diag(iXtX)
  stderr <- sqrt(sapply(MSE, function(xMSE)
    xMSE * diXtX
  ))
  
  # t-statistic
  tvals <- b/stderr
  
  return(list(
    dd=diXtX, 
    mse=MSE, 
    b=b,
    stderr=stderr,
    tvals=tvals,
    pvals=2*pt(abs(tvals), rdf, lower.tail=F), 
    rdf=rdf
  ))
}

qlm.contrasts <- function(qlm.res, contrasts) {
  contrast_coefs  <- contrasts %*% qlm.res$b
  
  contrast_dd     <- contrasts %*% qlm.res$dd
  contrast_se     <- sapply(qlm.res$mse, function(mse) {
    sqrt(mse %*% contrast_dd)
  })
  
  contrast_tvals  <- contrast_coefs/contrast_se
  
  return(list(
    dd=contrast_dd, 
    b=contrast_coefs, 
    se=contrast_se,
    tvals=contrast_tvals
  ))
}


###
# CLASSIFICATION FUNCTION
###
wrap.model <- function(x, obs, foldsI, alpha=1, nlambda=100, parallel=T) {
  # 1. Get a range of Lambdas
  prefit <- glmnet(x, as.numeric(obs), family="multinomial", alpha=alpha, 
                   nlambda=nlambda)
  lambdas <- prefit$lambda
  
  # 2. Fit the model using a leave-one-run out procedure
  fits <- llply(1:length(foldsI), function(k) {
    indsT   <- foldsI[[k]]
    classesT<- obs[indsT]
    xT      <- x[indsT,]
    fit     <- glmnet(xT, as.numeric(classesT), family="multinomial", 
                      alpha=alpha, lambda=lambdas)
    fit
  }, .parallel=parallel)
  
  # 3. Get the probabilities and predictions for the outcome (categories)
  probs <- array(NA, c(nrow(x), length(levels(obs)), length(lambdas)))
  preds <- matrix(NA, nrow(x), length(lambdas))
  for (k in 1:length(foldsI)) {
    indsT   <- foldsI[[k]]
    # note: if the fit didn't converge than some lambda fits weren't directly
    # calculated. and so glmnet will use linear interpolations to make the 
    # predictions
    probs0 <- predict(fits[[k]], newx=x[!indsT,], type="response", s=lambdas)
    preds0 <- predict(fits[[k]], newx=x[!indsT,], type="class", s=lambdas)
    #       # to check for which lambda's converged
    #       lambdas0 <- fits[[k]]$lambda
    #       linds    <- lambdas %in% lambdas0 # lambda inds
    oinds          <- which(!foldsI[[k]]) # obs inds
    probs[oinds,,] <- probs0
    preds[oinds,]  <- preds0
    rm(probs0, preds0)
  }
  # relabel the prediction outputs
  preds <- alply(preds, 2, function(xx) {
    factor(as.numeric(xx), levels=unique(as.numeric(obs)), 
           labels=levels(obs))
  })
  
  # 4. Model statistics like accuracy to help determine best lambda etc
  #    (based on caret)...compiles everything into one data frame
  all.stats <- ldply(1:length(preds), function(i) {
    oinds <- !is.na(preds[[i]]) # remove any observations without predictions
    row   <- get_stats(probs[oinds,,i], preds[[i]][oinds], obs[oinds])
    data.frame(alpha=alpha, lambda=lambdas[i], t(row))
  }, .parallel=parallel)
  
  # 5. Refit the model using the full data
  #    this way we can get the betas and feature information
  final.fits <- glmnet(x, as.numeric(obs), family="multinomial", alpha=alpha, 
                       lambda=lambdas)
  # Number of non-zero features in any class for each lambda
  betas <- coef(final.fits, s=lambdas)
  betas <- laply(betas, as.matrix) # make into array
  nfeats <- aaply(betas, .(3), function(x) sum(colSums(x[,-1]!=0)!=0))
  # add this to the statistics
  all.stats <- cbind(all.stats, nfeats)
  all.stats$lind <- 1:nrow(all.stats)
  #     # If you wanted to do the same thing as above with each fold fit
  #     tmp <- sapply(1:10, function(k) {
  #       betas <- laply(1:5, function(i) as.matrix(fits[[k]]$beta[[i]]))
  #       feats <- apply(betas, c(2,3), function(xx) sum(xx!=0)>0)
  #       nfeats <- colSums(feats)
  #       nfeats
  #     })
  #     rowMeans(tmp)
  
  # 6. Gather information for the best model fit
  # this gives me the average accuracy, logLoss, ROC, and nfeats
  lind <- which.max(all.stats$Accuracy)
  best.stat <- all.stats[lind,]
  # get the betas
  betas   <- coef(final.fits, s=lambdas[lind])
  betas   <- laply(betas, as.matrix) # make into array
  obetas  <- betas      # keep intercept
  betas   <- betas[,-1] # remove the intercept
  
  # 7. Gather information for the best model fit per category
  # get the accuracy, logloss, and ROC
  pred <- preds[[lind]]
  prob <- probs[,,lind]
  colnames(prob) <- levels(pred)
  cat_stats <- category_stats(obs, pred, prob)
  # get the # of non-zero features for each category
  nfeats <- aaply(betas, .(1), function(x) sum(x!=0))
  cat_stats <- cbind(cat_stats, nfeats)
  best.category.stat <- data.frame(alpha=alpha, lambda=lambdas[lind], cat_stats)
  ## We also want to save the overlap between categories
  feats   <- t((betas!=0)*1)
  overlap <- crossprod(feats)
  dimnames(overlap) <- list(levels(obs), levels(obs))
  perc.overlap <- overlap/nrow(feats)
  
  # Return
  list(alpha=alpha, lambdas=lambdas, fold.fits=fits, final.fits=final.fits, 
       stats=all.stats, obs=obs, preds=preds, probs=probs, 
       best=list(stats=best.stat, category.stats=best.category.stat, 
                 betas=betas, intercept=obetas[,1], 
                 overlap=overlap, perc.overlap=perc.overlap))
}

save.vars <- function(model.lst, rois) {
  list(res=model.lst, 
       acc=model.lst$best$stats$Accuracy, 
       pval=model.lst$best$stats$AccuracyPValue, 
       urois=sort(unique(rois)))
}

wrap.rcfe <- function(sdat, obs, foldsI, rois, alpha=1, nlambda=100) {
  if (is.factor(rois)) {
    urois <- levels(rois)
  } else {
    urois <- sort(unique(rois))
  }
  
  # Do everything
  cat("Full Model\n")
  x          <- sdat
  res        <- wrap.model(x, obs, foldsI, alpha, nlambda)
  full.model <- save.vars(res, rois) # list(res, acc, pval, urois)
  
  # Use accuracy of the full model as a reference
  # and save the unique ROI vals
  current.model <- full.model
  prev.acc      <- full.model$acc
  
  cat("Partial Model\n")
  for (iter in 1:(length(urois)-1)) {
    cat("+++\nIteration:", iter, "\n")
    cat("fitting\n")
    
    # Run each possible model with one of the regions removed
    roi.combns <- llply(1:length(current.model$urois), function(i) {
      current.model$urois[-i]
    })
    reduced.models <- llply(roi.combns, function(reduced.urois) {
      x <- sdat[, rois %in% reduced.urois] # take subset of data
      res <- wrap.model(x, obs, foldsI, alpha, nlambda)
      res <- save.vars(res, rois[rois %in% reduced.urois])
      res
    }, .progress="text")
    
    # Determine the model with the largest positive or least negative change 
    # in accuracy relative to the best model accuracy in the previous iteration
    cat("determining best\n")
    reduced.accs <- laply(reduced.models, function(mod) mod$acc)
    best.ind     <- which.max(reduced.accs - prev.acc)
    best.diff    <- max(reduced.accs - prev.acc)
    
    # If there's no improvement (<0) in accuracy with the best model, then end
    if (best.diff < 0) {
      cat("nothing to remove, stopping!!!\n")
      break
    }
    
    # Save the best model with the region that was removed
    current.model <- reduced.models[[best.ind]]
    prev.acc      <- current.model$acc
  }
  
  cat(sprintf("Best reduced model has %i regions\n", 
              length(current.model$urois)))
  
  # For the purpose of comparison
  cat("Individual ROI Models\n")
  individual.models <- llply(urois, function(uroi) {
    x   <- sdat[,rois==uroi,drop=F]
    res <- wrap.model(x, obs, foldsI, alpha, nlambda)
    res <- save.vars(res, rois[rois==uroi])
    res
  }, .progress="text")
  
  list(
    full = full.model, 
    reduced=current.model, 
    individual=individual.models
  )
}
