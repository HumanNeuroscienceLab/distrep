

#' This script will test out the hdi on the group average data using
#' stability selection
#' 

library(hdi)

hdi <- function (x, y, method = "multi.split", B = NULL, fraction = 0.5, 
          model.selector = NULL, EV = NULL, threshold = 0.75, gamma = seq(0.05, 
                                                                          0.99, by = 0.01), classical.fit = NULL, args.model.selector = NULL, 
          args.classical.fit = NULL, verbose = FALSE, ...) 
{
  if (method %in% c("multi.split", "multi-split")) {
    if (is.null(B)) 
      B <- 100
    if (is.null(model.selector)) 
      model.selector <- lasso.cv
    if (is.null(classical.fit)) 
      classical.fit <- lm.pval
    out <- multi.split(x = x, y = y, B = B, fraction = fraction, 
                       model.selector = model.selector, classical.fit = classical.fit, 
                       gamma = gamma, args.model.selector = args.model.selector, 
                       args.classical.fit = args.classical.fit, verbose = verbose, 
                       ...)
  }
  else if (method == "stability") {
    if (is.null(B)) 
      B <- 100
    if (is.null(EV)) 
      stop("EV not defined")
    if (is.null(model.selector)) 
      model.selector <- lasso.firstq
    out <- stability(x = x, y = y, B = B, fraction = fraction, 
                     threshold = threshold, model.selector = model.selector, 
                     EV = EV, args.model.selector = args.model.selector, 
                     verbose = verbose, ...)
  }
  else {
    stop("Method not (yet) defined")
  }
  out$call <- match.call()
  class(out) <- "hdi"
  out
}


stability <- function (x, y, EV, threshold = 0.75, B = 100, fraction = 0.5, 
                       model.selector = lasso.firstq, args.model.selector = NULL, 
                       parallel = FALSE, ncores = getOption("mc.cores", 2L), verbose = FALSE) 
{
  if (threshold > 1 | threshold < 0.5) 
    stop("threshold has to be in (0.5, 1)")
  n <- nrow(x)
  p <- ncol(x)
  col.nam <- colnames(x)
  q <- ceiling(sqrt(EV * p * (2 * threshold - 1)))
  sel.mat <- matrix(FALSE, nrow = B, ncol = p)
  sel.n <- floor(fraction * n)
  oneSample <- function(...) {
    sel <- sample(1:n, sel.n, replace = FALSE)
    x.sel <- x[sel, ]
    y.sel <- y[sel]
    sel.model <- do.call(model.selector, c(list(x = x.sel, 
                                                y = y.sel, q = q), args.model.selector))
    out <- logical(ncol(x))
    out[sel.model] <- TRUE
    out
  }
  if (parallel) {
    if (verbose) 
      cat("...starting parallelization of bootstrap samples\n")
    sel.mat <- matrix(unlist(mclapply(1:B, oneSample, mc.cores = ncores)), 
                      nrow = B, byrow = TRUE)
  }
  else {
    for (b in 1:B) {
      if (verbose) 
        cat("...Subsample", b, "\n")
      sel.mat[b, ] <- oneSample()
    }
  }
  freq <- colMeans(sel.mat)
  names(freq) <- col.nam
  out <- list()
  sel.current <- which(freq >= threshold)
  names(sel.current) <- col.nam[sel.current]
  if (length(sel.current) == 0) 
    sel.current <- NULL
  out <- sel.current
  out <- list(selected = sel.current, EV = EV, threshold = threshold, 
              freq = freq, q = q, method = "stability", call = match.call())
  class(out) <- "hdi"
  return(out)
}

# Easiest way will be to select one of the 
lasso.firstq.multi <- function (x, y, q, sel.level, ...) 
{
  fit   <- glmnet(x, y, dfmax = q, family="multinomial", ...)
  m     <- predict(fit, type = "nonzero")
  #m     <- m[[sel.level]]
  m     <- lapply(1:length(m[[1]]), function(i) {
    unique(unlist(sapply(m, function(mm) mm[[i]])))
  })
  delta <- q - unlist(lapply(m, length))
  delta[delta < 0] <- Inf
  take  <- which.min(delta)
  m[[take]]
}



#prefit <- glmnet(x, as.numeric(obs), family=family, alpha=alpha, 
#                 nlambda=nlambda)

ret <- hdi(as.matrix(subset.grpdat), as.numeric(ys), 
            method="stability", B=100, fraction=0.75, 
            model.selector=lasso.firstq.multi, 
            args.model.selector=list(nlambda=100, sel.level=1, alpha=1), 
            EV=10, threshold=1, 
            ncores=24, parallel=TRUE, verbose=TRUE)




# Ok so this is actually to get the p-values
lasso.cv2 <- function (x, y, nfolds = 10, grouped = nrow(x) > 3 * nfolds, sel.level=1, 
                       ...) 
{
  fit.cv <- cv.glmnet(x, y, nfolds = nfolds, grouped = grouped, 
                      ...)
  sel <- predict(fit.cv, type = "nonzero")
  sel[[sel.level]]
}

lm.pval2 <- function (x, y, exact = TRUE, ...) {
  fit.lm <- lm(y ~ x, ...)
  fit.summary <- summary(fit.lm)
  tstat <- coef(fit.summary)[-1, "t value"]
  setNames(2 * (if (exact) 
    pt(abs(tstat), df = fit.lm$df.residual, lower.tail = FALSE)
    else pnorm(abs(tstat), lower.tail = FALSE)), colnames(x))
}

ret <- hdi(as.matrix(subset.grpdat), as.numeric(ys), 
           method="multi.split", B=100, fraction=0.75, 
           model.selector=lasso.cv2, 
           args.model.selector=list(nlambda=100, sel.level=1, alpha=1), 
           EV=10, threshold=1, 
           ncores=24, parallel=TRUE, verbose=TRUE)

