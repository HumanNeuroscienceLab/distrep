#' # Computational Model for Faces
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

library(RColorBrewer)
library(igraph)
library(plyr)

library(doMC)
registerDoMC(24)

# Load the data
indir <- outdir <- "/data1/ffg05/analysis/comp_model/v1_model"

timing <- NULL; catlist <- NULL; layer1.feats <- NULL; layer2.feats <- NULL
load(file.path(indir, "layer1.rda"))
load(file.path(indir, "layer2.rda")) # 50 components
#load(file.path(indir, "layer2_25.rda")) # not that much different than 50...
#load(file.path(indir, "layer2_100.rda"))

layer2.feats <- t(layer2.feats)
#layer2.feats2 <- t(layer2.feats2)
#layer2.feats3 <- t(layer2.feats3)

# Rearrange the rows to match the timing information
stiming <- ddply(timing, .(category), function(x) x)
stim.inds <- stiming$ind # so the inds in this order should be right
layer1.feats <- layer1.feats[,stim.inds]
layer2.feats <- layer2.feats[,stim.inds] # this is the correct arrangement!

# Do the correlation business
system.time(cmat1 <- cor(layer1.feats))
system.time(cmat2 <- cor(layer2.feats))
#system.time(cmat3 <- cor(layer2.feats2))
#system.time(cmat4 <- cor(layer2.feats3))


#' Let's try for this canonical correlation analysis
#' 
#+ cancor
# first let's see how each feature relates to the categories
lfeats <- t(layer2.feats)
cats <- as.character(stiming$category)
tmp <- apply(lfeats, 2, function(x) tapply(x, cats, mean))
round(tmp, 3)*100
# also can do a type of regression
fit <- lm(lfeats ~ cats)
tmp2 <- summary(aov(fit))
pvals <- sapply(tmp2, function(x) x$`Pr(>F)`[1])
zvals <- qt(pvals, Inf, lower.tail=F)
round(zvals, 2)


cvfit <- cv.glmnet(grpdat, lfeats[,1:10], family="mgaussian", 
                   nlambda=100, nfolds=10, parallel=TRUE)
ii <- which(cvfit$lambda == cvfit$lambda.min)
gfit <- cvfit$glmnet.fit
gfit
dim(gfit$beta[[1]])
length(gfit$beta[[1]][,ii])

betas <- sapply(1:10, function(i) gfit$beta[[i]][,ii])
betas[betas==0] <- NA
heatmap(t(betas), Colv = NA, Rowv = NA, labCol = F)

preds <- predict(cvfit, newx=grpdat, s="lambda.min")
preds <- preds[,,1]
tmp3 <- apply(preds, 2, function(x) tapply(x, cats, mean))
round(tmp3, 3)*100
round(tmp[,1:10], 3)*100


cvfit <- cv.glmnet(grpdat, lfeats[,1:1], family="gaussian", 
                   nlambda=100, nfolds=10, parallel=TRUE)
ii <- which(cvfit$lambda == cvfit$lambda.min)
betas <- cvfit$glmnet.fit$beta[,ii]
table(roi.df$names[grouping[which(betas>0)]])
table(roi.df$names[grouping[which(betas<0)]])

cvfit <- cv.glmnet(grpdat, lfeats[,2], family="gaussian", 
                   nlambda=100, nfolds=10, parallel=TRUE)
ii <- which(cvfit$lambda == cvfit$lambda.min)
betas <- cvfit$glmnet.fit$beta[,ii]
table(roi.df$names[grouping[which(betas>0)]])
table(roi.df$names[grouping[which(betas<0)]])

cvfit <- cv.glmnet(grpdat, lfeats[,3], family="gaussian", 
                   nlambda=100, nfolds=10, parallel=TRUE)
ii <- which(cvfit$lambda == cvfit$lambda.min)
betas <- cvfit$glmnet.fit$beta[,ii]
table(roi.df$names[grouping[which(betas>0)]])
table(roi.df$names[grouping[which(betas<0)]])

cvfit <- cv.glmnet(grpdat, lfeats[,20], family="gaussian", 
                   nlambda=100, nfolds=10, parallel=TRUE)
ii <- which(cvfit$lambda == cvfit$lambda.min)
betas <- cvfit$glmnet.fit$beta[,ii]
table(roi.df$names[grouping[which(betas>0)]])
table(roi.df$names[grouping[which(betas<0)]])


require(CCA)
cc1 <- cc(grpdat, lfeats)

library(PMA)
perm.out <- CCA.permute(grpdat, lfeats, typex="standard", typez="standard", nperms=25)
out <- CCA(grpdat, lfeats, typex="standard", typez="standard", K=10, 
           penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz)
betas <- out$u[,8]
table(roi.df$names[grouping[which(betas>0)]])
table(roi.df$names[grouping[which(betas<0)]])
zz <- t(t(out$v) %*% t(lfeats))
tmp3 <- apply(zz, 2, function(x) tapply(x, cats, mean))
round(tmp3, 3)*100


library(nscancor)
ret2 <- nscancor(grpdat, lfeats)

# ok so what i want to know is when i predict each of those features
# what is it that i get (what brain areas)

# Plotting please
#' ### Graph Network
#' 
#' We plot the correlation matrix as a thresholded network
#' 
#+ network
#detach("package:network", unload=TRUE)
delete.isolates <- function(graph, mode = 'all') {
  isolates <- which(degree(graph, mode = mode) == 0)
  delete.vertices(graph, isolates)
}
g     <- graph.adjacency(cmat1*(cmat1>0.3), mode="undirected", weighted=T, diag=F)
g     <- delete.isolates(g)
ucols <- brewer.pal(8, "Set2")[c(3,5,2,1)] # Faces, Fruits, Letters, & Vehicles
ucats <- levels(factor(catlist))
#roinds<- unlist(lapply(ucats[uinds], function(ucat) which(catlist==ucat)))
cols  <- ucols[as.numeric(factor(catlist))]
plot.igraph(g, 
            edge.width=E(g)$weight/2, 
            vertex.size=3,
            vertex.color=cols, 
            vertex.label=NA, 
            vertex.frame.color = cols)



# with pca output
delete.isolates <- function(graph, mode = 'all') {
  isolates <- which(degree(graph, mode = mode) == 0)
  delete.vertices(graph, isolates)
}
## try to vary the threshold a bit
g     <- graph.adjacency(cmat2*(cmat2>0.1), mode="undirected", weighted=T, diag=F)
g     <- delete.isolates(g)
ucols <- brewer.pal(8, "Set2")[c(3,5,2,1)] # Faces, Fruits, Letters, & Vehicles
ucats <- levels(factor(catlist))
#roinds<- unlist(lapply(ucats[uinds], function(ucat) which(catlist==ucat)))
cols  <- ucols[as.numeric(factor(catlist))]
plot.igraph(g, 
            edge.width=E(g)$weight/2, 
            vertex.size=3,
            vertex.color=cols, 
            vertex.label=NA, 
            vertex.frame.color = cols)

# Now I want to know what to know how well the classifier can separete those
# categories. I should do a form of cross-validation
# So I could just use the weights 
# 
# Ok so let's try the multinomial logistic regression
# Actually I just need the catlist here
library(nnet)
cv.multinom <- function(x, y, runs, nrepeats=1, ...) {
  uruns    <- sort(unique(runs))
  nruns    <- length(uruns)
  runFolds <- caret::createMultiFolds(uruns, k=nruns, times = nrepeats)
  foldsI   <- lapply(runFolds, function(x) {
    runs %in% x
  })
  
  uy <- levels(y)
  n <- length(y)
  pred <- vector("character", n)
  prob <- matrix(0, n, length(uy))
  
  for (i in 1:length(foldsI)) {
    cat("Fold/rep", i, "\n")
    
    train.inds <- foldsI[[i]]
    train.x    <- as.data.frame(x[train.inds,])
    train.y    <- y[train.inds]
    
    test.inds  <- !foldsI[[i]]
    test.x     <- as.data.frame(x[test.inds,])
    #test.y     <- y[test.inds]
    
    trainModel <- multinom(train.y ~ ., data=train.x, ...)
    
    test.pred  <- predict(trainModel, newdata=test.x, type='class')
    test.prob  <- predict(trainModel, newdata=test.x, type='prob')
    
    pred[test.inds]  <- as.character(test.pred)
    prob[test.inds,] <- test.prob
    
    cat("\n")
  }
  
  pred <- factor(pred)
  dimnames(prob) <- list(trial=y, category=uy)
  
  list(pred=pred, prob=prob)
}

x    <- as.matrix(t(layer2.feats))
colnames(x) <- sprintf("comp%02i", 1:50)
y    <- factor(catlist)
runs <- stiming$run
mres <- cv.multinom(x, y, runs)
caret::confusionMatrix(mres$pred, y)

plot.ts(mres$prob)




library(glmnet)
wrap.model <- NULL
source("/data1/ffg05/scripts/repsim/scratch/stats_lab_funs.R")

nrepeats <- 1
runs     <- stiming$run
uruns    <- sort(unique(runs))
nruns    <- length(uruns)
runFolds <- caret::createMultiFolds(uruns, k=nruns, times = nrepeats)
foldsI   <- lapply(runFolds, function(x) {
  runs %in% x
})

gres <- wrap.model(x, y, foldsI, alpha=1, nlambda=100, parallel=T)

## get probs for best lambda
lind  <- gres$best$stats$lind
probs <- gres$probs[,,lind]
colnames(probs) <- unique(catlist)
plot.ts(probs, ylim=c(0,1))

caret::confusionMatrix(gres$preds[[lind]], gres$obs)

cvfit <- cv.glmnet(layer1.feats, y, nlambda=50)
?predict.nnet


# Okay but if I want to look at the trial by trial relationship than note
# that my timing file order is the order of my brain patterns
# 
# However, the order of this stimulus data is given the ind column
# I could just pass that ind column to get the rearranged stim data
# Oh wait but the betas are actually ordered by category
stim.inds <- ddply(timing, .(category), function(x) x)$ind # so the inds in this order should be right
layer2.feats[stim.inds,] # this is the correct arrangement!
## we could do a canonical correlation
## the idea would be to see which components differentiate the categories
## and where they might be found in the brain
## ...wonder if the sparse canonical one might be better

## like we said before, we can also run some analysis and get the probabilities
## of each of the category


nnet.res <- nnet(y ~ ., data=as.data.frame(x), size=2, trace=T)


library(neuralnet)
data <- as.data.frame(cbind(x, class.ind(catlist)))
f <- paste("Faces", "~", paste(colnames(x), collapse=" + "))
net <- neuralnet(f, data=data, hidden=c(16,8))
plot(net)
?neuralnet
