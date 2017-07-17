#' # Computational Model for Faces
if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))

#' ## Load Data
#' 
#' Let's load in the data. Note that the `stims` variable is a list with 
#' `stims` being an array with the photos and `catlist` having the corresponding 
#' list of category labels for each trial stimulus.
#' 
#+ load-data
library(R.matlab)
stims0  <- readMat('/data1/ffg05/data/stimuli/stim.mat')
catlist <- unlist(stims0$catlist, use.names=T)
stims   <- stims0$stims ## 3 dim array: x-axis x y-axis x stimuli
rm(stims0); gc(F,T)
# list of associated filenames for each stimulus
# this is for ordering the categories properly
fnames  <- lapply(unique(catlist), function(category) {
  list.files(file.path("/data1/ffg05/data/stimuli", category), pattern="tiff")
})
fnames  <- unlist(fnames)


#' ## Layer 1: ICA simple/complex cells
#' 
#' This is all run in matlab and read in here
library(R.matlab)
layer1.feats <- readMat("/data1/ffg05/analysis/comp_model/step1_v1_features.mat")
layer1.feats <- layer1.feats$features


#' ## Select Data
#' 
#' Some of the features in `layer1.feats` might be constant or 0 across the 
#' different subjects, suggesting there isn't any meaningful information in 
#' those features. So we remove those features here. 
#' 
#+ load-remove-feats
require(Rcpp)
require(RcppArmadillo)
require(inline)
require(bigmemory)
# we first need to load all the inline plugins to run the c code
l <- getPlugin("RcppArmadillo")
plugin_bigmemory <- Rcpp::Rcpp.plugin.maker(
  include.before = "#include <RcppArmadillo.h>",  
  include.after = '
  #include "bigmemory/BigMatrix.h"     
  #include "bigmemory/MatrixAccessor.hpp"     
  #include "bigmemory/bigmemoryDefines.h"     
  #include "bigmemory/isna.hpp"     
  ', 
  libs    = "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)", 
  LinkingTo = c("BH", "bigmemory", "RcppArmadillo", "Rcpp"), 
  Depends = c("bigmemory", "RcppArmadillo", "Rcpp"), 
  package = "bigmemory"
)
inline::registerPlugin("plug_bigmemory", plugin_bigmemory)
getPlugin("plug_bigmemory")
# now we can call our simple getting the SD of the rows in a matrix
fun <- '
/* R -> C++ -> Arma */
Rcpp::NumericMatrix r_x(R_x);
arma::mat X(r_x.begin(), r_x.nrow(), r_x.ncol(), false);
/* Get row sds */
arma::vec rX = arma::stddev(X,0,1);
return Rcpp::wrap(rX);'
cpp_rowsd <- cxxfunction(signature(R_x = "numeric"), fun, 
                         plugin = "plug_bigmemory")
## can double check
## mat <- matrix(rnorm(10*100), 10, 100)
## res1 <- cpp_rowsd(mat)
## res2 <- apply(mat, 1, sd)
## all.equal(as.numeric(res1),as.numeric(res2))
# finally can run and reduce the layers1.feat matrix
tmp.sd <- cpp_rowsd(layer1.feats) # quite a bit faster than the alternative
layer1.feats <- layer1.feats[tmp.sd!=0,]


#' Note that we won't be examining information related to tools or circles 
#' so we shall remove it here.
#' 
#' Should I also remove TOOLS???
#' 
#+ load-remove-circles
layer1.feats <- layer1.feats[,!(catlist %in% c("Circles", "Tools"))]
stims        <- stims[,,!(catlist %in% c("Circles", "Tools"))]
fnames       <- fnames[!(catlist %in% c("Circles", "Tools"))]
fprefixes    <- sub("[.]tiff", "", fnames)
catlist      <- catlist[!(catlist %in% c("Circles", "Tools"))]

#' We will also load the timing information for each trial with the associated
#' stimulus file prefix. Since the order of the stimuli in this timing matrix 
#' is different than that in our `stims` array, we will also add a column that 
#' indicates what index that particular trial's stimulus can be found in the 
#' `stims` array.
#' 
#+ load-timing
timing <- read.csv("/data1/ffg05/notes/timing.csv")[,-1]
timing <- timing[!(timing$category %in% c("Circles", "Objects")),] ## remove circle odd-ball trials
timing$ind <- sapply(as.character(timing$name), function(x) which(x==fprefixes)) ## index of stimulus for trial in stims
all.equal(as.character(timing$name), fprefixes[timing$ind]) ## double check


#' ## Layer 2: PCA grouping features
#' 
#' This is really a dimensionality reduction of the features (rows)
#' Here we take only the top # of features (here 50)
#' We use this fast version of SVD based off netflix competition
#' ...http://cran.r-project.org/web/packages/irlba/vignettes/irlba.pdf
#' stats:::prcomp.default
library(irlba)
system.time(layer2.svd <- irlba(layer1.feats, nv=50, nu=0))
layer2.feats <- layer2.svd$v
colnames(layer2.feats) <- sprintf("comp%02i", 1:ncol(layer2.feats))
layer2.sdev  <- layer2.svd$d/sqrt(max(1, nrow(layer1.feats) - 1))
# TODO: autoselect number of components

#' ## Save
#' 
#+ save
outdir <- "/data1/ffg05/analysis/comp_model/v1_model"
dir.create(outdir, showWarnings=F)
save(stims, catlist, timing, layer1.feats, file=file.path(outdir, "layer1.rda"))
save(stims, catlist, timing, layer2.sdev, layer2.feats, file=file.path(outdir, "layer2.rda"))

#' Try more/less PCA components
system.time(layer2.svd2 <- irlba(layer1.feats, nv=100, nu=0))
layer2.feats2 <- layer2.svd2$v
colnames(layer2.feats2) <- sprintf("comp%02i", 1:ncol(layer2.feats2))
layer2.sdev2  <- layer2.svd2$d/sqrt(max(1, nrow(layer1.feats) - 1))
save(stims, catlist, timing, layer2.sdev2, layer2.feats2, file=file.path(outdir, "layer2_100.rda"))

system.time(layer2.svd3 <- irlba(layer1.feats, nv=25, nu=0))
layer2.feats3 <- layer2.svd3$v
colnames(layer2.feats3) <- sprintf("comp%02i", 1:ncol(layer2.feats3))
layer2.sdev3  <- layer2.svd3$d/sqrt(max(1, nrow(layer1.feats) - 1))
save(stims, catlist, timing, layer2.sdev3, layer2.feats3, file=file.path(outdir, "layer2_25.rda"))
