#' #### Overlap Between Voxels
#' 
#' We want to get and visualize the amount of overlap observed. This should 
#' honestly just be a table.
#' 
#+ save-overlap
res$best$overlap

#' Ok so if we compare this to chance, we find that the fact that just by chance
#' we can get a pretty similar overlap structure. This would suggest then, that
#' it isn't too special to get lots of non-zero feature weights that are only 
#' for one of the categories.
#' 
#+ perm-overlap
ref.df <- res$final.fits$df[res$best$stats$lind]
perm.overlap <- laply(1:100, function(i) {
  mod <- wrap.model(subset.grpdat, sample(ys), foldsI, alpha=alpha, nlambda=100, parallel=F)
  ind <- which.min(abs(mod$final.fits$df - ref.df))
  betas <- sapply(mod$final.fits$beta, function(beta) beta[,ind])
  crossprod(betas!=0)
}, .parallel=T)

# So we want to compare the permuted results to what we get here
ref.overlap <- res$best$overlap
round(apply(perm.overlap, 2:3, mean)) # the average null
round(apply(perm.overlap, 2:3, sd)) # the sd null
ref.overlap
tt <- laply(1:4, function(i) {
  laply(1:4, function(j) {
    mean(ref.overlap[i,j]>=c(perm.overlap[,i,j], ref.overlap[i,j]))
  })
})
tt
1-tt

#' The one last thing that we can check is the location of the non-zero feature
#' weights. Do the features in a given region tend to clump? I'm not totally 
#' sure how to phrase this question. Maybe one way would be if 
#' 

hdr <- read.nifti.header("/mnt/nfs/share/fsl/current/data/standard/MNI152_T1_3mm_brain.nii.gz")
## we subset our whole brain mask to our ROIs
subset.mask <- mask; subset.mask[!(mask %in% urois[1:4])] <- 0; subset.mask[subset.mask!=0] <- 1
## convert it into a 3D array
mask3d <- subset.mask; dim(mask3d) <- hdr$dim
## then get the array indices of each non-zero value
inds <- which(mask3d==1, arr.ind=T)
## finally we want to get the xyz coordinate given the ijk coordinate
ijkt <- cbind(inds, rep(1, nrow(inds)))
xyzt <- hdr$qto.xyz %*% t(ijkt)
xyz  <- t(xyzt[1:3,])
## and more finally we get the distance of these different voxels
betas <- res$best$betas
d    <- dist(inds)
dmat <- as.matrix(d)

# So here we got how close the features in the FFA are to each other for faces
sdmat <- dmat[subset.grouping==1 & betas[1,]!=0, subset.grouping==1 & betas[1,]!=0]
sdmat <- dmat[subset.grouping==1, subset.grouping==1]
mean(sdmat[lower.tri(sdmat)])
hist(sdmat[lower.tri(sdmat)])
hist()

# Can ask how likely it is that if you are within 3 voxel neighborhood, then it's good
nei.dist <- dist(rbind(c(3,3,3), c(0,0,0)))[1]
sdmat <- dmat[betas[4,]!=0, betas[4,]!=0]
mean(colSums(sdmat<=nei.dist))
mean(sdmat[lower.tri(sdmat)])

dim(hdr$qto.xyz %*% t(inds))

perm.betas <- laply(1:100, function(i) {
  mod <- wrap.model(subset.grpdat, sample(ys), foldsI, alpha=alpha, nlambda=100, parallel=F)
  ind <- which.min(abs(mod$final.fits$df - ref.df))
  betas <- sapply(mod$final.fits$beta, function(beta) beta[,ind])
  betas
}, .parallel=T)


# Here we want to see if on average the distance between the non-zero features
# says something about them clumping together or being very far apart.
# It's interesting to note that on average the letters category is pretty spread
# out.
# 
# Faces are marginally significant in terms of being close together as well as
# fruits....fruits is a bit surprising.
betas <- t(res$best$betas)
ref.zztops <- sapply(1:4, function(i) {
  sdmat <- dmat[betas[,i]!=0, betas[,i]!=0]
  #mean(colSums(sdmat<=nei.dist))
  mean(sdmat[lower.tri(sdmat)])
})

zztops <- aaply(perm.betas, .(1), function(betas) {
  sapply(1:4, function(i) {
    sdmat <- dmat[betas[,i]!=0, betas[,i]!=0]
    #mean(colSums(sdmat<=nei.dist))
    mean(sdmat[lower.tri(sdmat)])
  })
}, .parallel=T)

ref.zztops
sapply(1:4, function(i) mean(ref.zztops[i]<=c(zztops[,i], ref.zztops[i])))
sapply(1:4, function(i) mean(ref.zztops[i]>=c(zztops[,i], ref.zztops[i])))


# This measure is significant! This also should make more sense
# So is it the case that the non-zero features clump together?
# If they do, then we might expect that the non-zero features will be more
# likely to be close-by within some radius. Here we set the radius to 3 voxels
# We find that only for the faces, case is there significant amount of clumping.
betas <- t(res$best$betas)
ref.zztops <- sapply(1:4, function(i) {
  sdmat <- dmat[betas[,i]!=0, betas[,i]!=0]
  mean(colSums(sdmat<=nei.dist))
})

zztops <- aaply(perm.betas, .(1), function(betas) {
  sapply(1:4, function(i) {
    sdmat <- dmat[betas[,i]!=0, betas[,i]!=0]
    mean(colSums(sdmat<=nei.dist))
  })
}, .parallel=T)

ref.zztops
colMeans(zztops)
sapply(1:4, function(i) mean(ref.zztops[i]<=c(zztops[,i], ref.zztops[i])))
sapply(1:4, function(i) mean(ref.zztops[i]>=c(zztops[,i], ref.zztops[i])))

