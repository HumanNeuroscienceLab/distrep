if (!any(.libPaths() == "/home/zshehzad/R_libs")) .libPaths(c("~/R_libs", .libPaths()))
library(methods)

library(plyr)
library(RColorBrewer)


subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])

# Test with one subject
subject <- subjects[2]
subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
res.rfces <- NULL; vox.grouping <- NULL; classes <- NULL; roi.df <- NULL
load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))

# regions found
use.rois <- res.rfces$alpha0.5$partial$rois
rinds    <- roi.df$vals %in% use.rois
as.character(roi.df$fullnames)[roi.df$vals %in% use.rois]

# get betas and rois
betas <- res.rfces$alpha0.5$partial$res$best$betas
rois  <- vox.grouping[vox.grouping %in% use.rois]

# rorder
oinds <- order(rois)
betas <- betas[,oinds]
rois  <- rois[oinds]

# plot 
dat.cols <- rev(brewer.pal(11, "RdYlBu"))
roi.df$cols <- c(brewer.pal(3, "Dark2")[1], 
                 brewer.pal(4, "Paired"), brewer.pal(3, "Dark2")[3], 
                 brewer.pal(10, "Paired")[5:10])
cc <- as.character(factor(rois, levels=roi.df$vals[rinds], 
                          labels=roi.df$cols[rinds]))

betas[betas==0] <- NA
rownames(betas) <- levels(classes)

heatmap(betas, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=dat.cols, zlim=c(-1,1))

x <- (betas>0) - (betas<0)
heatmap(x, Rowv=NA, Colv=NA, ColSideColors=cc, labCol=NA, 
        col=c("blue", "white","red"))



# More
rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")

accs <- laply(rets, function(x) {
  sapply(x$res, function(xx) xx$partial$acc)
})
table(apply(accs, 1, which.max))
t.test(accs[,2], accs[,3], paired=T)

rns <- llply(rets, function(x) {
  as.character(roi.df$fullnames[roi.df$vals %in% x$res$alpha0.5$partial$rois])
})
sort(table(unlist(rns)))





tmp <- laply(rets, function(sret) {
  g <- sret$grouping
  sbetas <- sign(sret$res$alpha0.5$full$res$best$betas)
  
  tmp <- laply(sort(unique(g)), function(ri) {
    laply(1:4, function(ci) {
      dat <- sbetas[ci,g==ri]
      cbind(zero=sum(dat==0), negative=sum(dat==(-1)), positive=sum(dat==1))
    })
  })
  rnames <- roi.df$names[roi.df$vals %in% unique(g)]
  dimnames(tmp) <- list(roi=rnames, category=categories, effect=dimnames(tmp)[[3]])
  
  return(tmp)
})
atmp <- apply(tmp, 2:4, mean) # average across subjects
