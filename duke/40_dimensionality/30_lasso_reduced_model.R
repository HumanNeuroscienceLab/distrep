
indir <- "/data1/ffg05/analysis/duke/classification"
load(file.path(indir, "subjects_4region_models.rda"))

snames   <- sub("[LR] ", "", sroi.names[1:4])
snames[4]<- "VWFA"


lst.res[[2]]$res$best$category.stats$nfeats/sum(table(dat.subs[[2]]$grouping))


# for the individual region results
subjects <- as.character(read.table("/data1/ffg05/scripts/repsim/sublist.txt")[,])
sub.rets <- llply(subjects, function(subject) {
  subdir <- sprintf("/data1/ffg05/analysis/classification/%s", subject)
  load(file.path(subdir, "z2_bs_neurosynthvoxs_rfce_4cats_glmnet.rda"))
  list(res=res.rfces, grouping=vox.grouping)
}, .progress="text")

table(sub.rets[[2]]$grouping)


# reorder cuz right now: faces, fruits, letters, vehicles

# Individual ROI analysis
# This will give me the number of features across subjects
indiv.nfeats <- laply(1:length(subjects), function(si) {
  grouping <- dat.subs[[si]]$grouping
  x        <- sub.rets[[si]]$res$alpha1$individual
  sub.regions <- x$res[1:4]
  #sub.tot     <- table(grouping)
  ret  <- sapply(1:length(sub.regions), function(i) {
    x <- sub.regions[[i]]
    ret <- x$res$best$category.stats$nfeats[c(1,3,2,4)]
    names(ret) <- snames
    ret
  })
  colnames(ret) <- categories.title
  ret
})
#apply(indiv.nfeats, c(2,3), mean)
apply(indiv.nfeats, c(2,3), median)

#ret2  <- sapply(1:length(sub.regions), function(i) {
#  x <- sub.regions[[i]]
#  ret <- x$res$best$category.stats$nfeats[c(1,3,2,4)]
#  names(ret) <- snames
#  ret/sub.tot[i]
#})
#colSums(ret)/sum(sub.tot) # faces, letters, fruits, vehicles


# Combined ROI
combined.nfeats <- laply(1:length(subjects), function(si) {
  rois     <- dat.subs[[si]]$rois
  grouping <- dat.subs[[si]]$grouping
  cres     <- lst.res[[si]]$res
  
  cres$best$category.stats$nfeats
  ind <- cres$best$stats$lind
  
  # the beta will be the # of categories
  # for each element in beta, you'll have a matrix of nvoxs x nlambdas
  nfeats.mat <- sapply(cres$final.fits$beta, function(beta) {
    nfeats <- sapply(rois, function(roi) {
      sum(beta[grouping == roi,ind]!=0)
    })
    names(nfeats) <- snames
    nfeats
  })
  colnames(nfeats.mat) <- categories.title
  
  nfeats.mat
})
apply(combined.nfeats, c(2,3), median)



# Compare
apply(indiv.nfeats, c(2,3), median)
apply(combined.nfeats, c(2,3), median)


# DOESN'T WORK!
# Want to graph the percentage of voxels
# I could also go back to the dfs and get the numbers in each
# Then can do better direct comparison
library(ggplot2)
library(grid)
library(RColorBrewer)

source("duke/40_dimensionality/qa_plot_measures.R")


df.indiv <- reshape2::melt(indiv.nfeats)
colnames(df.indiv) <- c("subject", "roi", "category", "nfeats")
df.comb <- reshape2::melt(combined.nfeats)
colnames(df.comb) <- c("subject", "roi", "category", "nfeats")
df <- data.frame(
  data.frame(df.indiv, roi.type="individual"), 
  data.frame(df.comb, roi.type="combined")
)

p1 <- ggplot(df.indiv, aes(x=roi, y=nfeats, color=roi, fill=roi)) +
  geom_boxplot(lwd=1.5, width=0.5, position="dodge") + 
  stat_summary(geom="crossbar", width=0.3, fatten=2, color="white", fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) + 
  ylab("# of Components") +
  xlab("") + 
  theme_bw() + 
  theme(panel.border=element_blank()) +
  theme(panel.grid.major.x= element_blank()) + 
  theme(panel.grid.minor  = element_blank()) + 
  theme(panel.grid.major.y= element_line(color="grey70")) + 
  theme(axis.title.x      = element_blank()) +  
  theme(axis.title.y      = element_text(family = "Century Gothic", face = "plain", 
                                         size=18, angle=90, vjust=1)) +  
  theme(axis.text.x       = element_text(family = "Century Gothic", face = "plain", 
                                         size=16, 
                                         margin=margin(0.15,0.15,0.15,0.15,"lines"))) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5, 
                                         margin=margin(0.15,0.15,0.15,0.15,"lines"))) + 
  theme(axis.ticks.y = element_line(color="grey")) +
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.x      = element_blank()) +
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines")) + 
  theme(legend.position="none")
plot(p1)

p2 <- ggplot(df.comb, aes(x=roi, y=nfeats, color=roi, fill=roi)) +
  geom_boxplot(lwd=1.5, width=0.5, position="dodge") + 
  stat_summary(geom="crossbar", width=0.3, fatten=2, color="white", fun.data = function(x) c(y=median(x), ymin=median(x), ymax=median(x))) + 
  ylab("# of Components") +
  xlab("") + 
  theme_bw() + 
  theme(panel.border=element_blank()) +
  theme(panel.grid.major.x= element_blank()) + 
  theme(panel.grid.minor  = element_blank()) + 
  theme(panel.grid.major.y= element_line(color="grey70")) + 
  theme(axis.title.x      = element_blank()) +  
  theme(axis.title.y      = element_text(family = "Century Gothic", face = "plain", 
                                         size=18, angle=90, vjust=1)) +  
  theme(axis.text.x       = element_text(family = "Century Gothic", face = "plain", 
                                         size=16, 
                                         margin=margin(0.15,0.15,0.15,0.15,"lines"))) + 
  theme(axis.text.y       = element_text(family = "Times", face = "plain", 
                                         size=16, angle=0, hjust=0.5, 
                                         margin=margin(0.15,0.15,0.15,0.15,"lines"))) + 
  theme(axis.ticks.y = element_line(color="grey")) +
  theme(axis.ticks.length = unit(.15, "lines")) + 
  theme(axis.ticks.x      = element_blank()) +
  theme(plot.margin       = unit(c(1, 1, 0.25, 1), "lines")) + 
  theme(legend.position="none")
plot(p2)
#ggsave("anat_icc_btw.png", width=4, height=3, dpi=100)
