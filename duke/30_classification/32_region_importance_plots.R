library(ggplot2)
library(RColorBrewer)

bar_theme <- function() {
  library(grid)
  # Colors
  color.background = "#F0F0F0"
  color.grid.major = "#D0D0D0"
  color.axis.text = "#535353"
  color.axis.title = "#535353"
  color.title = "#3C3C3C"
  
  # Begin construction of chart
  theme_bw() +
    
    # Set the entire chart region to a light gray color
    #theme(panel.background=element_rect(fill=color.background, color=color.background)) +
    #theme(plot.background=element_rect(fill=color.background, color=color.background)) +
    theme(panel.border=element_rect(color="grey")) +
    
    # Format the grid
    theme(panel.grid.major.y=element_line(color=color.grid.major,size=.75)) +
    theme(panel.grid.major.x=element_blank()) +
    theme(panel.grid.minor=element_blank()) +
    theme(axis.ticks.x=element_blank()) +
    
    # Format the legend, but hide by default
    #theme(legend.position="none") +
    #theme(legend.background = element_rect(fill=color.background)) +
    theme(legend.text = element_text(size=16,color=color.axis.title)) +
    
    # Set title and axis labels, and format these and tick marks
    theme(plot.title=element_text(face="bold",color=color.title,size=24,vjust=2)) + #hjust=-0.08
    theme(axis.text.x=element_text(face="bold",size=17,color=color.axis.text)) +
    theme(axis.text.y=element_text(face="bold",size=17,color=color.axis.text)) +
    theme(axis.title.x=element_text(face="bold",size=18,color=color.axis.title, vjust=-.5)) +
    theme(axis.title.y=element_text(face="bold",size=18,color=color.axis.title, vjust=1.5)) +
    
    # Plot margins
    theme(plot.margin = unit(c(1, 1, .7, .7), "cm"))
}

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## Norms the data within specified groups in a data frame; it normalizes each
## subject (identified by idvar) so that they have the same mean, within each group
## specified by betweenvars.
##   data: a data frame.
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   na.rm: a boolean that indicates whether to ignore NA's
normDataWithin <- function(data=NULL, idvar, measurevar, betweenvars=NULL,
                           na.rm=FALSE, .drop=TRUE) {
  library(plyr)
  
  # Measure var on left, idvar + between vars on right of formula.
  data.subjMean <- ddply(data, c(idvar, betweenvars), .drop=.drop,
                         .fun = function(xx, col, na.rm) {
                           c(subjMean = mean(xx[,col], na.rm=na.rm))
                         },
                         measurevar,
                         na.rm
  )
  
  # Put the subject means with original data
  data <- merge(data, data.subjMean)
  
  # Get the normalized data in a new column
  measureNormedVar <- paste(measurevar, "_norm", sep="")
  data[,measureNormedVar] <- data[,measurevar] - data[,"subjMean"] +
    mean(data[,measurevar], na.rm=na.rm)
  
  # Remove this subject mean column
  data$subjMean <- NULL
  
  return(data)
}

## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {
  
  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))
  
  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }
  
  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL
  
  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)
  
  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")
  
  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)
  
  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )
  
  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor
  
  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}


#' ## Number of Subjects with Region
#' 
#' For this plot, I want to simply indicate the number of subjects with
#' that region.
#' 
#+ regions-in-reduced-models
indir   <- "/data1/ffg05/analysis/duke/classification"
infile  <- file.path(indir, "region_importance_reduced_model.csv")
imp <- read.csv(infile)
df <- reshape2::melt(imp/20*100) # convert number of subjects to a percentage
colnames(df) <- c("roi", "nsubjects")
df$roi <- sub("\\.", " ", df$roi)
df$roi <- factor(df$roi, levels=df$roi, labels=df$roi)

# plot
cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
p <- ggplot(df, aes(x=roi, y=nsubjects, fill=roi)) +
  geom_bar(position="dodge", stat="identity") + 
  coord_cartesian(ylim = c(0, 100)) + 
  scale_fill_manual(values=cols) + 
  ylab("Percent of Subjects") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank(), 
    legend.position="none"
  )
p

# save
outdir  <- "/data1/ffg05/analysis/duke/classification"
plotdir <- file.path(outdir, "plots")
outfile <- file.path(plotdir, "region_importance_reduced_model.pdf")
ggsave(outfile, p, width=5, height=3)


#' ## Region Importance
#' 
#' And now for this plot, we want to get the increase in accuracy by including
#' that region in the model.
#' 
#+ accuracy-boost-with-region

## NOTE: SEE 36_region_imp_by_cat_plot.R where this is redone (reimagined)

indir   <- "/data1/ffg05/analysis/duke/classification"
infile <- file.path(indir, "region_importance_leave_one_out.csv")
accuracies <- read.csv(infile)

change.accs <- as.matrix(sweep(accuracies[,-1], 1, accuracies[,1]))
change.accs <- change.accs*-1 # flip it  

srois <- sub("\\.", " ", colnames(change.accs))
srois <- factor(srois, levels=srois, labels=srois)

dimnames(change.accs) <- list(subject=1:nrow(change.accs), roi=srois) # this is so reshape uses the correct labels
dfw_long <- reshape2::melt(change.accs*100, value.name="change.accuracy")
dfwc <- summarySEwithin(dfw_long, measurevar="change.accuracy", withinvars="roi",
                        idvar="subject", na.rm=FALSE, conf.interval=.95)
dfwc

max.ylim <- (max(dfwc$change.accuracy)+max(dfwc$se))*1.02
# plot
cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
p <- ggplot(dfwc, aes(x=roi, y=change.accuracy, fill=roi)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=change.accuracy-se, ymax=change.accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim = c(-2, max.ylim)) + 
  scale_fill_manual(values=cols) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank(), 
    legend.position="none"
  )
# plot
plot(p)
# save
outfile <- file.path(plotdir, "region_importance_leave_one_out.pdf")
ggsave(outfile, p, width=5, height=3)


#' ## Individual Region Accuracy
#' 
#+ accuracy-individual-region
sroi.names <- c("R FFA", "PPA", "LOC", "L VWF", "V1", "M1", "A1")

indir   <- "/data1/ffg05/analysis/duke/classification"
infile <- file.path(indir, "region_importance_individual_accuracy.csv")
indiv.df <- read.csv(infile)
indiv.df$roi <- factor(as.character(indiv.df$roi), levels=sroi.names, 
                       labels=sroi.names)
indiv.df$accuracy <- indiv.df$accuracy*100

# Summarize input subjects => group average
dfwc <- summarySEwithin(indiv.df, measurevar="accuracy", withinvars="roi", 
                        idvar="subject", na.rm=F, conf.interval=.95)
dfwc

max.ylim <- (max(dfwc$accuracy)+max(dfwc$se))*1.02
# plot
cols <- brewer.pal(12, "Set3")[c(6:7,10,1,11,8,3)] # ok choose this finally
p <- ggplot(dfwc, aes(x=roi, y=accuracy, fill=roi)) +
  geom_bar(position="dodge", stat="identity") + 
  geom_errorbar(aes(ymin=accuracy-se, ymax=accuracy+se),
                width=.1,
                position=position_dodge(.9)) + 
  coord_cartesian(ylim = c(25, max.ylim)) + 
  scale_fill_manual(values=cols) + 
  ylab("Classification Accuracy") + 
  bar_theme() + 
  theme(
    axis.title.x = element_blank(), 
    legend.position="none"
  )
# plot
plot(p)
# save
outdir  <- "/data1/ffg05/analysis/duke/classification"
plotdir <- file.path(outdir, "plots")
outfile <- file.path(plotdir, "region_importance_individual_accuracy.pdf")
ggsave(outfile, p, width=5, height=3)

