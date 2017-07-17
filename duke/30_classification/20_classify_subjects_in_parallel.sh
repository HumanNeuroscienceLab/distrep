#!/usr/bin/env bash

subjs=( $( cat ../sublist.txt ) )
runtype="Categories"
nthreads=5
njobs=6

parallel --no-notice -j $njobs --eta \
Rscript --vanilla xx_4cats_manualroi.R {} ${nthreads} ::: ${subjs[@]}
