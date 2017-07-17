#!/usr/bin/env bash

# This script will register all the zstats for each subject into 3mm space

subjects=( $( cat /data1/ffg05/scripts/distrep/duke/sublist.txt ) )
categories=( faces letters fruits vehicles )

basedir="/data1/ffg05/analysis/task_activity"
outdir="/data1/ffg05/analysis/combined/task_std_3mm"
mkdir ${outdir} 2> /dev/null

reffile="/data1/ffg05/analysis/rois/group/neurosynth_visual_rois_3mm.nii.gz"

for subject in ${subjects[@]}; do
  for category in ${categories[@]}; do
    echo "$subject - $category"
    indir="${basedir}/${subject}/Categories/task_activity_spmg1.reml"
    infile="${indir}/stats/zstat_${category}_gt_all.nii.gz"
    regdir="${indir}/reg"
    outfile="${outdir}/${subject}_zstat_${category}_gt_all.nii.gz"
    gen_applywarp.rb -i ${infile} -r ${regdir} -w 'exfunc-to-standard' -m ${reffile} -o ${outfile}
  done
  echo
done

#/data1/ffg05/analysis/task_activity/44172/Categories/task_activity_spmg1.reml/stats/zstat_faces_gt_all.nii.gz
