#!/usr/bin/env bash

# This will transform the neurosynth functional rois from standard to functional space
#
# If you want to run in parallel:
# subjs=( $(cat ../sublist.txt) )
# parallel --no-notice -j 6 --eta ./10_neurosynth_std2sub.bash {} ::: ${subjs[@]}


run() {
  echo "$@"
  eval "$@"
  return $?
}


###
# Read in args
###

if [[ $# -ne 1 ]]; then
  echo "usage: $0 subject"
  echo "consider:"
  head -n 5 ../sublist.txt
  exit 2
fi

subject="$1"


###
# PATHS
###

base="/data1/ffg05"
roidir="${base}/analysis/rois"
runtype='Categories'


###
# Transform from standard to functional space
###

echo "=== TRANSFORM ==="

infile="${roidir}/group/neurosynth_visual_rois.nii.gz"
funcdir="${base}/analysis/preprocessed/${subject}/${runtype}"

outdir="${roidir}/${subject}"

mkdir -p ${outdir} 2> /dev/null

run "gen_applywarp.rb --overwrite -i ${infile} -r ${funcdir}/reg -w 'standard-to-exfunc' -o ${outdir}/neurosynth_visual_funcspace.nii.gz --interp nn --mask ${funcdir}/mask.nii.gz"

