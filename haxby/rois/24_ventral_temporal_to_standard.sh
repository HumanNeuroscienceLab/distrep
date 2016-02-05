#!/usr/bin/env bash

run() {
  echo "$@"
  eval "$@"
  return $?
}

mni="/mnt/nfs/psych/faceMemoryMRI/analysis/groups/mni152"
std="${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz"

prebase="/mnt/nfs/psych4/haxby01/analysis/preprocessed"
indir="/mnt/nfs/psych4/haxby01/rois/mni152"
outbase="/mnt/nfs/psych4/haxby01/rois/mni2func"
mkdir ${outbase} 2> /dev/null

hemis="lh rh"
rois="ventral_temporal lateral_temporal ventrolateral_occipital"
subjects=$( cat ../sublist.txt )
for subject in ${subjects}; do
  echo "$subject"
  
  fdir="${prebase}/${subject}/func"
  rdir="${fdir}/reg"
  
  for roi in ${rois}; do
    for hemi in ${hemis}; do
      i_roifile="${indir}/${hemi}_${roi}_maskslices.nii.gz"
      o_roifile="${outbase}/${subject}/${hemi}_${roi}.nii.gz"
      mkdir ${outbase}/${subject} 2> /dev/null
      run "gen_applywarp.rb --overwrite -i ${i_roifile} -r ${rdir} -w 'standard-to-exfunc' -o ${o_roifile} -m ${fdir}/mask.nii.gz --interp nn"
    done
    run "3dMean -overwrite -mask_union -prefix ${outbase}/${subject}/${roi}.nii.gz ${outbase}/${subject}/*_${roi}.nii.gz"
  done
  
  run "ln -sf ${fdir}/mean_func.nii.gz ${outbase}/${subject}/mean_func.nii.gz"
done
