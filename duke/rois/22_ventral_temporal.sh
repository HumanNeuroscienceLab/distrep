#!/usr/bin/env bash

# From the haxby paper
# --------------------
# Volumes of interest (VOI) were drawn on the high-resolution structural images to identify 
# ventral temporal, lateral temporal, and ventrolateral occipital cortex.
# 
# The VOI for ventral temporal cortex extended from 70 to 20 mm 
# and consisted of lingual, parahippocampal, fusiform, & inferior temporal gyri
# 
# The VOI for lateral temporal cortex also extended from 70 to 20 mm 
# and consisted of middle temporal gyrus & both banks of the superior temporal sulcus
# 
# The VOI for ventrolateral occipital cortex extended from the occipital pole to 70 mm 
# and consisted of lingual, fusiform, inferior occipital, & middle occipital gyri


run() {
  echo "$@"
  eval "$@"
  return $?
}

mni="/mnt/nfs/psych/faceMemoryMRI/analysis/groups/mni152"
std="${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz"
outbase="/data1/ffg05/analysis/rois/group"
outdir="${outbase}/from_aparc_DKTatlas40"
mkdir ${outdir} 2> /dev/null
ow=" -overwrite"

echo

indir="${mni}/freesurfer/aparc_DKTatlas40"

echo
echo "For the VTC"
echo "First use freesurfer segmented MNI152 brain to get each region"
echo "Then dilate each region and fill holes"
regions="lingual parahippocampal fusiform inferiortemporal"
lh_outs=""; rh_outs=""
for region in $regions; do
run "3dresample${ow} -inset ${indir}/lh_${region}.nii.gz -master ${std} -prefix ${outdir}/lh_${region}.nii.gz"
run "3dresample${ow} -inset ${indir}/rh_${region}.nii.gz -master ${std} -prefix ${outdir}/rh_${region}.nii.gz"

run "3dmask_tool -overwrite -inputs ${outdir}/lh_${region}.nii.gz -dilate_input 2 -1 -fill_holes -prefix ${outdir}/lh_${region}.nii.gz"
run "3dmask_tool -overwrite -inputs ${outdir}/rh_${region}.nii.gz -dilate_input 2 -1 -fill_holes -prefix ${outdir}/rh_${region}.nii.gz"

lh_outs="${lh_outs} ${outdir}/lh_${region}.nii.gz"
rh_outs="${lh_outs} ${outdir}/rh_${region}.nii.gz"
done
voi="ventral_temporal"
run "3dMean -mask_union -prefix ${outbase}/lh_${voi}.nii.gz${lh_outs}"
run "3dMean -mask_union -prefix ${outbase}/rh_${voi}.nii.gz${rh_outs}"
run "3dcalc -a ${outbase}/lh_${voi}.nii.gz -b ${outbase}/slices_72to20_mask.nii.gz -expr 'step(a)*step(b)' -prefix ${outbase}/lh_${voi}_maskslices.nii.gz"
run "3dcalc -a ${outbase}/rh_${voi}.nii.gz -b ${outbase}/slices_72to20_mask.nii.gz -expr 'step(a)*step(b)' -prefix ${outbase}/rh_${voi}_maskslices.nii.gz"
