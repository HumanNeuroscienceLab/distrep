#!/usr/bin/env bash

# In the Haxby paper, 
# the volumes of interest (VOIs) are all in Talairach coordinates
# and posterior to the AC (or up to the AC).
# 
# For each VOI:
# * ventral temporal cortex is 70 to 20mm
# * lateral temporal cortex is also 70 to 20mm
# * ventrolateral occipital cortex is occipital pole to 70mm

# I found out that -70 from tal to mni is -72
# and -20 from tal to mni is still -20

# First I need to figure out the equivalence from talairach to mni
# Second I need to find out how to get these slices blocks

run() {
  echo "$@"
  eval "$@"
  return $?
}

#mni="/mnt/nfs/psych/faceMemoryMRI/analysis/groups/mni152"
std="${FSLDIR}/data/standard/MNI152_T1_2mm_brain.nii.gz"
outdir="/data1/ffg05/analysis/rois/group"
mkdir ${outdir} 2> /dev/null
ow=" -overwrite"

echo
echo "Only take y values above -70 to -20"
run "3dcalc${ow} -a ${std} -expr 'step(y-19)*step(72-y)' -prefix ${outdir}/slices_72to20_mask.nii.gz"
echo
echo "Only take y values above -72 to -108"
run "3dcalc${ow} -a ${std} -expr 'step(y-71)*step(108-y)' -prefix ${outdir}/slices_108to72_mask.nii.gz"
