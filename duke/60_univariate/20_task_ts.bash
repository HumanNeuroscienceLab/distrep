#!/usr/bin/env bash

#subject="44172"
subjects=( $( cat /data1/ffg05/scripts/distrep/duke/sublist.txt ) )

# Basic paths
basedir="/data1/ffg05/analysis"
roidir="${basedir}/rois"
taskdir="${basedir}/task_activity"
predir="${basedir}/preprocessed"
tdir="/data1/ffg05/command/timing/afni"

for subject in ${subjects[@]}; do
  echo
  echo ${subject}
  
  # Get the input files
  roi="${roidir}/${subject}/neurosynth_visual_funcspace.nii.gz"
  mask="${predir}/${subject}/Categories/mask.nii.gz"
  func="${predir}/${subject}/Categories/filtered_func_data.nii.gz"
  motion="${predir}/${subject}/Categories/motion.1D"

  # Get the input/output files
  iodir="${taskdir}/${subject}/Categories/task_activity.rois"

  # Run
  #-input1D neurosynth_rois.1D \
  #-TR_1D 1.5 \
  cd ${iodir}
  3dDeconvolve \
      -local_times \
      -input ${predir}/${subject}/Categories/preproc/filtered_func_run*.nii.gz \
      -force_TR 1.5 \
      -polort 0 \
      -num_stimts 12 \
      -stim_times 1 "${tdir}/all_circles.1D" 'SPMG1' \
      -stim_label 1 circles \
      -stim_times 2 "${tdir}/all_faces.1D" 'SPMG1' \
      -stim_label 2 faces \
      -stim_times 3 "${tdir}/all_fruits.1D" 'SPMG1' \
      -stim_label 3 fruits \
      -stim_times 4 "${tdir}/all_letterstrings.1D" 'SPMG1' \
      -stim_label 4 letters \
      -stim_times 5 "${tdir}/all_tools.1D" 'SPMG1' \
      -stim_label 5 objects \
      -stim_times 6 "${tdir}/all_vehicles.1D" 'SPMG1' \
      -stim_label 6 vehicles \
      -stim_file 7 ${motion}'[0]' \
      -stim_base 7 \
      -stim_label 7 roll \
      -stim_file 8 ${motion}'[1]' \
      -stim_base 8 \
      -stim_label 8 pitch \
      -stim_file 9 ${motion}'[2]' \
      -stim_base 9 \
      -stim_label 9 yaw \
      -stim_file 10 ${motion}'[3]' \
      -stim_base 10 \
      -stim_label 10 dS \
      -stim_file 11 ${motion}'[4]' \
      -stim_base 11 \
      -stim_label 11 dL \
      -stim_file 12 ${motion}'[5]' \
      -stim_base 12 \
      -stim_label 12 dP \
      -num_glt 6 \
      -glt_label 1 circles_gt_all \
      -gltsym 'SYM: +4*circles -faces -fruits -letters -vehicles' \
      -glt_label 2 faces_gt_all \
      -gltsym 'SYM: +4*faces -circles -fruits -letters -vehicles' \
      -glt_label 3 fruits_gt_all \
      -gltsym 'SYM: +4*fruits -circles -faces -letters -vehicles' \
      -glt_label 4 letters_gt_all \
      -gltsym 'SYM: +4*letters -circles -fruits -faces -vehicles' \
      -glt_label 5 objects_gt_all \
      -gltsym 'SYM: +5*objects -circles -fruits -letters -faces -vehicles' \
      -glt_label 6 vehicles_gt_all \
      -gltsym 'SYM: +4*vehicles -circles -fruits -letters -faces' \
      -noFDR \
      -nobucket \
      -x1D xmat.1D \
      -xjpeg xmat.jpg \
      -x1D_stop
  
  3dREMLfit -matrix xmat.1D \
      -input neurosynth_rois.1D\' \
      -tout -noFDR \
      -Rbuck stats_bucket.1D \
      -verb
  
  echo
done
