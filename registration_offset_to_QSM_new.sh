#!/bin/bash

fmap_path=$1
QSM_path=$2
#offsets_real=$3
#offsets_imaginary=$4

N4BiasFieldCorrection -d 3 -i $fmap_path/src/mag_1_mean.nii -o $fmap_path/src/mag_1_mean_N4.nii
bet2 $fmap_path/src/mag_1_mean_N4.nii $fmap_path/src/mag_1_mean_N4_bet -f 0.7 -m
N4BiasFieldCorrection -d 3 -i $fmap_path/src/mag_1_mean.nii -x $fmap_path/src/mag_1_mean_N4_bet_mask.nii.gz -o $fmap_path/src/mag_1_mean_N4_N4.nii
   
N4BiasFieldCorrection -d 3 -i $QSM_path/src/mag_1_mean.nii -o $QSM_path/src/mag_1_mean_N4.nii
bet2 $QSM_path/src/mag_1_mean_N4.nii $QSM_path/src/mag_1_mean_N4_bet -f 0.2 -m
N4BiasFieldCorrection -d 3 -i $QSM_path/src/mag_1_mean.nii -x $QSM_path/src/mag_1_mean_N4_bet_mask.nii.gz -o $QSM_path/src/mag_1_mean_N4_N4.nii
   


antsRegistration -d 3 \
-o [$fmap_path/src/mag_1_mean_N4_N4_ants, $fmap_path/src/mag_1_mean_N4_N4_ants_Warped.nii.gz] \
--interpolation Linear --winsorize-image-intensities [0.005,0.995] --use-histogram-matching 0 \
--transform Translation[0.05] \
--metric MI[$QSM_path/src/mag_1_mean_N4_N4.nii, $fmap_path/src/mag_1_mean_N4_N4.nii,1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox \
--transform Rigid[0.05] \
--metric MI[$QSM_path/src/mag_1_mean_N4_N4.nii, $fmap_path/src/mag_1_mean_N4_N4.nii,1,32,Regular,0.25] \
--convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox -v

#echo $offsets_real
#echo $offsets_imaginary


#antsApplyTransforms -d 4 -e 3 \
#-i $fmap_path/$offsets_real \
#-r $QSM_path/src/mag_1_mean.nii.gz \
#-o $QSM_path/${offsets_real%%.nii}_Warped.nii.gz \
#-n Linear \
#-t [$fmap_path/src/mag_1_mean_N4_ants0GenericAffine.mat]


#antsApplyTransforms -d 4 -e 3\
#-i $fmap_path/$offsets_imaginary \
#-r $QSM_path/src/mag_1_mean.nii.gz \
#-o $QSM_path/${offsets_imaginary%%.nii}_Warped.nii.gz \
#-n Linear \
#-t [$fmap_path/src/mag_1_mean_N4_ants0GenericAffine.mat]
