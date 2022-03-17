#!/bin/bash

fmap_path=$1
QSM_path=$2
part=$3
channels=$4



#for i in {01..$channels}
#i=01
#while [ "$i" -le "$channels" ];
#for ((j = 01; j <= $channels; j++ ));
for i in $(seq -f "%02g" 1 $channels)
    do

    echo $i
    offset="offset_"${part}"00"
    warp="_Warped"
    ending=".nii"

    antsApplyTransforms -d 3 \
    -i $fmap_path/$offset"$i"$ending \
    -r $QSM_path/src/mag_1_mean.nii \
    -o $QSM_path/$offset"$i"$warp$ending".gz" \
    -n BSpline \
    -t [$fmap_path/src/mag_1_mean_N4_N4_ants0GenericAffine.mat]

#     i=$((i + 01))
 done

ending=".nii.gz"
fslmerge -t $QSM_path/offset_${part}_Warped.nii.gz $QSM_path/$offset*$warp$ending

gunzip -f $QSM_path/offset_${part}_Warped.nii.gz

rm -r $fmap_path/$offset*
rm -r $QSM_path/$offset*
