# QSM_coil_combination_with_low_res_fieldmaps
Matlab code to combine multi channel phase array using low resolution fieldmaps as reference scan. The offset map of each channel is calculated from the low resolution fieldmaps. The fieldmaps are registered to the high resolution QSM magnitude images. The deformation fields of the registration are applied to the offset maps. Afterwards the phase images of the high resolution QSM images are corrected with the registered offset maps. 

Quality control images of the phase combination are saved in addition to the combined phase and magnitude images. 

The code requires two toolboxes:
Nifti images are read in with the nifti_utils toolbox from Justin Blaber:
https://github.com/justinblaber/nifti_utils.git

The offset maps are calcluated with the peom method from Hongfu Sun. However, the peom method had to be adapted to save the offset maps from the low resolution fieldmaps (peom_bi3_new.m and poem_new.m). These methods call functions from Hongfu Sun's toolbox:
https://github.com/sunhongfu/QSM.git
