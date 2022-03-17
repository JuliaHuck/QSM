addpath(genpath('~/Documents/QSM_coil_combination_with_low_res_fieldmaps'))

options.readout = 'bipolar';

niftiFolder = '/NAS/home/jhuck/Documents/2021_preventAD/nifti_new/'; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
% Get a list of all files and folders in this folder.
files = dir(niftiFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subj_ID = {subFolders(3:end).name} % Start at 3 to skip . and ..

for i = 1:length(subj_ID)
   % if ( ~ exist(strcat('/NAS/home/jhuck/Documents/2021_preventAD/nifti/', 'sub-', num2str(subj_ID(i)), '/ses-' , sessions{i},'/anat', '/QSM/corrected_phase/mag_cbm_smoothed_ppm.nii')))
    subj_ID{i}
        
    subjFolder = strcat(niftiFolder,'/',subj_ID{i}); % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
    % Get a list of all files and folders in this folder.
    files = dir(subjFolder);
    % Get a logical vector that tells which is a directory.
    dirFlags = [files.isdir];
    % Extract only those that are directories.
    subFolders = files(dirFlags); % A structure with extra info.
    % Get only the folder names into a cell array.
    sessionFolderNames = {subFolders(3:end).name} % Start at 3 to skip . and ..
    %% only run if combined phase doesn't excist
    [pathstr, ~, ~] = fileparts(which(strcat(subj_ID{i},'_',sessionFolderNames{1},'_ph-cbm.nii')));
    combined_phase = strcat(pathstr,'/', subj_ID{i},'_',sessionFolderNames{1},'_ph-cbm.nii')
    if (~ exist(combined_phase))
        coil_combination(niftiFolder, num2str(subj_ID{i}), sessionFolderNames{1}, options)
    end
    
end
