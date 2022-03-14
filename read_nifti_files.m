function [ph, mag, TE, vox, matrix_size,imsize, NumEcho] = read_nifti_files(path,  subject, session, anat_or_fmap, file_ending, options)
addpath(genpath('/NAS/home/jhuck/Documents/QSM_toolboxes/Hongfu_Sun/QSM'))
addpath(genpath('/NAS/home/jhuck/Documents/2021_preventAD/scripts'))
addpath(genpath('/NAS/home/jhuck/Documents/QSM_toolboxes/nifti_utils/'))


if options.QSM_folder == 1
    path_nifti= strcat(path, subject, '/' , session, '/',anat_or_fmap, '/', options.QSM_folder_name,'/');
else
    path_nifti= strcat(path, subject, '/' , session, '/',anat_or_fmap, '/');
end
%% Phase and magnitude images should be in separate folders!!!!
phase_nifti_path=strcat(path_nifti)
magnitude_nifti_path=strcat(path_nifti)


%%
wrapped_list  = dir(fullfile(phase_nifti_path, '*.nii.gz'));
wrapped_name = wrapped_list(1).name;
ph1 =  load_untouch_nii([phase_nifti_path wrapped_name]).hdr.dime;

%ph_TE_list = dir(fullfile(unwrapped_path, '*_ph.nii.gz'));

NumEcho = ph1.dim(1);

if length(wrapped_list) > 10
    channels = 'uncombined';
    NumEcho = round(length(wrapped_list)/64);
else
    channels = 'combined';
    NumEcho = round(length(wrapped_list)/2);
end


    if strcmpi('uncombined',channels)   
            %%% read all phase echos
            rctr = 0; ictr=0;
            if NumEcho>1
                matrix_size = [ph1.dim(2);ph1.dim(3);ph1.dim(4)];
                vox = [ph1.pixdim(2);ph1.pixdim(3);ph1.pixdim(4)];

                n =  [ph1.dim(2) ph1.dim(3) ph1.dim(4)];
                N_std=ones(n)*(.1); 
                TE = single(zeros([NumEcho 1]));

                imsize=[matrix_size(1) matrix_size(2) matrix_size(3) NumEcho 32];
                ph = (zeros(imsize));
                mag = (zeros(imsize));
                for echo = 1:NumEcho
                    for channel = 1:32
                        sdir_phase = strcat('*echo-',num2str(echo),'_cH-',num2str(channel),'_part-phase_',file_ending,'.nii.gz')
                        phase_file = dir(fullfile(phase_nifti_path, sdir_phase));
                        wrapped_name = phase_file.name;
                        phase_tmp = nifti_utils.load_untouch_nii_vol_RAS([phase_nifti_path wrapped_name],'double');
                        
                        max_value = max(phase_tmp(:));
                        min_value = min(phase_tmp(:));
                        sdir = strcat('*echo-',num2str(echo),'_cH-',num2str(channel),'_part-mag_',file_ending,'.nii.gz');
                        mag_file=dir(fullfile(magnitude_nifti_path, sdir));
                        mag_name = mag_file.name;
                        mag_tmp = nifti_utils.load_untouch_nii_vol_RAS([magnitude_nifti_path mag_name],'double');
                        
                        for slice = 1:matrix_size(3)
                            ph(:,:,slice,echo,channel) = 2*pi.*((phase_tmp(:,:,slice)-min_value)/(max_value-min_value))-pi;
                            mag(:,:,slice,echo,channel) = single(mag_tmp(:,:,slice));
                        end
                    end

                    ph_tmp = load_untouch_nii([phase_nifti_path wrapped_name]).hdr.hist;
                    TE_tmp = ph_tmp.descrip;
                    Te = strsplit(TE_tmp,';');
                    T2_2 = strsplit(Te{1},'=');

                    TE (echo)= str2num(T2_2{2}) * 1e-3;     
                end
            else
                matrix_size = [ph1.dim(2);ph1.dim(3);ph1.dim(4)];
                vox = [ph1.pixdim(2);ph1.pixdim(3);ph1.pixdim(4)];

                n =  [ph1.dim(2) ph1.dim(3) ph1.dim(4)];
                N_std=ones(n)*(.1); 
                TE = single(zeros([NumEcho 1]));

                imsize=[matrix_size(1) matrix_size(2) matrix_size(3) 32];
                ph = (zeros(imsize));
                mag = (zeros(imsize));
                for echo = 1:NumEcho
                    for channel = 1:32
                        sdir_phase = strcat('*_cH-',num2str(channel),'_part-phase_',file_ending,'.nii.gz')
                        phase_file = dir(fullfile(phase_nifti_path, sdir_phase));
                        wrapped_name = phase_file.name;
                        phase_tmp = nifti_utils.load_untouch_nii_vol_RAS([phase_nifti_path wrapped_name],'double');
                        
                        max_value = max(phase_tmp(:));
                        min_value = min(phase_tmp(:));
                        
                        sdir = strcat('*_cH-',num2str(channel),'_part-mag_',file_ending,'.nii.gz');
                        mag_file=dir(fullfile(magnitude_nifti_path, sdir));
                        mag_name = mag_file.name;
                        mag_tmp = nifti_utils.load_untouch_nii_vol_RAS([magnitude_nifti_path mag_name],'double');
                        
                        for slice = 1:matrix_size(3)
                            ph(:,:,slice,echo,channel) = 2*pi.*((phase_tmp(:,:,slice)-min_value)/(max_value-min_value))-pi;
                            mag(:,:,slice,channel) = single(mag_tmp(:,:,slice));
                        end
                    end

                    ph_tmp =  load_untouch_nii([phase_nifti_path wrapped_name]).hdr.hist;
                    TE_tmp = ph_tmp.descrip;
                    Te = strsplit(TE_tmp,';');
                    T2_2 = strsplit(Te{1},'=');

                    TE (echo)= str2num(T2_2{2}) * 1e-3; 
                end
            end

    elseif strcmpi('combined',channels)   
            matrix_size = [ph1.dim(2);ph1.dim(3);ph1.dim(4)];
            vox = [ph1.pixdim(2);ph1.pixdim(3);ph1.pixdim(4)];

            n =  [ph1.dim(2) ph1.dim(3) ph1.dim(4)];
            N_std=ones(n)*(.1); 
            TE = single(zeros([NumEcho 1]));

            imsize=[matrix_size(1) matrix_size(2) matrix_size(3) NumEcho];
            ph = (zeros(imsize));
            mag = (zeros(imsize));

            %%% read all phase echos
            rctr = 0; ictr=0;
            for echo = 1:NumEcho
                
                idx = matrix_size(3)*echo
                sdir_phase = strcat('*TE',num2str(echo),'.nii.gz')
                phase_file = dir(fullfile(phase_nifti_path, sdir_phase));
                wrapped_name = phase_file.name;
                phase_tmp = nifti_utils.load_untouch_nii_vol_RAS([phase_nifti_path wrapped_name],'double');

                max_value = max(phase_tmp(:));
                min_value = min(phase_tmp(:));

                sdir = strcat('*TE',num2str(echo),'.nii.gz');
                mag_file=dir(fullfile(magnitude_nifti_path, sdir));
                mag_name = mag_file.name;
                mag_tmp = nifti_utils.load_untouch_nii_vol_RAS([magnitude_nifti_path mag_name],'double');


                for slice = 1:matrix_size(3)
                    ph(:,:,slice,echo) = (single(phase_tmp(:,:,slice))*slope+intercept)/(max_value)*pi;
                    mag(:,:,slice,echo) = single(mag_tmp(:,:,slice));
                end

                ph_tmp =  load_untouch_nii([phase_nifti_path wrapped_name]).hdr.hist;
                TE_tmp = ph_tmp.descrip;
                Te = strsplit(TE_tmp,';');
                T2_2 = strsplit(Te{1},'=');

                TE (echo)= str2num(T2_2{2}) * 1e-3;     

            end
    end
    
