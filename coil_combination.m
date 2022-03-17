function coil_combination(path, subject, session, options)
%parfor ex=1:32
addpath(genpath('~/Documents/nifti_utils/'))
addpath(genpath('~/Documents/QSM_toolboxes/Hongfu_Sun/QSM'))

%QSM_R2S_PRISMA Quantitative susceptibility mapping from R2s sequence at PRISMA (3T).
%   QSM_R2S_PRISMA(PATH_MAG, PATH_PH, PATH_OUT, OPTIONS) reconstructs susceptibility maps.
%
%   Re-define the following default settings if necessary
%
%   PATH         - directory of niftis
%   OPTIONS      - parameter structure including fields below
%    .readout    - multi-echo 'unipolar' or 'bipolar'        : 'unipolar'


if ~ exist('path','var') || isempty(path)
    error('Please input the directory of nifti files')
end

if ~ exist('options','var') || isempty(options)
    options = [];
end

if ~ isfield(options,'readout')
    options.readout = 'unipolar';
end

readout    = options.readout;

%% read in NIFTIs of both magnitude and raw unfiltered phase images
% low resolution phase images    
anat_or_fmap= 'fmap';
file_ending = 'fieldmap';
options.QSM_folder = 0;

files = dir(strcat(path, subject, '/' , session,'/' ,anat_or_fmap));

if length(files)>1

    [ph, mag, TE, vox, matrix_size, imsize, NumEcho] = read_nifti_files(path,  subject, session, anat_or_fmap, file_ending, options);


    % define output directories

    fmap_out = strcat(path, subject, '/' , session,'/' ,anat_or_fmap, '/QSM/');
    mkdir(fmap_out);
    fmap_out = strcat(path, subject, '/' , session,'/' ,anat_or_fmap, '/QSM/corrected_phase/');
    mkdir(fmap_out);
    init_dir = pwd;
    cd(fmap_out);

    %% save magnitude and raw phase nifties for each echo
    mkdir('src')
    for echo = 1:imsize(4)
        mag_tmp = permute(mag, [1,2,3,5,4]);
        nii = make_nii(mag_tmp(:,:,:,:,echo),vox);
        save_nii(nii,['src/mag_' num2str(echo) '.nii']);
        ph_tmp = permute(ph, [1,2,3,5,4]);
        nii = make_nii(ph_tmp(:,:,:,:,echo),vox);
        save_nii(nii,['src/ph_' num2str(echo) '.nii']);
    end

    nii = make_nii(mag(:,:,:,:,:),vox);
    save_nii(nii,['src/mag_all.nii']);
    nii = make_nii(ph(:,:,:,:,:),vox);
    save_nii(nii,['src/ph_all.nii']);

    unix('fslmaths src/mag_1.nii -Tmean src/mag_1_mean.nii.gz');
    unix('gunzip -f src/mag_1_mean.nii.gz')

    %%

    anat_or_fmap= 'anat'
    path_out = strcat(path, subject, '/' , session,'/' ,anat_or_fmap, '/QSM/');
    mkdir(path_out);
    path_out = strcat(path, subject, '/' , session,'/' ,anat_or_fmap, '/QSM/','corrected_phase/');
    mkdir(path_out);
    path_out_tmp = strcat(path,  subject, '/' , session,'/' ,anat_or_fmap, '/QSM/','corrected_phase/src');
    mkdir(path_out_tmp);
    init_dir = pwd;
    cd(path_out);
    
    files = dir(strcat(path, subject, '/' , session,'/' ,anat_or_fmap, '/QSM_uncombined'));

    if length(files)>1

        anat_or_fmap= 'anat';
        file_ending = 'GRE';
        options.QSM_folder = 1;
        options.QSM_folder_name = 'QSM_uncombined'

        [ph, mag, TE_anat, vox, matrix_size, imsize, NumEcho] = read_nifti_files(path,  subject, session, anat_or_fmap, file_ending, options);

        mkdir('src')

        nii = make_nii(mag(:,:,:,:),vox);
        save_nii(nii,['src/mag_1.nii']);
        nii = make_nii(ph(:,:,:,:),vox);
        save_nii(nii,['src/ph_1.nii']);

        %% brain extraction
        % generate mask from magnitude of the first echo

        unix('fslmaths src/mag_1.nii -Tmean src/mag_1_mean.nii.gz');
        unix('gunzip -f src/mag_1_mean.nii.gz')

        %% registration
        [pathstr, ~, ~] = fileparts(which('registration_offset_to_QSM_new.sh'));
        setenv('pathstr',pathstr);
        
        cmdStr = ['sh ${pathstr}/registration_offset_to_QSM_new.sh' ' ' fmap_out ' ' path_out ];
        system(cmdStr);

        %% phase offset correction
        % if unipolar 
        cd(fmap_out);

        name = strcat(fmap_out,'src/ph_all.nii');
        nii = load_nii(name);
        ph = double(nii.img);

        name = strcat(fmap_out, 'src/mag_all.nii');
        nii = load_nii(name);
        mag = double(nii.img);
        vox_new= [nii.hdr.dime.pixdim(2);nii.hdr.dime.pixdim(3);nii.hdr.dime.pixdim(4)];

        unix('gunzip -f src/mag_1_mean_N4_bet_mask.nii.gz')
        name = strcat(fmap_out,'/src/mag_1_mean_N4_bet_mask.nii');
        nii = load_nii(name);
        mask = double(nii.img);

        if strcmpi('unipolar',readout)
            [ph_corr, mag_cmb] = geme_cmb2(mag.*exp(1j*ph),vox,TE,mask);
            % if bipolar
        elseif strcmpi('bipolar',readout)
            ph_corr = zeros(imsize);

            if mod(NumEcho,2) == 0 
                ph_corr(:,:,:,1:2:end) = geme_cmb(mag(:,:,:,1:2:end).*exp(1j*ph(:,:,:,1:2:end)),vox,TE(1:2:end),mask);
                ph_corr(:,:,:,2:2:end) = geme_cmb(mag(:,:,:,2:2:end).*exp(1j*ph(:,:,:,2:2:end)),vox,TE(2:2:end),mask);
            else
                [ph_corr, mag_cmb, offset_real, offset_imag] = poem_bi3_new(mag, ph, vox_new, TE, mask);
            end
        else
            error('is the sequence unipolar or bipolar readout?')
        end

        imsize = size(ph_corr);
        nii = make_nii(ph_corr(:,:,:,:),vox_new);
        save_nii(nii,['ph_corr.nii']);

        nii = make_nii(mag_cmb(:,:,:,:),vox_new);
        save_nii(nii,['mag_cmb.nii']);

        imsize = size(offset_imag);
        nii = make_nii(offset_imag(:,:,:,:),vox_new);
        save_nii(nii,['offset_imaginary.nii']);

        nii = make_nii(offset_real(:,:,:,:),vox_new);
        save_nii(nii,['offset_real.nii']);

        %% register offset to QSM and apply defomration fields

        offset_real = load_nii('offset_real.nii').img;
        channelNbr = size(offset_real);
        for ch = 1:channelNbr(4)

            nii = make_nii(offset_real(:,:,:,ch),vox_new);
            if ch<10
                save_nii(nii,[strcat('offset_real000',num2str(ch),'.nii')]);
            else
                save_nii(nii,[strcat('offset_real00',num2str(ch),'.nii')]);
            end
        end
        [pathstr, ~, ~] = fileparts(which('apply_ants_reg.sh'));
        setenv('pathstr',pathstr);
        cmdStr = ['sh ${pathstr}/apply_ants_reg.sh' ' ' fmap_out ' ' path_out ' real ' num2str(channelNbr(4))];
        system(cmdStr);

        offset_imaginary= load_nii('offset_imaginary.nii').img;
        for ch = 1:channelNbr(4)

            nii = make_nii(offset_imaginary(:,:,:,ch),vox_new);
            if ch<10
                save_nii(nii,[strcat('offset_imaginary000',num2str(ch),'.nii')]);
            else
                save_nii(nii,[strcat('offset_imaginary00',num2str(ch),'.nii')]);
            end
        end
        cmdStr = ['sh ${pathstr}/apply_ants_reg.sh' ' ' fmap_out ' ' path_out ' imaginary ' num2str(channelNbr(4))];
        system(cmdStr);

        %% subtract offset from high resolution phase images

        offsets_real = nifti_utils.load_untouch_nii4D_vol_RAS([path_out '/offset_real_Warped.nii'], 'double');
        offsets_imaginary = nifti_utils.load_untouch_nii4D_vol_RAS([path_out '/offset_imaginary_Warped.nii'],'double');
        
        offsets = complex(offsets_real, offsets_imaginary);

        mag_tmp = load_nii([path_out_tmp '/mag_1.nii']);
        mag = double(mag_tmp.img);

        ph_tmp = load_nii([path_out_tmp '/ph_1.nii']);
        ph_corr = double(ph_tmp.img);
        vox = [ph_tmp.hdr.dime.pixdim(2);ph_tmp.hdr.dime.pixdim(3);ph_tmp.hdr.dime.pixdim(4)];

        img_cmb = mag.*exp(1j*ph_corr)./offsets;
        img_cmb(isnan(img_cmb)) = 0;
        ph_cmb = angle(sum(img_cmb,4));
        ph_cmb(isnan(ph_cmb)) = 0;
        mag_cmb = abs(mean(img_cmb,4));
        mag_cmb(isnan(mag_cmb)) = 0;

        %% save offset corrected and combined high resolution QSM images
        nii = make_nii(ph_cmb(:,:,:,:),vox);
        name = strcat(path_out, subject,'_',session,'_ph-cbm.nii')
        save_nii(nii,name);
        cd(path_out)

        nii = make_nii(mag_cmb(:,:,:,:),vox);
        name = strcat(path_out, subject,'_',session,'_mag-cmb_mean.nii')
        save_nii(nii,name);
        cd(path_out)


        %% QUality control of coil combination (set contrast to 0 -100)

        mag_sum = sum(mag,4);

        Q = ((abs(sum(img_cmb,4)))./mag_sum).*100;
        
        Q(Q>100) = 100;

        nii = make_nii(Q(:,:,:),vox);
        name = strcat(path_out, subject,'_',session,'_Quality-Control','.nii')
        save_nii(nii,name);

        unix('gunzip -f src/mag_1_mean_N4_bet_mask.nii.gz')
        movefile(strcat(path_out,'src/','mag_1_mean_N4_bet_mask.nii'),strcat(path_out,'src/',subject,'_',session,'_mag_1_mean_N4_bet_mask','.nii'));
        movefile(strcat(path_out,'src/','mag_1_mean_N4_bet.nii.gz'),strcat(path_out,'src/',subject,'_',session,'_mag_1_mean_N4_bet','.nii.gz'));

    end
    
    fprintf('No high reolution QSM images')

end
fprintf('No low reolution fieldmaps')

clear all;
