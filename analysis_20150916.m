%% MRI Map Analysis
% 
%  This code is for finding means of one or more regions of interest (ROIs)
%  in an image. It is intended to analyze a T2, T1, M0 Map. This is not for
%  processing spin echo data (i.e. raw images that still need to be plotted
%  and fitted to obtain the relaxation time and amplitude values). 
% 
%  EDITS
%  Adapted from t2mapanalysis_SE_MC_Shorter_better - Lina A. Colucci 
%  16 Sept 2015 - Adapted from gre_analysis.m - Lina A. Colucci 
%  ------------------------------------------------------------------------

clear; close all; 

%% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')

%% Load  NIFTI files
dir = ['/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Niftis/t2_map.nii']; % Where is image of interest located? 

% Shorter TE's
imageShortTE_struct = load_nifti(dir);
imageShortTE = squeeze(imageShortTE_struct.vol); 
size_short = size(imageShortTE);
nSlices_short = size_short(3);


%% Load Masks
mask_tube1 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/tube1_closest-to-skin.nii');
mask_tube2 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/tube2.nii');
mask_tube3 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/tube3.nii');
mask_tube4 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/tube4.nii');
mask_tube5 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/tube5_furthest-from-skin.nii');
mask_noise = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/noise.nii');
mask_subcu = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/subcu.nii');
mask_muscle = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/muscle.nii');

% *************************************************
slice=81; % which slice # do you want to analyze? 
% ^WHY IS IT in slice 81 instead of 80 like it was created in?????????????
% *************************************************

roi_tube1 = logical(mask_tube1.vol(:,:,slice));
roi_tube2 = logical(mask_tube2.vol(:,:,slice));
roi_tube3 = logical(mask_tube3.vol(:,:,slice));
roi_tube4 = logical(mask_tube4.vol(:,:,slice));
roi_tube5 = logical(mask_tube5.vol(:,:,slice));
roi_noise = logical(mask_noise.vol(:,:,slice));
roi_subcu = logical(mask_subcu.vol(:,:,slice));
roi_muscle = logical(mask_muscle.vol(:,:,slice));

% *************************************************
roi1 = roi_tube1; %which roi do you want to analyze?
roi1name = ['Tube 1 Slice 80'];
% *************************************************

% set te array --> NEXT TIME CREATE A LOOP TO DO THIS
% delta_TE for shorter se_mc is 7.9ms
te_short = [7.9 15.8 23.7 31.6 39.5 47.4 55.3 63.2 71.1 79 86.9 94.8 102.7 110.6 118.5 126.4 134.3 142.2 150.1 158 165.9 173.8 181.7 189.6 197.5 205.4 213.3 221.2 229.1 237 244.9 252.8];
% delta_TE for longer se_mc is 32ms
te_long = [32 64 96 128 160 192 224 256 288 320 352 384 416 448 480 512 544 576 608 640 672 704 736 768 800 832 864 896 928 960 992 1024];

%% Visualize

figure()    
imagesc(squeeze(imageShortTE(:,:,slice)), 'AlphaData', 0.8); colormap(gray); hold on; 
contour(squeeze(mask_tube1.vol(:,:,slice)), 'r','LineWidth',1); 
contour(squeeze(mask_tube2.vol(:,:,slice)), 'r','LineWidth',1);
contour(squeeze(mask_tube4.vol(:,:,slice)), 'r','LineWidth',1); 
contour(squeeze(mask_tube5.vol(:,:,slice)), 'r','LineWidth',1);
contour(squeeze(mask_tube3.vol(:,:,slice)), 'r','LineWidth',1);
contour(squeeze(mask_noise.vol(:,:,slice)), 'r','LineWidth',1);
contour(squeeze(mask_subcu.vol(:,:,slice)), 'r','LineWidth',1);
contour(squeeze(mask_muscle.vol(:,:,slice)), 'r','LineWidth',1);
title('Image Overlaid with All Masks')

%% Save Results
