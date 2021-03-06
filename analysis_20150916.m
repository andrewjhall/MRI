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
workingDirectory = pwd; 

% *************************************************
dir = ['/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Niftis/m0_map_despot2.nii']; % Path to image you want to analyze 

slice=80; % which slice # do you want to analyze? 
% ^WHY IS IT in slice 81 instead of 80 like it was created in?????????????

saveFilename = sprintf('20150706-M0MapDESPOT2-Slice%d',slice); 

% % set te array --> NEXT TIME CREATE A LOOP TO DO THIS
% % delta_TE for shorter se_mc is 7.9ms
% te_short = [7.9 15.8 23.7 31.6 39.5 47.4 55.3 63.2 71.1 79 86.9 94.8 102.7 110.6 118.5 126.4 134.3 142.2 150.1 158 165.9 173.8 181.7 189.6 197.5 205.4 213.3 221.2 229.1 237 244.9 252.8];
% % delta_TE for longer se_mc is 32ms
% te_long = [32 64 96 128 160 192 224 256 288 320 352 384 416 448 480 512 544 576 608 640 672 704 736 768 800 832 864 896 928 960 992 1024];

% *************************************************

%% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')
path(path, '../Dialysis/functions')

%% Load  NIFTI files

% Shorter TE's
imageShortTE_struct = load_nifti(dir);
imageShortTE = squeeze(imageShortTE_struct.vol); 
size_short = size(imageShortTE);
nSlices_short = size_short(3);

%% Load Masks
mask_tube1 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_tube1_closest-to-skin.nii');
mask_tube2 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_tube2.nii');
mask_tube3 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_tube3.nii');
mask_tube4 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_tube4.nii');
mask_tube5 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_tube5_furthest-from-skin.nii');
mask_noise = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_noise.nii');
mask_subcu = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_subcu.nii');
mask_muscle = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/map_slice80_muscle.nii');

%% Define ROIs (logicals)
roi_tube1 = logical(mask_tube1.vol(:,:,slice));
roi_tube2 = logical(mask_tube2.vol(:,:,slice));
roi_tube3 = logical(mask_tube3.vol(:,:,slice));
roi_tube4 = logical(mask_tube4.vol(:,:,slice));
roi_tube5 = logical(mask_tube5.vol(:,:,slice));
roi_noise = logical(mask_noise.vol(:,:,slice));
roi_subcu = logical(mask_subcu.vol(:,:,slice));
roi_muscle = logical(mask_muscle.vol(:,:,slice));

%% Visualize Masks
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

%% Analyze 
imageSlice = imageShortTE(:,:,slice);
resultTube1 = mean(imageSlice(roi_tube1)); 
resultTube2 = mean(imageSlice(roi_tube2)); 
resultTube3 = mean(imageSlice(roi_tube3)); 
resultTube4 = mean(imageSlice(roi_tube4)); 
resultTube5 = mean(imageSlice(roi_tube5)); 
resultSubcu = mean(imageSlice(roi_subcu)); 
resultMuscle = mean(imageSlice(roi_muscle)); 
resultNoise = mean(imageSlice(roi_noise)); 

results = [resultTube1 resultTube2 resultTube3 resultTube4 resultTube5 resultSubcu resultMuscle resultNoise]; 
headers = {'Tube1', 'Tube2', 'Tube3', 'Tube4', 'Tube5', 'Subcu', 'Muscle', 'Noise'}; 

%% Save Results
print('../outputs-mri/m0-despot2-slice80.png', '-f1','-dpng')
csvwrite_with_headers('../outputs-mri/m0-despot2-slice80.csv',results, headers)

