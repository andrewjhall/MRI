%% SE-MC MRI Analysis
% 
%  This code is meant for analyzing T2 fittings from the se-mc (spin echo) scans of 
%  a Siemens scanner. It finds the means of a region of interest (ROI) in
%  an image, plots those means vs. TE (echo time array), and then performs
%  a fitting on that decay curve to find the T2 relaxation time(s). 
% 
%  Take a 
%
%
%  EDITS
%  Adapted from t2mapanalysis_SE_MC_Shorter_better - Lina A. Colucci 
%  17 Sept 2015 - Adapted from gre_analysis.m - Lina A. Colucci 
%  29 Sept 2015 - Pull in file paths from config.csv file  (LAC)
 

% ------------ LOAD CONFIG AND SET UP ENVIRONMENT ----------------------
clear; close all; 
workingDir = pwd; 
% Connect to folders that contain scripts I need
path(path, fullfile(workingDir, 'functions')); 
% Import Config File 
% ***!!!!*** Any time there are additional columns added to config.csv, this line must change ****
[id,date,folderName,goal,contents1,pathNiftis,pathMasks,nScans,scan001,scan002,scan003,scan004,scan005,scan006,scan007,scan008,scan010,scan011,scan012,scan013,scan014,scan015,scan016,scan017] = importConfigFile('config.csv');
%MAYBE TURN THESE INTO A STRUCTURE config.[all these variables]

%% ======================= USER INPUTS ============================== %%
% Which MRI data do you want to analyze? (Enter in the ID number)
mriID = 2; 

% Which se_mc scripts do you want to analyze? (Enter in their scan numbers, i.e. column headers on config.csv file)
whichShort = scan001;
whichLong = scan002; 

% Which masks to use? (choose the image type and slice)
pathSpecificMasks = 'semc/slice1'; 

slice=1; % which slice # do you want to analyze? 
% ^WHY IS IT in slice 81 instead of 80 like it was created in?????????????
% b/c FSLView starts on slice 0 so I need to draw masks on layer 79 if I
% want slice 80, for example. 

timept = 1; %which TE do you want to visualize? 

% set te arrays --> NEXT TIME DO THIS AUTOMATICALLY (PULL FROM DATA FILE). from info-dump
te_short = [8.6 17.2 25.8 34.4 43 51.6 60.2 68.8 77.4 86 94.6 103.2 111.8 120.4 129 137.6 146.2 154.8 163.4 172 180.6 189.2 197.8 206.4 215 223.6 232.2 240.8 249.4 258 266.6 275.2]; 
te_long = [50:50:1600]; % 50ms increments from 50-1600ms. 

%% -------------------------------------------------------------------- %%

%% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')
path(path, '../Dialysis/functions')
%Add dialysis functions to path
path(path, '/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/Dialysis/functions'); 

%% Load  NIFTI files

% Load Short TE Image and Find Its Dimensions
loadShortNifti = load_nifti(char(fullfile(pathNiftis(mriID), whichShort(mriID)))); %structure with 48 fields
imageShort = loadShortNifti.vol; %image is located in 'vol' field
sizeImageShort = size(imageShort); nSlicesShort = sizeImageShort(3);nTimePtsShort = sizeImageShort(4); %pull out dimensions

% Load Long TE Image and Find Its Dimensions
loadLongNifti = load_nifti(char(fullfile(pathNiftis(mriID), whichLong(mriID)))); %structure with 48 fields
imageLong = loadLongNifti.vol; %image is located in 'vol' field
sizeImageLong = size(imageLong); nSlicesLong = sizeImageLong(3);nTimePtsLong = sizeImageLong(4); %pull out dimensions


%% Load Masks 
%      Masks are 3D arrays of 0's and 1's that have a dimension (x,y,nSlices)
%      Each row in the cell array 'masks' contains a different mask, they are
%      in the order that that they were read during the 'dir' command (see 'maskFiles')
maskFiles = dir(char(fullfile(pathMasks(mriID), pathSpecificMasks, '*.nii'))); %directory with names of all .nii files (structure)
maskFiles = {maskFiles.name}; %put .nii filenames in cell array
masks = cell(size(maskFiles)); %initialize 
for ii=1:length(maskFiles)
    tempLoad = load_nifti(char(fullfile(pathMasks(mriID), pathSpecificMasks, maskFiles{ii}))); 
    masks{ii} = tempLoad.vol; %cell array
end


%mask_tube1 = load_nifti('/Volumes/cimalab/lcolucci/MRI/20150919_COLUCCI_PHANTOMS/Masks/tubeC_slice3.nii');
% mask_tube2 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube2.nii');
% mask_tube3 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube3.nii');
% mask_tube4 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube4.nii');
% mask_tube5 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube5_furthest-from-skin.nii');
% mask_noise = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_noise.nii');
% mask_subcu = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_subcu.nii');
% mask_muscle = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_muscle.nii');


%% Define ROIs (logicals)
%      Take masks from above and turn them into logicals just for our slice of interest
%     'rois' is a cell array where each row contains the logical array for a particular mask. 
%      The rows are in the same order that they were read during the 'dir' command (see 'maskFiles')
rois = cell(size(masks)); 
for jj=1:length(rois)
    rois{jj} = logical(masks{jj}(:,:,slice)); %cell array of logical matrices
end

%roi_tube1 = logical(mask_tube1.vol(:,:,slice));
% roi_tube2 = logical(mask_tube2.vol(:,:,slice));
% roi_tube3 = logical(mask_tube3.vol(:,:,slice));
% roi_tube4 = logical(mask_tube4.vol(:,:,slice));
% roi_tube5 = logical(mask_tube5.vol(:,:,slice));
% roi_noise = logical(mask_noise.vol(:,:,slice));
% roi_subcu = logical(mask_subcu.vol(:,:,slice));
% roi_muscle = logical(mask_muscle.vol(:,:,slice));

rois = [roi_tube1]; %, 'roi_tube2','roi_tube3', 'roi_tube4', 'roi_tube5', 'roi_noise', 'roi_subcu', 'roi_muscle'}; % !!! Need to make rois a structure
nROIs = length(rois); 

%% Visualize Masks
%  Show the slice we are analyzing with the ROIs on top
%     NOTE: If you drew masks on an image of different resolution (even if
%     exact same image), the masks won't line up right. Need to draw mask on
%     an image of the same size that you will analyze. 

figure()    
imagesc(squeeze(imageShort(:,:,slice,timept)), 'AlphaData', 0.8); colormap(gray); hold on; 
for kk=1:length(masks)
    contour(masks{kk}(:,:,slice),'r','Linewidth',1);
end
title('Image Overlaid with All Masks')

%contour(squeeze(mask_tube1.vol(:,:,slice)), 'r','LineWidth',1); 
% contour(squeeze(mask_tube2.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_tube4.vol(:,:,slice)), 'r','LineWidth',1); 
% contour(squeeze(mask_tube5.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_tube3.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_noise.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_subcu.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_muscle.vol(:,:,slice)), 'r','LineWidth',1);

%% Calculate signal average for each ROI (Raw Results Matrix)
%  Result is a matrix with nRows = # of TE time points, nColumns = # ROIs.
%  You can plot each column vs. TE to get the decay for each ROI. 

for j= 1:length(rois) %loop through each roi
    tempROI = rois{j}; 
    for i=1:nTimePtsShort %loop through all time points   
                          %NOTE: this assumes nTimePtsShort = nTimePtsLong
        tempImageShort = (imageShort(:,:,slice,i));
        tempImageLong  = (imageLong(:,:,slice,i)); 
        meansShort(i,j) = mean(tempImageShort(tempROI)); 
        meansLong(i,j) = mean(tempImageLong(tempROI)); 
        clear tempImageShort tempImageLong % Clear temp variables
    end
    clear tempROI % Clear temp variables
end
    
% throw out 1st point
meansShort = meansShort(2:end, :); 
meansLong = meansLong(2:end, :);
te_short = te_short(:, 2:end);
te_long = te_long(:, 2:end); 

%% Plot the T2 decays
for j=1:length(rois)
    figure; 
    plot(te_short, meansShort(:,j), 'k-', 'LineWidth',2.5)
    hold on
    plot(te_long, meansLong(:,j),'k-', 'LineWidth',2.5)
    title(sprintf('%s',maskFiles{j}),'Interpreter','none')
    xlabel('TE (ms)'); ylabel('Intensity (A.U.)');
end

% styles = {'bo-', 'gs-','r+-','cd-','mx-','k^-','bs:','g+:','rd:','cx:','m^:','ko:'};
% styles2 = {'b-','g-','r-','c-','m-','k-','k--','m--','c--','r--','g--','b--','b-.','g-.','r-.','c-.','m-.','k-.','k:','m:','c:','r:','g:','b:'};

%% Multiexponential Fitting

% combine data from shorter and longer TEs
te = [te_short, te_long]'; 
means = [meansShort, meansLong];
% te_union = sort(te_union); means_union = sort(means_union, 'descend'); 

% Exponential fitting of data
%[fit_data,fresult,gof,out] = exponential_fit_no_offset_2exp([te_nopt1' means_nopt1']);
teForFitting = repmat(te, 1, length(rois)); %analysis function needs a te column for each ROI. Repeat TE vector N times. 
noiseIndex = find(~cellfun(@isempty, strfind(maskFiles, 'noise'))); %use strmatch to find which mask is the noise mask file
                                                                    %use 'find' to find its index in the 'maskFiles' cell array                               
noise = std(means(:,noiseIndex));
[results1Exp, results2Exp, results3Exp] = analysisT2Decay(length(rois), teForFitting, means, 0, noise, [100 10], [1000 10 5 50], []);

for 

%% SNR calculation

% Use this to manually define ROI regions instead of loading a mask
% roi_sd = roipoly; 
% roi_mean = roipoly;

noise = std(squeeze(img1_te1(roi_noise)));
signal_mean = mean(squeeze(img1_te1(roi1)));
signal_max = max(squeeze(img1_te1(roi1)));
sprintf('These are results for %s and slice %d', roi1name, slice) 
snr_mean = signal_mean/noise
snr_max = signal_max/noise

