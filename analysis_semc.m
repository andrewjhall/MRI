%% SE-MC MRI Analysis
% 
%  This code is meant for analyzing T2 fittings from the se-mc (spin echo) scans of 
%  a Siemens scanner. It finds the means of a region of interest (ROI) in
%  an image, plots those means vs. TE (echo time array), and then performs
%  a fitting on that decay curve to find the T2 relaxation time(s). 
% 
%    - Takes a long TE and a short TE sequence  
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

% Where will results be saved? 
savePath = '/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/outputs-mri/semc';
saveFilename = strcat(folderName(mriID),'.csv'); 

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

%% Define ROIs (logicals)
%      Take masks from above and turn them into logicals just for our slice of interest
%     'rois' is a cell array where each row contains the logical array for a particular mask. 
%      The rows are in the same order that they were read during the 'dir' command (see 'maskFiles')
rois = cell(size(masks)); 
for jj=1:length(rois)
    rois{jj} = logical(masks{jj}(:,:,slice)); %cell array of logical matrices
end

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
means = [meansShort; meansLong];

% Calculate noise
noiseIndex = find(~cellfun(@isempty, strfind(maskFiles, 'noise'))); %use strmatch to find which mask is the noise mask file
                                                                    %use 'find' to find its index in the 'maskFiles' cell array                               
noise = std(means(:,noiseIndex)); %Noise = std deviation of noise voxel
                                  %DOUBLE CHECK THIS!!!! Also try std deviation of the entire voxel (image[noise roi]), not just mean across different TEs
                                  %Also compare to std deviation of noise floor of T2 decays

% Perform 1-exp fitting 
fittings = []; 
for i=1:length(rois)
    fittings(i,:) = createFit(te, means(:,i), [], noise); %Columns = [Amp, 95% CI of Amp, T2 (ms), 95% CI of T2, Rsquare, SSE, RMSE, SNR]
end

%% Print Results 
% Prepare results into cell array
columnHeader = [{'ROI'},{'Amp'}, {'95% Conf Int of A'}, {'T2 Relax Time'}, {'95% Conf Int of T2 Relax Time'} ,{'Rsquare'}, {'SSE'}, {'RMSE'}, {'SNR'}]; 
description = sprintf('This is 1-exp fitting result for data from %s. Scans are: %s and %s. Masks come from: %s. Images are visualizations of time point %d.', ...
    char(folderName(mriID)), char(scan001(mriID)), char(scan002(mriID)), char(pathSpecificMasks), timept); 
result = cell(length(rois)+3, 9); 
result{1,1}=description;
result(3:end,:) = vertcat(columnHeader, horzcat(maskFiles', num2cell(fittings))); %result as cell array. ROI's as rows. Fitting results as column variables. 

%Save
cell2csv(char(fullfile(savePath, saveFilename)), result, ',')

%Other things to print to the .csv file later
% - sequence parameters: TEs, TRs, resolution, FOV, etc.
% - date that analysis was performed & by whom
% - 
