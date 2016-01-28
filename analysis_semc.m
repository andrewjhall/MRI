%% SE-MC MRI Analysis
%  This code is meant for analyzing T2 fittings from the se-mc (spin echo) scans of 
%  a Siemens scanner. It finds the means of a region of interest (ROI) in
%  an image, plots those means vs. TE (echo time array), and then performs
%  a fitting on that decay curve to find the T2 relaxation time(s). 
% 
%    - Takes a long TE and a short TE sequence 
%
%  Outputs
%  Creates folder in outputs-mri and saves plot of all T2 decays for each
%  ROI, photo of masks overlayed on slice, and csv file with fitting
%  summaries for each slice. 
%
%
%  EDITS
%  Adapted from t2mapanalysis_SE_MC_Shorter_better - Lina A. Colucci 
%  17 Sept 2015 - Adapted from gre_analysis.m - Lina A. Colucci 
%  29 Sept 2015 - Pull in file paths from config.csv file  (LAC)
%  27 Oct 2015 - char() not working. Need mat2str(cell2mat()) instead. I removed nScans from importconfig function (was giving errors). 
%  31 Oct 2015 - Added capability to analyze 1 slice or several at a time. Put slice # in the plot titles. Saved figure and .csv 
%                results with slice # in their name. 
%  28 Jan 2016 - Andrew Hall in here now


% ------------ LOAD CONFIG AND SET UP ENVIRONMENT ----------------------
clear; close all; 
workingDir = pwd; 
% Connect to folders that contain scripts I need
path(path, fullfile(workingDir, 'functions')); 
% Import Config File 
% ***!!!!*** Any time there are additional columns added to config.csv, this line must change ****
[id,date,folderName,goal,contents1,pathNiftis,pathMasks,scan001,scan002,scan003,scan004,scan005,scan006,scan007,scan008,scan010,scan011,scan012,scan013,scan014,scan015,scan016,scan017] = importConfigFile('config_mri.csv');
%MAYBE TURN THESE INTO A STRUCTURE config.[all these variables]

%% ================================================================== %%
%  ======================= USER INPUTS ==============================  %
% Which MRI data do you want to analyze? (Enter in the ID number)
mriID = 7; 

% Which masks to use? (i.e. the folder in which all the masks of interest are located)
% The script will analyze all masks within this folder 
pathSpecificMasks = ''; 

% Which se_mc scripts do you want to analyze? (Enter in their scan numbers, i.e. column headers on config.csv file)
whichShort = scan003;
whichLong = scan004; 
whichIntermed = scan003; 

% Which Analysis to perform? 
%   1 = just short/long
%   2 = just intermediate
%   NOT YET IMPLEMENTED[3 = both short/long + intermediate]
whichAnalysis = 1; 

slice=[1:3]; % which slice(s) # do you want to ANALYZE? 
         % eg. 1, 6, 1:3, etc. 
         % REMEMBER! FSLView starts on slice 0 so I need to draw masks on layer 79 in FSL if I want slice 80, for example. 

timept = 1; %which TE do you want to visualize?
sliceVis = 2; %which slice do you want to VISUALIZE? 

user = 'LAC'; %who is running this analysis? 

save = 'No'; %'Yes' to save these results

% set te arrays --> NEXT TIME DO THIS AUTOMATICALLY (PULL FROM DATA FILE). from info-dump
te_short = [8:8:256];
te_long = [25.5:25.5:816]; % 50ms increments from 50-1600ms. %25.5:25.5:816
te_intermed = [10:10:320]; %10ms incredments from 10-320ms 
%% ================================================================== %%

%% Define Connections 
% match desired mriID to position in the loaded csv file
mriID = find(id == mriID);
% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')
path(path, '../Dialysis/functions')
%Add dialysis functions to path
path(path, '/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/Dialysis/functions'); 
% Where will results be saved? 
savePath = char(fullfile('/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/outputs-mri/semc',cell2mat(folderName(mriID))));
saveFilename = sprintf('%s_1exp_Slices%s.csv', pathSpecificMasks, num2str(slice)); 
% Create folder on savePath if it doesn't already exist
if ~exist(savePath, 'dir')
    mkdir(savePath); 
end
%% Load  NIFTI files

switch whichAnalysis
    case 1 
        % Load Short TE Image and Find Its Dimensions
        loadShortNifti = load_nifti(char(fullfile(pathNiftis(mriID), whichShort(mriID)))); %structure with 48 fields
        imageShort = loadShortNifti.vol; %image is located in 'vol' field
        sizeImageShort = size(imageShort); nSlicesShort = sizeImageShort(3);nTimePtsShort = sizeImageShort(4); %pull out dimensions

        % Load Long TE Image and Find Its Dimensions
        loadLongNifti = load_nifti(char(fullfile(pathNiftis(mriID), whichLong(mriID)))); %structure with 48 fields
        imageLong = loadLongNifti.vol; %image is located in 'vol' field
        sizeImageLong = size(imageLong); nSlicesLong = sizeImageLong(3);nTimePtsLong = sizeImageLong(4); %pull out dimensions
    case 2
        % Load Intermed TE Image and Find Its Dimensions
        loadIntermedNifti = load_nifti(char(fullfile(pathNiftis(mriID), whichIntermed(mriID)))); %structure with 48 fields
        imageIntermed = loadIntermedNifti.vol; %image is located in 'vol' field
        sizeImageIntermed = size(imageIntermed); nSlicesIntermed = sizeImageIntermed(3);nTimePtsIntermed = sizeImageIntermed(4); %pull out dimensions
end

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
        rois{jj} = logical(masks{jj}); %cell array of logical matrices
end


%% Visualize Masks
%  Show the slice we are analyzing with the ROIs on top
%     NOTE: If you drew masks on an image of different resolution (even if
%     exact same image), the masks won't line up right. Need to draw mask on
%     an image of the same size that you will analyze. 

figure()
switch whichAnalysis
    case 1
        imagesc(squeeze(imageShort(:,:,sliceVis,timept)), 'AlphaData', 0.8); colormap(gray); hold on; 
    case 2 
        imagesc(squeeze(imageIntermed(:,:,sliceVis,timept)), 'AlphaData', 0.8); colormap(gray); hold on;
end
for kk=1:length(masks)
    contour(masks{kk}(:,:,sliceVis),'r','Linewidth',1);
end
title(sprintf('Slice %s Overlaid with All Masks', num2str(sliceVis)))
print(char(fullfile(savePath, sprintf('Slice%s_Masks.jpg',num2str(sliceVis)))), gcf, '-dpng')

%% Calculate signal average for each ROI (Raw Results Matrix)
%  Result is a matrix with nRows = # of TE time points, nColumns = # ROIs.
%  You can plot each column vs. TE to get the decay for each ROI. 

for i = 1:length(rois)
    tempROI = rois{i};
    for j = 1:nTimePtsShort
        for k=slice
            tempImageShort = imageShort(:,:,k,j); %one slice and one time point
            sumShort(k) = sum(sum(tempImageShort(tempROI(:,:,k)))); %sum of all pixel intensities in this ROI
            nPxShort(k) = sum(sum(tempROI(:,:,k)));
            
            tempImageLong = imageLong(:,:,k,j); %one slice and one time point
            sumLong(k) = sum(sum(tempImageLong(tempROI(:,:,k)))); %sum of all pixel intensities in this ROI
            nPxLong(k) = sum(sum(tempROI(:,:,k)));
        end
        
        sumShortAllSlices = sum(sumShort); 
        nPxShortAllSlices = sum(nPxShort); 
        meansShort(j,i) = sumShortAllSlices/nPxShortAllSlices; %Signal = (Sum of All Pixel Intensities at a Particular ROI)/(Number Pixels at that ROI)
    
        sumLongAllSlices = sum(sumLong); 
        nPxLongAllSlices = sum(nPxLong); 
        meansLong(j,i) = sumLongAllSlices/nPxLongAllSlices; %Signal = (Sum of All Pixel Intensities at a Particular ROI)/(Number Pixels at that ROI)
    
        % !!! Need to add "Intermediate TE Only" (case 2) option to this loop. Copy from code below that does this for 1 slice only. 
    end
end


for j= 1:length(rois) %loop through each roi
    tempROI = rois{j}; 
    switch whichAnalysis
        case 1
            for i=1:nTimePtsShort %loop through all time points   
                                  %NOTE: this assumes nTimePtsShort = nTimePtsLong
                tempImageShort = (imageShort(:,:,sliceVis,i)); 
                tempImageLong  = (imageLong(:,:,sliceVis,i)); 
                meansShort1Slice(i,j) = mean(tempImageShort(tempROI(:,:,sliceVis))); % i=nTimePts  j=nROIs  k=nSlices
                meansLong1Slice(i,j) = mean(tempImageLong(tempROI(:,:,sliceVis))); 
                %clear tempImageShort tempImageLong % Clear temp variables
            end
        case 2 
            for i=1:nTimePtsIntermed %loop through all time pts
                tempImageIntermed = (imageIntermed(:,:,sliceVis,i)); 
                meansIntermed(i,j) = mean(tempImageIntermed(tempROI)); 
                clear tempImageIntermed %clear temp variable
            end
    end
    %clear tempROI % Clear temp variables
end


% throw out 1st point
switch whichAnalysis
    case 1
        meansShort1Slice = meansShort1Slice(2:end, :); 
        meansLong1Slice = meansLong1Slice(2:end, :);
        te_short = te_short(:, 2:end);
        te_long = te_long(:, 2:end); 
        meansShort = meansShort(2:end,:);
        meansLong = meansLong(2:end, :); 
    case 2
        meansIntermed = meansIntermed(2:end, :); 
        te_intermed = te_intermed(:,2:end); 
end

%% Plot the T2 decays
for j=1:length(rois)
    figure; 
    switch whichAnalysis
        case 1
            plot(te_short, meansShort1Slice(:,j), 'k-', 'LineWidth',2.5)
            hold on
            plot(te_long, meansLong1Slice(:,j),'k-', 'LineWidth',2.5)
            hold on
            plot(te_short, meansShort(:,j), 'r:', 'LineWidth', 2.5)
            plot(te_long, meansLong(:,j), 'r:','LineWidth',2.5)
            legend('1 Slice (short)', '1 Slice (long)', 'All Slices (short)', 'All Slices (long)')
        case 2
            plot(te_intermed, meansIntermed(:,j), 'k-', 'LineWidth', 2.5)
    end
    name = sprintf('%s_slices%s',strrep(maskFiles{j},'.nii',''),num2str(slice)); 
    title(name,'Interpreter','none')
    xlabel('TE (ms)'); ylabel('Intensity (A.U.)');
end

% styles = {'bo-', 'gs-','r+-','cd-','mx-','k^-','bs:','g+:','rd:','cx:','m^:','ko:'};
% styles2 = {'b-','g-','r-','c-','m-','k-','k--','m--','c--','r--','g--','b--','b-.','g-.','r-.','c-.','m-.','k-.','k:','m:','c:','r:','g:','b:'};

%% Multiexponential Fitting

% combine data from shorter and longer TEs
switch whichAnalysis 
    case 1
        te = [te_short, te_long]'; 
        means = [meansShort; meansLong]; %Only do following analysis for the "All Slices" means data. If I wanted to do the single slice, change to "meansShort1Slice"
    case 2
        te = te_intermed'; 
        means = meansIntermed; 
end

% Calculate noise
noiseIndex = find(~cellfun(@isempty, strfind(maskFiles, 'noise'))); %use strmatch to find which mask is the noise mask file
                                                                    %use 'find' to find its index in the 'maskFiles' cell array                               
noise = std(means(:,noiseIndex)); %Noise = std deviation of noise voxel
                                  %DOUBLE CHECK THIS!!!! Also try std deviation of the entire voxel (image[noise roi]), not just mean across different TEs
                                  %Also compare to std deviation of noise floor of T2 decays

% Perform 1-exp fitting 
fittings = []; 
for i=1:length(rois)
    if i == noiseIndex
        continue
    end
    [fittings(i,:), h] = createFit(te, means(:,i), [800 100 0], noise); %Columns = [Amp, 95% CI of Amp, T2 (ms), 95% CI of T2, Rsquare, SSE, RMSE, SNR]
    name = sprintf('%s_slice%s',strrep(maskFiles{i},'.nii',''),num2str(slice)); 
    title(sprintf('%s',name),'Interpreter','none')
    print(h, '-dpng', fullfile(savePath, strcat(name,'.jpeg')))
end

%% Print Results 
% Prepare results into cell array
result = cell(length(rois)+8, 9); 
result{1,1} = sprintf('This is 1-exp fitting result for data from %s.', mat2str(cell2mat((folderName(mriID))))); 
switch whichAnalysis
    case 1
        result{2,1} = sprintf('Scans are: %s and %s.', char(whichShort(mriID)), char(whichLong(mriID))); 
    case 2
        result{2,1} = sprintf('Scans are: %s.', char(whichIntermed(mriID))); 
end
result{3,1} = sprintf('Masks come from: %s.', char(pathSpecificMasks)); 
result{4,1} = sprintf('Images are visualizations of time point %d and slice %d.', timept, sliceVis); 
result{5,1} = sprintf(strcat('The means array comes from analyzing slices: ',num2str(slice)));
result{6,1} = sprintf('Analysis run on %s by %s', datestr(now), user); 
columnHeader = [{'ROI'},{'Amp'}, {'95% Conf Int of A'}, {'T2 Relax Time'}, {'95% Conf Int of T2 Relax Time'} ,{'Rsquare'}, {'SSE'}, {'RMSE'}, {'SNR'}]; 
maskFiles(noiseIndex) = []; %remove noise from 'maskFiles' since there is no fitting data for 'fittings' and cell arrays must match in size for the following row. 
result(9:end,:) = vertcat(columnHeader, horzcat(maskFiles', num2cell(fittings))); %result as cell array. ROI's as rows. Fitting results as column variables. 

%Save
switch save
    case 'Yes'
        cell2csv(char(fullfile(savePath, saveFilename)), result, ',')
end

%Other things to print to the .csv file later
% - sequence parameters: TEs, TRs, resolution, FOV, etc.
% - date that analysis was performed & by whom
% - 

%% IMPORT THE SE MEANS AND TE VARIABLES AND PLOT SEMC VS SE


for i=1:length(rois)
    if i == noiseIndex
        continue
    end
    figure; 
    plot(te, means(:,i),'bo'); hold on; plot(teSE, meansSE(:,i), 'r*')
    name = sprintf('%s',strrep(maskFiles{i},'.nii','')); 
    title(sprintf('%s',name),'Interpreter','none')    
    legend('SEMC', 'SE')
%     [fittings(i,:), h] = createFit(te, means(:,i), [800 100 0], noise); %Columns = [Amp, 95% CI of Amp, T2 (ms), 95% CI of T2, Rsquare, SSE, RMSE, SNR]
%     name = sprintf('%s_slice%s',strrep(maskFiles{i},'.nii',''),num2str(slice)); 
%     title(sprintf('%s',name),'Interpreter','none')
    print('-dpng', fullfile(savePath, strcat('SEvsSEMC',name,'.jpeg')))
end


figure; plot(te, means(:,1)); hold on; plot(teSE, meansSE(:,1), 'r*')

