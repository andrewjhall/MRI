% ------------ LOAD CONFIG AND SET UP ENVIRONMENT ----------------------
clear; close all; 
workingDir = pwd; 
% Connect to folders that contain scripts I need
path(path, fullfile(workingDir, 'functions')); 

% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')
path(path, '../Dialysis/functions')
%Add dialysis functions to path
path(path, '/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/Dialysis/functions'); 

%% Load  NIFTI files
pathNiftis = '/Volumes/cimalab/lcolucci/MRI/2016-01-13_SE_SEMC_Phantoms/Niftis/'; 
pathMasks = '/Volumes/cimalab/lcolucci/MRI/2016-01-13_SE_SEMC_Phantoms/Masks';
savePath = '/Users/linacolucci/Documents/1Lab/Dialysis Study/Analysis/outputs-mri/se /2016-01-13_SE_SEMC_Phantoms'; 
sliceVis=1; timept=10; slice=1:3; 
te = [8,11, 32, 40, 51, 64, 76, 88, 102, 112,128,144,153,168,178,192,204,255,281,306,331,357,382,...
    434,484,535,586,637,688,739,790,816,840,890,940,990]; 

allSEs = dir(char(fullfile(pathNiftis, 'se_te*.nii'))); %Find all se files with .nii extension
images=[]; 
for ii = 1:size(allSEs,1) %allSEs has one row per file, so loop over those
    loadNifti = load_nifti(char(fullfile(pathNiftis, allSEs(ii).name))); %store in cell array 
    images = cat(4,images,loadNifti.vol); %image is located in 'vol' field
end 
nTimePts = length(allSEs); 
nSlices = 1:size(images, 3); 

%% Load Masks 
%      Masks are 3D arrays of 0's and 1's that have a dimension (x,y,nSlices)
%      Each row in the cell array 'masks' contains a different mask, they are
%      in the order that that they were read during the 'dir' command (see 'maskFiles')
maskFiles = dir(char(fullfile(pathMasks, '*.nii'))); %directory with names of all .nii files (structure)
maskFiles = {maskFiles.name}; %put .nii filenames in cell array
masks = cell(size(maskFiles)); %initialize 
for ii=1:length(maskFiles)
    tempLoad = load_nifti(char(fullfile(pathMasks, maskFiles{ii}))); 
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

figure; imagesc(squeeze(images(:,:,sliceVis,timept)), 'AlphaData', 0.8); colormap(gray); hold on; 
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
    for j = 1:nTimePts
        for k=slice
            tempImage = images(:,:,k,j); %one slice and one time point
            tempSum(k) = sum(sum(tempImage(tempROI(:,:,k)))); %sum of all pixel intensities in this ROI
            nPixels(k) = sum(sum(tempROI(:,:,k))); %number of pixels 

        end
        
        sumAllSlices = sum(tempSum); 
        nPxAllSlices = sum(nPixels); 
        means(j,i) = sumAllSlices/nPxAllSlices; %Signal = (Sum of All Pixel Intensities at a Particular ROI)/(Number Pixels at that ROI)
    end
end

%% Plot the T2 decays
for j=1:length(rois)
    figure; 
            plot(te, means(:,j), 'k-', 'LineWidth',2.5)
            hold on
            legend('1 Slice (short)', '1 Slice (long)', 'All Slices (short)', 'All Slices (long)')
    name = sprintf('%s_slices%s',strrep(maskFiles{j},'.nii',''),num2str(slice)); 
    title(name,'Interpreter','none')
    xlabel('TE (ms)'); ylabel('Intensity (A.U.)');
end


%% Fittings

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
result{1,1} = sprintf('This is 1-exp fitting result for data from %s.', '2016-01-13_SE_SEMC_Phantoms'); 
result{2,1} = sprintf('Scans are: %s and %s.', 'something', 'something'); 
result{3,1} = sprintf('Masks come from: %s.', 'regularMasks'); 
result{4,1} = sprintf('Images are visualizations of time point %d and slice %d.', timept, sliceVis); 
result{5,1} = sprintf(strcat('The means array comes from analyzing slices: ',num2str(slice)));
result{6,1} = sprintf('Analysis run on %s by %s', datestr(now), 'LAC'); 
columnHeader = [{'ROI'},{'Amp'}, {'95% Conf Int of A'}, {'T2 Relax Time'}, {'95% Conf Int of T2 Relax Time'} ,{'Rsquare'}, {'SSE'}, {'RMSE'}, {'SNR'}]; 
maskFiles(noiseIndex) = []; %remove noise from 'maskFiles' since there is no fitting data for 'fittings' and cell arrays must match in size for the following row. 
result(9:end,:) = vertcat(columnHeader, horzcat(maskFiles', num2cell(fittings))); %result as cell array. ROI's as rows. Fitting results as column variables. 

        cell2csv(char(fullfile(savePath, 'results')), result, ',')

