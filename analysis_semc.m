%% SE-MC MRI Analysis
% 
%  This code is meant for analyzing T2 fittings from the se-mc (spin echo) scans of 
%  a Siemens scanner. It finds the means of a region of interest (ROI) in
%  an image, plots those means vs. TE (echo time array), and then performs
%  a fitting on that decay curve to find the T2 relaxation time(s). 
% 
% 
%  EDITS
%  Adapted from t2mapanalysis_SE_MC_Shorter_better - Lina A. Colucci 
%  17 Sept 2015 - Adapted from gre_analysis.m - Lina A. Colucci 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

clear; close all; 

% ************ User Inputs *************************************
dirShort = ['/Volumes/cimalab/lcolucci/MRI/20150919_COLUCCI_PHANTOMS/Niftis/se_mc_shorter_te.nii']; % Path to image you want to analyze 
dirLong = ['/Volumes/cimalab/lcolucci/MRI/20150919_COLUCCI_PHANTOMS/Niftis/se_mc_longer_te.nii']; % Path to 2nd image you want to analyze


slice=3; % which slice # do you want to analyze? 
% ^WHY IS IT in slice 81 instead of 80 like it was created in?????????????
% b/c FSLView starts on slice 0 so I need to draw masks on layer 79 if I
% want slice 80, for example. 


timept = 1; %which TE do you want to visualize? 

% set te array --> NEXT TIME DO THIS AUTOMATICALLY (PULL FROM DATA FILE)
te_short = [8.6 17.2 25.8 34.4 43 51.6 60.2 68.8 77.4 86 94.6 103.2 111.8 120.4 129 137.6 146.2 154.8 163.4 172 180.6 189.2 197.8 206.4 215 223.6 232.2 240.8 249.4 258 266.6 275.2]; 
te_long = [50:50:1600]; % 50ms increments from 50-1600ms. 

%te_short= [8600 17200 25800 34400 43000 51600 60200 68800 77400 86000 94600 103200 111800 120400 129000 137600 146200 154800 163400 172000 180600 189200 197800 206400 215000 223600 232200 240800 249400 258000 266600 275200]/1000;


% *************************************************
% roi1 = roi_tube1; %which roi do you want to analyze?
% roi1name = ['topright slice 1'];
% *************************************************
% *************************************************


%% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')
path(path, '../Dialysis/functions')

%% Load  NIFTI files

% Short TE's
loadShortTE = load_nifti(dirShort); 
imagesShortTE = loadShortTE.vol;
sizeImageShortTE = size(imagesShortTE); 
nSlicesShort = sizeImageShortTE(3);
nTimePtsShort = sizeImageShortTE(4);

% Long TE's
loadLongTE = load_nifti(dirLong); 
imagesLongTE = loadLongTE.vol; 
sizeImageLongTE = size(imagesLongTE); 
nSlicesLong = sizeImageLongTE(3); 
nTimePtsLong = sizeImageLongTE(4); 

%% Load Masks
mask_tube1 = load_nifti('/Volumes/cimalab/lcolucci/MRI/20150919_COLUCCI_PHANTOMS/Masks/tubeC_slice3.nii');
% mask_tube2 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube2.nii');
% mask_tube3 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube3.nii');
% mask_tube4 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube4.nii');
% mask_tube5 = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_tube5_furthest-from-skin.nii');
% mask_noise = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_noise.nii');
% mask_subcu = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_subcu.nii');
% mask_muscle = load_nifti('/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Masks/semc_slice1_muscle.nii');


%% Define ROIs (logicals)
roi_tube1 = logical(mask_tube1.vol(:,:,slice));
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
% NOTE: If you drew masks on an image of different resolution (even if
% exact same image), the masks won't line up right. Need to draw mask on
% an image of the same size that you will analyze. 

figure()    
imagesc(squeeze(imagesShortTE(:,:,slice,timept)), 'AlphaData', 0.8); colormap(gray); hold on; 
contour(squeeze(mask_tube1.vol(:,:,slice)), 'r','LineWidth',1); 
% contour(squeeze(mask_tube2.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_tube4.vol(:,:,slice)), 'r','LineWidth',1); 
% contour(squeeze(mask_tube5.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_tube3.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_noise.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_subcu.vol(:,:,slice)), 'r','LineWidth',1);
% contour(squeeze(mask_muscle.vol(:,:,slice)), 'r','LineWidth',1);
title('Image Overlaid with All Masks')

%% Calculate average of signal for each ROI (Raw Results Matrix)
%  Result is a matrix with nRows = # of TE time points, nColumns = # ROIs.
%  You can plot each column vs. TE to get the decay for each ROI. 

for j= 1:nROIs
    tempROI = rois; 
    disp(tempROI); 
    for i=1:nTimePtsShort
        tempImageShort = (imagesShortTE(:,:,slice,i));
        tempImageLong  = (imagesLongTE(:,:,slice,i)); 
        meansShort(i,j) = mean(tempImageShort(tempROI)); 
        meansLong(i,j) = mean(tempImageLong(tempROI)); 
        disp(i)
        disp(mean(tempImageShort(tempROI)))
        %clear tempImageTimePt % Clear temp variables
    end
    
    %clear tempROI % Clear temp variables
end
    
%% Calculate average of the signal with each ROI

% I don't think this commented out section worked 
% Calculate averages of the ROIs
% means = zeros(1,10);
% for i=1:10
%     means(i) = mean(img_slice1.te{i}(roi1)); 
% end

roi1 = roi_tube1; 

mean1 = mean(img1_te1(roi1)); 
mean2 = mean(img1_te2(roi1)); 
mean3 = mean(img1_te3(roi1)); 
mean4 = mean(img1_te4(roi1)); 
mean5 = mean(img1_te5(roi1)); 
mean6 = mean(img1_te6(roi1)); 
mean7 = mean(img1_te7(roi1)); 
mean8 = mean(img1_te8(roi1)); 
mean9 = mean(img1_te9(roi1)); 
mean10 = mean(img1_te10(roi1)); 
mean11 = mean(img1_te11(roi1)); 
mean12 = mean(img1_te12(roi1)); 
mean13 = mean(img1_te13(roi1)); 
mean14 = mean(img1_te14(roi1)); 
mean15 = mean(img1_te15(roi1)); 
mean16 = mean(img1_te16(roi1)); 
mean17 = mean(img1_te17(roi1)); 
mean18 = mean(img1_te18(roi1)); 
mean19 = mean(img1_te19(roi1)); 
mean20 = mean(img1_te20(roi1));
mean21 = mean(img1_te21(roi1)); 
mean22 = mean(img1_te22(roi1)); 
mean23 = mean(img1_te23(roi1)); 
mean24 = mean(img1_te24(roi1)); 
mean25 = mean(img1_te25(roi1)); 
mean26 = mean(img1_te26(roi1)); 
mean27 = mean(img1_te27(roi1)); 
mean28 = mean(img1_te28(roi1)); 
mean29 = mean(img1_te29(roi1)); 
mean30 = mean(img1_te30(roi1)); 
mean31 = mean(img1_te31(roi1)); 
mean32 = mean(img1_te32(roi1)); 

means_short = [mean1 mean2 mean3 mean4 mean5 mean6 mean7 mean8 mean9 mean10 mean11 mean12 mean13 mean14 mean15 mean16 mean17 mean18 mean19 mean20 mean21 mean22 mean23 mean24 mean25 mean26 mean27 mean28 mean29 mean30 mean31 mean32];

% Means for Longer TEs
means_long = zeros(1, length(te_long));
for n = 1:length(te_long)
    means_long(n) = mean(images_long.(sprintf('te%g',te_long(n)))(roi1)); % 'Means_long' is an array of the ROI intensity avg at each TE for 1 slice
end


%% Plot
% throw out 1st point
te_short_nopt1 = te_short(1,2:nTimePtsShort);
te_long_nopt1 = te_long(1,2:timepts_long);
means_short_nopt1 = means_short(1,2:nTimePtsShort);
means_long_nopt1 = means_long(1,2:timepts_long); 

figure; plot(te_short_nopt1, means_short_nopt1,'-o'); hold on; plot(te_long_nopt1, means_long_nopt1,'-s')
xlabel('TE (ms)'); ylabel('Intensity (A.U.)'); title(roi1name); legend('Shorter TEs', 'Longer TEs')

% combine data from shorter and longer TEs
te_union = [te_short_nopt1, te_long_nopt1]; 
means_union = [means_short_nopt1, means_long_nopt1];
te_union = sort(te_union); means_union = sort(means_union, 'descend'); 

% Exponential fitting of data
%[fit_data,fresult,gof,out] = exponential_fit_no_offset_2exp([te_nopt1' means_nopt1']);


%% Try doing this for each pixel --> BAD IDEA THAT TOOK ALL NIGHT TO RUN
% sz_image = size(img1_te1);
% npixels = sz_image(1)*sz_image(2); 
% 
% img1 = squeeze(img(:,:,slice,:)); %this is 1 slice with all te's concatenated
% amps = zeros(sz_image(1),sz_image(2));
% r2 = zeros(sz_image(1),sz_image(2));
% 
% 
% pixel_t2 = zeros(sz_image(1),sz_image(2));
% 
% for i=1:sz_image(1) %size of image in x-direction
%     for j=1:sz_image(2) % size of image in y-direction
%         fitting = exponential_fit_no_offset_1exp([te' squeeze(img1(i,j,1:timepts))]);
%         amps(i,j) = fitting.A1(1);
%         %r2(i,j) = gof.rsquare; 
%         pixel_t2(i,j) = fitting.tau1(1); %the single-exponential T2 fitting result 
%         clear fitting fresult gof out
%     end
% end
% 
% %[fit_data,fresult,gof,out] = exponential_fit_no_offset_1exp([te' squeeze(img1(1,1,1:timepts))]);    

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

