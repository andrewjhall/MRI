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
dir = ['/Volumes/cimalab/lcolucci/MRI/2015-07-06_COLUCCI_PHANTOMS_MEAT/Niftis/se_mc_shorter_TEs.nii']; % Path to image you want to analyze 

slice=81; % which slice # do you want to analyze? 
% ^WHY IS IT in slice 81 instead of 80 like it was created in?????????????

% set te array --> NEXT TIME DO THIS AUTOMATICALLY (PULL FROM DATA FILE)
te_short = [8.6 17.2 25.8 34.4 43 51.6 60.2 68.8 77.4 86 94.6 103.2 111.8 120.4 129 137.6 146.2 154.8 163.4 172 180.6 189.2 197.8 206.4 215 223.6 232.2 240.8 249.4 258 266.6 275.2]; 
te_long = [50:50:1600]; % 50ms increments from 50-1600ms. 

% *************************************************


%% Add Freesurfer files to path
path(path, '/Applications/freesurfer')
path(path, '/Applications/freesurfer/bin')
path(path, '/Applications/freesurfer/matlab')
path(path, '../Dialysis/functions')

%% Load  NIFTI files
dir = ['/Volumes/Transcend/COLUCCI_PHANTOMS_20150617_mgh/FLASH_3deg_30deg'];
hdr = load_nifti('/Volumes/Transcend/COLUCCI_PHANTOMS_20150617/mgh/T1.nii');

% Shorter TE's
filename_short=['t1']; %20150121_semc_shorterTE_002']; 
hdr_short = load_nifti([dir filename_short '.mgz']);
img_short_4darray = squeeze(hdr_short.vol); 
sz_4darray_short = size(img_short_4darray);
slices_short = sz_4darray_short(3);
timepts_short = sz_4darray_short(4);



%% Load Masks
mask_botleft = load_nifti('/Volumes/cimalab/lcolucci/Data/MRI/2015-01-21_PhantomsAE_MIT/Masks/se_mc_shorter-mask_botleft.nii');
mask_topleft = load_nifti('/Volumes/cimalab/lcolucci/Data/MRI/2015-01-21_PhantomsAE_MIT/Masks/se_mc_shorter-mask_topleft.nii');
mask_botmid = load_nifti('/Volumes/cimalab/lcolucci/Data/MRI/2015-01-21_PhantomsAE_MIT/Masks/se_mc_shorter-mask_botmid.nii');
mask_topright = load_nifti('/Volumes/cimalab/lcolucci/Data/MRI/2015-01-21_PhantomsAE_MIT/Masks/se_mc_shorter-mask_topright2.nii');
mask_botright = load_nifti('/Volumes/cimalab/lcolucci/Data/MRI/2015-01-21_PhantomsAE_MIT/Masks/se_mc_shorter-mask_botright.nii');
mask_noise = load_nifti('/Volumes/cimalab/lcolucci/Data/MRI/2015-01-21_PhantomsAE_MIT/Masks/se_mc_shorter-mask_noise.nii');

% *************************************************
slice=1; % which slice # do you want to analyze? 
% *************************************************

roi_botleft = logical(mask_botleft.vol(:,:,slice));
roi_topleft = logical(mask_topleft.vol(:,:,slice));
roi_botmid = logical(mask_botmid.vol(:,:,slice));
roi_topright = logical(mask_topright.vol(:,:,slice));
roi_botright = logical(mask_botright.vol(:,:,slice));
roi_noise = logical(mask_noise.vol(:,:,slice));

% *************************************************
roi1 = roi_topright; %which roi do you want to analyze?
roi1name = ['topright slice 1'];
% *************************************************

% set te array --> NEXT TIME CREATE A LOOP TO DO THIS
% delta_TE for shorter se_mc is 7.9ms
te_short = [7.9 15.8 23.7 31.6 39.5 47.4 55.3 63.2 71.1 79 86.9 94.8 102.7 110.6 118.5 126.4 134.3 142.2 150.1 158 165.9 173.8 181.7 189.6 197.5 205.4 213.3 221.2 229.1 237 244.9 252.8];
% delta_TE for longer se_mc is 32ms
te_long = [32 64 96 128 160 192 224 256 288 320 352 384 416 448 480 512 544 576 608 640 672 704 736 768 800 832 864 896 928 960 992 1024];

%% Visualize
% visualize all slices
figure;  
for n=1:slices_short
    subplot(2,4,n); 
    imagesc(squeeze(img_short_4darray(:,:,n,3))); 
    colormap(gray)
    title('Shorter TEs All Slices. 1 Timepoint')
end
 
figure;  
for n=1:slices_long
    subplot(2,4,n); 
    imagesc(squeeze(img_long_4darray(:,:,n,3))); 
    colormap(gray)
    title('Longer TEs. All Slices. 1 Timepoint')
end
 
% Overlay all masks with the images
figure
for n=1:slices_short
    subplot(2,4,n)
    imagesc(squeeze(img_short_4darray(:,:,n,1)), 'AlphaData', 0.5); hold on; contour(squeeze(mask_botleft.vol(:,:,n))); 
%     contour(squeeze(mask_botright.vol(:,:,n))); 
%     contour(squeeze(mask_topright.vol(:,:,n)));
    title('Shorter TEs. All Slices. Overlay Masks.')
end

figure
for n=1:slices_long
    subplot(2,4,n)
    imagesc(squeeze(img_long_4darray(:,:,n,1)), 'AlphaData', 0.5); hold on; contour(squeeze(mask_botleft.vol(:,:,n))); 
    contour(squeeze(mask_topleft.vol(:,:,n)));
    contour(squeeze(mask_topright.vol(:,:,n))); 
    contour(squeeze(mask_botright.vol(:,:,n)));
    contour(squeeze(mask_botmid.vol(:,:,n)));
    contour(squeeze(mask_noise.vol(:,:,n)));
    title('Longer TEs. All Slices. Overlay Masks.')
end


% visualize all time points for 1 slice & name them
figure
for n=1:timepts_short
    subplot(4,8,n)
    imagesc(squeeze(img_short_4darray(:,:,1,n)));
    colormap(gray)
    title('Shorter TEs. 1 Slice. All Timepoints')
    %img1_te{n}=img(:,:,1,n);
    % try structures
end

% % Visualize masks to ensure that they are same orientation as the image
% 
% figure; 
% for n=1:slices
%     subplot(2,4,n); 
%     imagesc(squeeze(mask_botright.vol(:,:,n))); 
%     colormap(gray)
%     title('Bot. Right - Mask through all Slices')
% end
% 
% figure; 
% for n=1:slices
%     subplot(2,4,n); 
%     imagesc(squeeze(mask_noise.vol(:,:,n))); 
%     colormap(gray)
%     title('Mask through all Slices')
% end


%% Flip masks accordingly
% I haven't figured out how to flip
% mask_subcu.vol(:,:,1) = flipup(mask_subcu.vol(:,:,1));
% mask_subcu.vol(:,:,2) = flipup(mask_subcu.vol(:,:,2));
% mask_subcu.vol(:,:,5) = flipup(mask_subcu.vol(:,:,5));
% mask_subcu.vol(:,:,6) = flipup(mask_subcu.vol(:,:,6));

%% Define images for each TE

img1_te1 = img_short_4darray(:,:,slice,1);
img1_te2 = img_short_4darray(:,:,slice,2);
img1_te3 = img_short_4darray(:,:,slice,3);
img1_te4 = img_short_4darray(:,:,slice,4);
img1_te5 = img_short_4darray(:,:,slice,5);
img1_te6 = img_short_4darray(:,:,slice,6);
img1_te7 = img_short_4darray(:,:,slice,7);
img1_te8 = img_short_4darray(:,:,slice,8);
img1_te9 = img_short_4darray(:,:,slice,9);
img1_te10 = img_short_4darray(:,:,slice,10);
img1_te11 = img_short_4darray(:,:,slice,11);
img1_te12 = img_short_4darray(:,:,slice,12);
img1_te13 = img_short_4darray(:,:,slice,13);
img1_te14 = img_short_4darray(:,:,slice,14);
img1_te15 = img_short_4darray(:,:,slice,15);
img1_te16 = img_short_4darray(:,:,slice,16);
img1_te17 = img_short_4darray(:,:,slice,17);
img1_te18 = img_short_4darray(:,:,slice,18);
img1_te19 = img_short_4darray(:,:,slice,19);
img1_te20 = img_short_4darray(:,:,slice,20);
img1_te21 = img_short_4darray(:,:,slice,21);
img1_te22 = img_short_4darray(:,:,slice,22);
img1_te23 = img_short_4darray(:,:,slice,23);
img1_te24 = img_short_4darray(:,:,slice,24);
img1_te25 = img_short_4darray(:,:,slice,25);
img1_te26 = img_short_4darray(:,:,slice,26);
img1_te27 = img_short_4darray(:,:,slice,27);
img1_te28 = img_short_4darray(:,:,slice,28);
img1_te29 = img_short_4darray(:,:,slice,29);
img1_te30 = img_short_4darray(:,:,slice,30);
img1_te31 = img_short_4darray(:,:,slice,31);
img1_te32 = img_short_4darray(:,:,slice,32);

% Images for Longer TEs
for n=1:length(te_long)
    value1 = img_long_4darray(:,:,slice,n);
    field1 = sprintf('te%g',te_long(n)); 
    images_long.(field1) = value1 % The structure 'images_long' stores the image arrays at each TE for 1 slice
end 

%% Define an ROI manually 
% define an roi
% figure; imagesc(squeeze(img1_te1));
%roi1 = roipoly; 

%% Calculate average of the signal with each ROI

% I don't think this commented out section worked 
% Calculate averages of the ROIs
% means = zeros(1,10);
% for i=1:10
%     means(i) = mean(img_slice1.te{i}(roi1)); 
% end

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
te_short_nopt1 = te_short(1,2:timepts_short);
te_long_nopt1 = te_long(1,2:timepts_long);
means_short_nopt1 = means_short(1,2:timepts_short);
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

