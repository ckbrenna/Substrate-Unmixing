%% HELPFUL CODE 
% Input into Command Window

% Normalized, pre-unmixed images (replace # with image number)      figure; imagesc(luc_normalized{#})
% Images of ROIs (column is each input image, row is each ROI)      figure; imagesc(Preunmix_ROI)
% Images of ROIs after unmixing                                     figure; imagesc(Unmix_ROI)               
% Reference Matrix                                                  figure; imagesc(A)


%% INPUTS 
% Upload images and assign ROIs to reference populations 

clear all 

% CHANGE BASED ON NUMBER OF SUBSTRATES (add ",imread('[filename].TIF')" for
% each additional substrate. List images in order of substrate addition.
luc_images = {imread('4-BrLuc.TIF'),imread('D-Luc.TIF'),imread('AkaLumine.TIF')};

% CHANGE THIS BASED ON NUMBER OF SUBSTRATES AND ROI COORDINATES 
% define coordinates of ROI for each luciferase
% X = column coordinate of left pixel, Y = row coordinate of top pixel

luc1_X = 30; luc1_Y = 80;
luc2_X = 70; luc2_Y = 80;
luc3_X = 110; luc3_Y = 80;

% define dimensions of the ROI 
ROI_width = 50;
ROI_height = 15;

% CHANGE THIS BASED ON NUMBER OF SUBSTRATES (add " luc#_X luc#_Y" for each
% additional substrate. Variable name must match the name of the
% coordinate. 
ROI_coord = [luc1_X luc1_Y; luc2_X luc2_Y; luc3_X luc3_Y]; % array containing X & Y minimum coord. 


%% IMAGE NORMALIZATION
% Transform images into vector and normalize pixel values from 0-65536

luc_normalized = {};
luc_normalized_vec = {};

for i = 1:length(luc_images)
    
    % filter & normalize pixel values of image vector 
    working_image = luc_images{i};
    filter = medfilt2(working_image, [5 5]);
    luc_vec = filter(:);
    bkgsub = im2double(luc_vec - 1362); % subtract minimum pixel value (average is 1362)
    maximum_pixel = max(bkgsub); % find maximum pixel value after background subtraction
    norm_factor = im2double(65536/maximum_pixel);
    normalized_image = bkgsub * norm_factor; % multiply subtracted image by normalization factor 
    
    % store normalized image vector in a new cell array
    luc_normalized_vec{i} = normalized_image;
    
    % reshape normalized image vector into 256x256 image and store 
    luc_normalized{i} = reshape(normalized_image,[size(luc_images{i})]);
    
end


%% GENERATE MATRICES FOR UNMIXING
% I = Intensity matrix, A = Reference Matrix 

% create intensity matrix , I (pixel x substrate)
I = double(cell2mat(luc_normalized_vec)); 

% create array of all ROI coordinates
ROI = zeros(length(luc_normalized),4);
for i = 1:length(luc_normalized)
    ROI(i,:) = [ROI_coord(i,2), ROI_coord(i,2) + ROI_height, ROI_coord(i,1), ROI_coord(i,1) + ROI_width];
end

%create images to verify ROIs

ROI_images_cellarray = {};
for i = 1:length(luc_images)
    working_image = luc_images{i};
    for j = 1:size(ROI,1)
        ROI_crop = imcrop(luc_normalized{i}, [ROI_coord(j,1) ROI_coord(j,2) ROI_width ROI_height]);
        ROI_images_cellarray{j,i} = ROI_crop;
    end
end
Preunmix_ROI = cell2mat(ROI_images_cellarray);

% create positive control matrix, A (luciferase x substrate). Contains
% average pixel value in a 5x5 ROI. 
A = zeros(length(luc_normalized),length(luc_normalized));
for i = 1:length(luc_normalized) % number of images 
    working_image = luc_normalized{i};
    for j = 1:size(ROI,1) % number of ROIs 
        ROI_matrix = working_image([ROI(j,1):ROI(j,2)],[ROI(j,3):ROI(j,4)]);
        ROI_vector = ROI_matrix(:);
        ROI_max = max(ROI_vector);
        A(j,i) = ROI_max;
    end
end


%% UNMIX 
% Generate C (Contribution matrix) 

C=zeros(size(I));
for i=1:size(C,1)  
    C(i,:)=I(i,:)*inv(A);
end

% reshape C into unmixed images and collect in cell array
unmix = {};
for i = 1:size(C,2)
    unmix_matrix = reshape(C(:,i),[size(luc_images{i})]);
    unmix{i} = unmix_matrix;
    figure; imagesc(unmix{i}); colormap(gray)
end

%create images to verify Unmixing at ROIs
Unmix_images_cellarray = {};
for i = 1:length(luc_images)
    working_image = luc_images{i};
    for j = 1:size(ROI,1)
        ROI_crop = imcrop(unmix{i}, [ROI_coord(j,1) ROI_coord(j,2) ROI_width ROI_height]);
        Unmix_images_cellarray{j,i} = ROI_crop;
    end
end
Unmix_ROI = cell2mat(Unmix_images_cellarray);


%% GENERATE OUTPUTS

% Normalized pre-unmixed images 
% ADD/REMOVE LINES BASED ON NUMBER OF IMAGES
% write to .txt file for ImageJ 
dlmwrite('Pre-Unmix 1.txt',luc_normalized{1},'delimiter','\t');
dlmwrite('Pre-Unmix 2.txt',luc_normalized{2},'delimiter','\t');
dlmwrite('Pre-Unmix 3.txt',luc_normalized{3},'delimiter','\t');

% Unmixed images (pixel values are 0-1)
% ADD/REMOVE LINES BASED ON NUMBER OF IMAGES
% write to .txt file for ImageJ 
dlmwrite('Luciferin 1.txt',unmix{1},'delimiter','\t');
dlmwrite('Luciferin 2.txt',unmix{2},'delimiter','\t');
dlmwrite('Luciferin 3.txt',unmix{3},'delimiter','\t');

% Conditional value
Cond = cond(A)