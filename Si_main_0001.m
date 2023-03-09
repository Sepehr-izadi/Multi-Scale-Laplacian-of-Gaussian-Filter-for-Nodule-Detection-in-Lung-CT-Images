% Step 1: Read the CT image
I = dicomread('E:\MATLAB\Project\CT-Scsn-rieh\DICOM FILE\SR_1\IM00040.dcm');

% Step 2: Preprocessing
% Convert the image to grayscale
I_gray = mat2gray(I);

% Apply median filtering to denoise the image
I_filtered = medfilt2(I_gray);

% Step 3: Nodule detection
% Use a multi-scale Laplacian of Gaussian (LoG) filter to detect potential nodules
scales = 3:3:30;
log_responses = zeros(size(I_filtered,1), size(I_filtered,2), length(scales));
for k = 1:length(scales)
    sigma = scales(k)/sqrt(2);
    filter_size = 2*ceil(3*sigma)+1;
    LoG = sigma^2*fspecial('log', filter_size, sigma);
    log_responses(:,:,k) = abs(imfilter(I_filtered, LoG, 'same', 'replicate'));
end
log_sum = sum(log_responses, 3);
threshold = graythresh(log_sum);
BW_nodules = imbinarize(log_sum, threshold);

% Use morphological operations to remove small objects and fill in gaps
se = strel('disk', 2);
BW_nodules = imopen(BW_nodules, se);
BW_nodules = imfill(BW_nodules, 'holes');
BW_nodules = bwareaopen(BW_nodules, 50);

% Step 4: Visualization
% Overlay the detected nodules on the original image
imshow(I_gray);
hold on;
B = bwboundaries(BW_nodules);
for k = 1:length(B)
    boundary = B{k};
    plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 1);
end

% Step 5: Output
% Save the output image with the detected nodules marked
imwrite(I_gray, 'lung_nodules.png');
