files = {'plate 1_4x_350us_1.png'
    'plate 1_4x_450us_1.png'
    'plate 1_4x_450us_2.png'
    'plate 1_4x_450us_3.png'
    'plate 1_4x_450us_4.png'
    'plate 1_4x_350us_1.png'
    'plate 2_4x_340us_1.png'
    'plate 2_4x_340us_2.png'};
n_files = length(files);

area_range = [1200 15000];
ecc_range = [0 0.95];
min_dist_range = [100 900];
edge_border = 80;
se = strel('disk',12);

for i = 1:n_files
    im1 = imread(fullfile('Octopi',files{i}));
    
    % convert to 8-bit and invert
    im1_gray = im2uint8(rgb2gray(im1)); 
    im1_inv = imcomplement(im1_gray);
    
    % crop
    im1_crop = imcrop(im1_inv);
    if isempty(im1_crop)
        break
    end

    % apply median filter to remove some speckle
    im1_filt = medfilt2(im1_crop,[10 10]);
    % apply top-hat filter to remove more speckle
    im1_filt = im1_filt - imtophat(im1_crop,se); 
    % apply gaussian blur
    im1_filt = imgaussfilt(im1_filt,10);
    
    % binarize image
    im1_bw = imbinarize(im1_filt,'adaptive','Sensitivity',0.3);

    % generate image labels and measure region properties
    im1_label = bwlabel(im1_bw);
    stats = regionprops(im1_bw,'area','centroid','eccentricity');
    blob_area = [stats.Area];
    blob_ecc = [stats.Eccentricity];
    
    % keep regions within the area and eccentricity specs
    keep_area = (blob_area>area_range(1)) & (blob_area<area_range(2));
    keep_ecc = (blob_ecc>ecc_range(1)) & (blob_ecc<ecc_range(2));
    keep_idx = find(keep_area&keep_ecc);

    % remove regions close to the image boundary
    cent = reshape([stats(keep_idx).Centroid],2,length(keep_idx))';
    idx1 = any(cent<edge_border,2);
    idx2 = any([cent(:,1)>size(im1_bw,2)-edge_border, cent(:,2)>size(im1_bw,1)-edge_border],2);
    keep_idx = keep_idx(~any([idx1, idx2],2));
    
    % create masks for the accepted and rejected blobs
    reject_idx = setdiff(1:length(stats),keep_idx);
    keep_mask = ismember(im1_label,keep_idx);
    reject_mask = ismember(im1_label,reject_idx);
    
    % plot results
    fig = figure;
    subplot(2,2,1); imagesc(im1_crop);
    subplot(2,2,2); imagesc(im1_filt);
    subplot(2,2,3); imagesc(im1_bw);
    for j = 1:length(stats)
        text(stats(j).Centroid(1),stats(j).Centroid(2),sprintf('%.f',stats(j).Area),'FontSize',12);
        text(stats(j).Centroid(1),stats(j).Centroid(2)+20,sprintf('%.2f',stats(j).Eccentricity),'FontSize',12);
    end
    subplot(2,2,4); imagesc(keep_mask*3 + double(reject_mask));
    drawnow
%     waitforbuttonpress; close all
    
    img = imshowpair(im1_crop,keep_mask);
    
    % save images
    imwrite(im1_crop,[files{i}(1:end-4) '_crop.png']);
    imwrite(im1_bw,[files{i}(1:end-4) '_mask.png']);
    imwrite(img.CData,[files{i}(1:end-4) '_overlay.png']);
    
    close all
end