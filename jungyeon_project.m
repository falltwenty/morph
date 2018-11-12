%% Clear before starting
clear all
close all

%% Load image of cell (in tif)
cd('/Users/jungyeonpark/Documents/Sophomore Fall/CS50/cells')
files = dir('*.tif');
image = imread(files(1).name);

%% Modify image
% Find Threshold for image
[~, threshold] = edge(image, 'sobel');
fudgeFactor = .5;
bw_image = edge(image, 'sobel', threshold * fudgeFactor);
figure
subplot(2, 3, 1)
imshow(bw_image)
title('Black and White Image');

% Dilate image with line
dil90 = strel('line', 3, 90);
dil0 = strel('line', 3, 0);
dil_image = imdilate(bw_image, [dil90 dil0]);
subplot(2, 3, 2)
imshow(dil_image)
title('Dilated Cells Image');

% Fill interior of holes
filled_image = imfill(dil_image, 'holes');
subplot(2, 3, 3)
imshow(filled_image)
title('Filled Disk Cells Image');

% Remove items near border
nobord_image = imclearborder(filled_image, 4);
subplot(2, 3, 4)
imshow(nobord_image)
title('Border Cleared Image');

% Smoothen cell shape
diamond = strel('diamond', 1);
final_image = imerode(nobord_image, diamond);
final_image = imerode(final_image, diamond);
subplot(2, 3, 5)
imshow(final_image)
title('Final Cells Image');

% Eliminate Small Light Dots
clean_image = bwareaopen(final_image, 300);
subplot(2, 3, 6)
imshow(clean_image)
title('Cleaner Image')

%% Show final image as outline
close all
outline = bwperim(clean_image);
outlined_image = image;
outlined_image(outline) = 255;
% figure
% imshow(outlined_image, [])
% title('Outlined Original Image');

%% Get information about cells
[labeled, n_labels] = bwlabel(outline);
cells_data = cell([n_labels, 7]);
mult_slp = zeros(1, n_labels);
el_area = zeros(1, n_labels);

for i = 1 : n_labels
    % First index is the coordinates of perimeter
    [row, col] = find(labeled == i);
    cells_data{i, 1} = [row, col]';
    perim = cells_data{i, 1};
    
    % Second index is coordinate of center
    cells_data{i, 2} = [mean(row), mean(col)]';
    centr = cells_data{i, 2};
    
    % Third index is distance from each point on perimeter to center
    dist = zeros(1, length(row));
    for j = 1 : length(perim)
        dist(1, j) = sqrt((perim(1, j) - centr(1, 1)).^2 + (perim(2, j) - centr(2, 1)).^2);
    end
    cells_data{i, 3} = dist;
    
    % Fourth index is coordinates of point closest from center
    [min_val, min_index] = min(dist);
    cells_data{i, 4} = perim(:, min_index);
    coord_clos = cells_data{i, 4};
    
    % Fifth index is coordinates of point furthest from center
    [max_val, max_index] = max(dist);
    cells_data{i, 5} = perim(:, max_index);
    coord_far = cells_data{i, 5};
    
    % Sixth index is values of area if ellipse (pi*a*b)
    cells_data{i, 6} = pi * min_val * max_val;
    el_area(1, i) = cells_data{i, 6};
    
    % Seventh index is product of slopes to see if it's perpendicular (-1)
    temp = coord_clos - centr;
    slp_clos = temp(2, 1) / temp(1, 1);
    temp = coord_far - centr;
    slp_far = temp(2, 1) / temp(1, 1);
    mult_slp(1, i) = slp_clos * slp_far;
    cells_data{i, 7} = mult_slp(1, i);
end

%% Eliminate non-oval shaped objects        
[labeled, n_labels] = bwlabel(outline);
stats = regionprops('table', clean_image, 'Centroid',...
    'MajorAxisLength', 'MinorAxisLength', 'FilledArea', 'Orientation');

el_area_2 = pi .* stats.MajorAxisLength .* stats.MinorAxisLength ./ 4;
area_diff = el_area_2 - stats.FilledArea;

new_outline = zeros(size(image));
k = 1;
for i = 1 : n_labels
    if area_diff(i) < 60
        for j = 1 : length(cells_data{i, 1})
            perim = cells_data{i, 1};
            loc = perim(:, j);
            r = loc(1, 1);
            c = loc(2, 1);
            new_outline(r, c) = 1;
        end
        good_table(k, :)= stats(i, :);
        k = k + 1;
    end
end
imshow(new_outline)

% Make outline into filled cells
new_filled_image = imfill(new_outline, 'holes');
imshow(new_filled_image)
title('Image with only proper cells')
  
% Show final image and original image
figure
imshowpair(image, new_filled_image, 'montage')
title('Original image and found cells side by side')
figure
imshowpair(image, new_filled_image, 'blend')
title('Found cells on top of original image')

%% Best not (yet) achieved: find dividing cells
[new_labeled, new_n_labels] = bwlabel(new_outline);
new_cells_data = cell(new_n_labels, 2);
coord_far = zeros(new_n_labels, 2);

for i = 1 : new_n_labels
    % Coordinates of the perimeter of new cells
    [new_row, new_col] = find(new_labeled == i);
    new_cells_data{i, 1} = [new_col, new_row]';
    new_perim = new_cells_data{i, 1}';
    
    % Coordinate of pole on major axis
    coord_far(i, 1) = good_table.Centroid(i, 1) - 0.5 * good_table.MajorAxisLength(i) * cosd(good_table.Orientation(i));
    coord_far(i, 2) = good_table.Centroid(i, 2) - 0.5 * good_table.MajorAxisLength(i) * sind(good_table.Orientation(i));
    
    % Divying up the major axis
    mton_1 = (coord_far + good_table.Centroid * 2) / 3;
    mton_2 = (mton_1 + good_table.Centroid * 2) / 3;
    
    % Finding the slopes for comparison
    sl_val = - atand(good_table.Orientation(i));
    size_gtable = size(good_table);
    length_gtable = size_gtable(1, 1);
    sl_mat = zeros(length_gtable, 1);
    for j = 1 : length(new_perim)
        sl_mat(j, 1) = (mton_1(i, 2) - new_perim(j, 2)) / (mton_1(i, 1) - new_perim(j, 1));
    end
    new_cells_data{i, 2} = sl_mat;
    
    % Find the slope most similar to sl_val
    diff_sl = sl_mat - sl_val;
    [temp, sl_simind] = min(abs(diff_sl));
    new_cells_data{i, 3} = [temp, sl_val, sl_simind, sl_mat(sl_simind)];
end
    
    
%% Sources Used
% https://www.mathworks.com/help/images/examples/detecting-a-cell-using-image-segmentation.html