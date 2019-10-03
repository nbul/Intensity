name = [b_dir, '/', num2str(k)];
cd(name);

I = imread('handCorrection.tif');
I = imbinarize(rgb2gray(I));
[im_x, im_y] = size(I);

%% data about borders
I(1,:) = 0;
I(end,:) = 0;
I(:,1) = 0;
I(:,end) = 0;
[B,L,N,A] = bwboundaries(I,'holes');

im_cells_data=regionprops(L,'Centroid', 'Area', 'Eccentricity','Orientation');

%% Recording cells and their boundaries for cells larger than 200 px
cell_data = zeros(1, 4);
b_valid = cell(0);
cell_counter = 0;
for i=1:numel(im_cells_data)
        if sum(A(:,i)) == 0 && isempty(find(A(i,:), 1)) == 0 
            cell_counter = cell_counter + 1;
            cell_data(cell_counter,1) = cell_counter;
            cell_data(cell_counter,2) = i;
            cell_data(cell_counter,3) = im_cells_data(i).Area;
            cell_data(cell_counter,4) = im_cells_data(i).Eccentricity;
            cell_data(cell_counter,5) = im_cells_data(i).Orientation;
            b_valid(cell_counter) = B(i);
        end
end