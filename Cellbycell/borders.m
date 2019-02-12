%% Collect data about cells and boundaries

clear image_borders checknumbers borders_bin B L b_valid
borders_bin = imbinarize(rgb2gray(Image_borders),0);
borders_bin(1,:) = 0;
borders_bin(end,:) = 0;
borders_bin(:,1) = 0;
borders_bin(:,end) = 0;
[B,L,N,A] = bwboundaries(borders_bin,'holes');

im_cells_data=regionprops(L,'Centroid', 'Area', 'Eccentricity','PixelList','Orientation');

%% Recording cells and their boundaries for cells larger than 200 px
cell_data = zeros(1, 4);
b_valid = cell(0);
cell_counter = 0;
for i=1:numel(im_cells_data)
    if im_cells_data(i).Area > 1000 
        if sum(A(:,i)) == 0
            cell_counter = cell_counter + 1;
            cell_data(cell_counter,1) = cell_counter;
            cell_data(cell_counter,2) = i;
            cell_data(cell_counter,3) = im_cells_data(i).Area;
            cell_data(cell_counter,4) = im_cells_data(i).Eccentricity;
            cell_data(cell_counter,5) = im_cells_data(i).Orientation;
            b_valid(cell_counter) = B(i);
        end
    end
end

cell_data(cell_data(:,5)<0,5) = cell_data(cell_data(:,5)<0,5) + 180;