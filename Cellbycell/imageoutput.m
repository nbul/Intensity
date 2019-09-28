%% Drowing overlay of selected cell borders on image of cytoskeleton
image1=figure;
I2 = imdilate(I, [se90I se0I]);
imshow(I2);
hold on;
for l = 1:length(b_valid)
    clear boundary_valid
    boundary = b_valid{l};
    c = im_cells_data(cell_data(l,2)).Centroid;
    c_labels = text(c(1), c(2), sprintf('%d', l),'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'middle', 'Fontsize', 10);
    set(c_labels,'Color',[1 1 0])
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
end

%% Saving the image
cd(cells_dir);
image_filename = [num2str(k),'_analysed_image.tif'];
print(image1, '-dtiff', '-r150', image_filename);
close all


