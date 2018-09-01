clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8_dir =[filedir, '/borders'];
tif16_dir = [filedir, '/tifs_original'];
%Folder to save summarised information
if exist([filedir, '/SummaryIntensity'],'dir') == 0
    mkdir(filedir,'/SummaryIntensity');
end
sum_dir = [filedir, '/SummaryIntensity'];

%Folder to save information about cells
if exist([filedir,'/images_analysed'],'dir') == 0
    mkdir(filedir,'/images_analysed');
end
im_dir = [filedir, '/images_analysed'];

%Folder to save information about cells
if exist([filedir,'/total_intensity'],'dir') == 0
    mkdir(filedir,'/total_intensity');
end
res_dir = [filedir, '/total_intensity'];
%dilation of borders
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);
%Reading 16-bit average intensity projection files
cd(tif16_dir);
files_tif = dir('*.tif');
average = zeros(1,12);
headers1 = {'cell', 'Eccentricity', 'Area', 'Perimeter', 'Mean Intensity', 'Total intensity'};
headers2 = {'Embryo', 'Eccentricity', 'sem', 'Area', 'sem', 'Perimeter', 'sem',...
    'Mean Intensity', 'sem', 'Total intensity', 'sem', 'Number of cells'};

for g=1:numel(files_tif)
    
    %% Open images and modify
    Cad = [num2str(g),'.tif'];
    cd(tif16_dir);
    Cad_im = imread(Cad);
    bd_dir = [tif8_dir,'/', num2str(g)];
    cd(bd_dir);
    %% Collect borders and info
    
    I=imread('handCorrection.tif');
    borders_bin = imbinarize(rgb2gray(I),0);
    borders_bin(1,:) = 0;
    borders_bin(end,:) = 0;
    borders_bin(:,1) = 0;
    borders_bin(:,end) = 0;
    [B,L,N,A] = bwboundaries(borders_bin,'holes');
    
    im_cells_data=regionprops(L,'Centroid', 'Area', 'Eccentricity');
    
    cell_data = zeros(1, 3);
    b_valid = cell(0);
    cell_counter = 0;
    for i=1:numel(im_cells_data)
        if sum(A(:,i)) == 0
            cell_counter = cell_counter + 1;
            cell_data(cell_counter,1) = im_cells_data(i).Area;
            cell_data(cell_counter,2) = im_cells_data(i).Eccentricity;
            cell_data(cell_counter,3) = i;
            b_valid(cell_counter) = B(i);
        end
    end
    image1=figure;
    imshow(imadjust(Cad_im)), title('Adjusted MTs Image');
    hold on;
    for k = 1:length(b_valid)
        clear boundary_valid
        boundary = b_valid{k};
        c = im_cells_data(cell_data(k,3)).Centroid;
        c_labels = text(c(1), c(2), sprintf('%d', k),'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'middle', 'Fontsize', 10);
        set(c_labels,'Color',[1 1 0])
        plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
    end
    
    %% Saving the image
    cd(im_dir);
    image_filename = [num2str(g),'_analysed_image.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    
    [y, x] = size(Cad_im);
    
    s_cells = zeros(1,1);
    for k = 1:length(b_valid)
        boundary = b_valid{k};
        temp = zeros(y,x);
        for l = 1:length(boundary)
            temp(boundary(l,1), boundary(l,2)) = 1;
        end
        temp = imdilate(temp, [se90I se0I]);
        cc_cells=bwconncomp(temp);
        s_cells2=regionprops(cc_cells, Cad_im, 'MeanIntensity');
        s_cells(k) = s_cells2.MeanIntensity;
    end
    result = zeros(1,6);
    for k = 1:length(b_valid)
        result(k,1) = k;
        result(k,2) = cell_data(k,2);
        result(k,3) = cell_data(k,1);
        result(k,4) = length(b_valid{k});
        result(k,5) = s_cells(k);
        result(k,6) = result(k,4) * result(k,5);
    end
    
    cd(res_dir);
    
    Otput = [num2str(g),'_total_protein.csv'];
    csvwrite_with_headers(Otput,result,headers1);
    
    
    average(g,1) = g;
    average(g,2) = mean(result(:,2));
    average(g,3) = sqrt(var(result(:,2))/length(result(:,2)));
    average(g,4) = mean(result(:,3));
    average(g,5) = sqrt(var(result(:,3))/length(result(:,3)));
    average(g,6) = mean(result(:,4));
    average(g,7) = sqrt(var(result(:,4))/length(result(:,4)));
    average(g,8) = mean(result(:,5));
    average(g,9) = sqrt(var(result(:,5))/length(result(:,5)));
    average(g,10) = mean(result(:,6));
    average(g,11) = sqrt(var(result(:,6))/length(result(:,6)));
    average(g,12) = length(result(:,6));
end

cd(sum_dir);
csvwrite_with_headers('Total_protein.csv',average,headers2);

cd(currdir);
clc
clear variables
close all
