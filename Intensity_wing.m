%% Clear all and initial parameters
clc
clear all
close all

cd('../');
mkdir('data/Cells');
mkdir('data/Distribution no BG');
mkdir('data/Distribution with BG');
mkdir('data/Distribution with BG2');
mkdir('data/Intensity');
mkdir('data/Summary');

cd('data/tifs_original');
files_tif = dir('*.tif');
cd('../../');
cd('embryo-intensity/');

Cells_total = zeros(1,6); % Information about all cells from all images
Cells_number = zeros(1,2); % Number of cells in all images
n_all = 0; % Counter for all cells

Intensity_total = zeros(1,4); % Intensity of all individual borders in all wings
Intensity_average = zeros(1,6);
n_0_90 = 0; % Number of all individual borders in all wings


%% Analysis of individual images
for g=1:numel(files_tif)
    
    %% Resetting counters
    cell_orientation = 0; %Final average orientation of all cells
    BG = 0; %Average intensity of background
    BG2=0;
    q=0; % First complete cell adjacent to the border
    q2=0; %Second complete cell adjacent to the border
    celldata = zeros(2,7); %Array with data of individual cells
    cell_temp = zeros(2,8); % Temporary cell data file
    intensity = zeros(2,2); %Array with data of individual borders
    
    
    %% Open images and modify
    Cad = [num2str(g),'.tif'];
    cd('../');
    cd('data/tifs_original');
    Cad_im=imread(Cad);
    Cad_im2=uint16(Cad_im);
    cd('../../');
    cd('embryo-intensity/');
    
    cd('../');
    Folder = ['data/tifs_8bit/', num2str(g)];
    cd(Folder);
    V=imread('vertices.png');
    V=im2bw(V,1/255);
    se90 = strel('line', 10, 90);
    se0 = strel('line', 9, 0);
    % Vdil - dilated image of all vertices
    Vdil = imdilate(V, [se90 se0]);
    % Individual vertices as objects
    cc_Vdil = bwconncomp(Vdil);
    L_Vdil = labelmatrix(cc_Vdil);
    
    I=imread('tracked_bd.png');
    I2=im2bw(I,1/255);
    se90 = strel('line', 5, 90);
    se0 = strel('line', 5, 0);
    cd('../../../');
    cd('embryo-intensity/');
    
    % I_cells - inverted image of all cells that are completely in frame;
    % s_cells - individual cells as objects
    I_cells = imdilate(I2, [se90 se0]);
    I_cells=imcomplement(I_cells);
    I_cells=imclearborder(I_cells);
    cc_cells=bwconncomp(I_cells);
    L_cells=labelmatrix(cc_cells);
    s_cells=regionprops(cc_cells, Cad_im, 'MeanIntensity', 'Orientation', 'PixelIdxList', 'Eccentricity', 'Area');
    
    % I_borders - all cell-cell borders of the cells that are completely in
    % frame; Junction_cad - individual borders as objects
    I_borders = imdilate(I2, [se90 se0]);
    I_borders = imsubtract(I_borders, Vdil);
    I_borders = im2bw(I_borders, 1/255);
    cc_borders=bwconncomp(I_borders);
    Junction_cad=regionprops(cc_borders,Cad_im,'MeanIntensity','Orientation','Perimeter', 'centroid', 'Extrema', 'PixelList');
    
%     BWoutlineall = bwperim(I_borders);
%     Segoutall = ceil(Cad_im2/16);
%     Segoutall(BWoutlineall) = 255;
%     figure, imshow(Segoutall, [0 255]), title('segmented image')
%     
    %% Image data: BG and orientation
    % Background and image orientation: measuring summarized background and
    % summarized orietnations of positively and negatively oriented cells
    for i=1:numel(s_cells);
        BG2=BG2+s_cells(i).MeanIntensity;
    end
    
    BG2=BG2/numel(s_cells);
    BG = graythresh(Cad_im2)*256*16;
    
    %% Cell data
    k=0;
    for n=1:numel(s_cells)
        pixel_value = 0;
        newImg = im2bw(I,255/255);
        newImg( vertcat( s_cells(n).PixelIdxList ) ) = 1; % getting cell n only
        se90 = strel('line', 4, 90);
        se0 = strel('line', 4, 0);
        newImg = imdilate(newImg, [se90 se0]);
        cc_newImg=bwconncomp(newImg);
        s_newImg=regionprops(cc_newImg, newImg, 'Area', 'PixelList');
        for i=1:s_newImg(1).Area
            pixel_value(i) = L_Vdil(s_newImg(1).PixelList(i,2),s_newImg(1).PixelList(i,1));
        end
        pixel_value = unique(pixel_value);
        pixel_value(pixel_value==0) = [];
        % Writing data for cell n in celldata matrix
        k=k+1;
        celldata(k,1)=n;
        celldata(k,2)=s_cells(n).Area;
        celldata(k,3)=s_cells(n).MeanIntensity;
        celldata(k,4)=s_cells(n).Orientation;
        celldata(k,5) = s_cells(n).Eccentricity;
        celldata(k,6) = 1/sqrt(1-celldata(k,5)^2);
        celldata(k,7)= length(pixel_value);
    end
    
    %% Summary of cells in the images
    cell_temp = celldata;
    Vertices3 = celldata(:,7);
    k=0;
    for n=1:length(Vertices3)
        if Vertices3(n)>8
            cell_temp(n-k,:) = [];
            k=k+1;
        end;
        
        if Vertices3(n)<3
            cell_temp(n-k,:) = [];
            k=k+1;
        end;
    end
    Cells_number(g,1) = g;
    Cells_number(g,2) = length(cell_temp);
    for n=1:length(cell_temp)
        n_all = n_all + 1;
        Cells_total(n_all,1) = g;
        Cells_total(n_all,2) = cell_temp(n,1);
        Cells_total(n_all,3) = cell_temp(n,2);
        Cells_total(n_all,4) = cell_temp(n,3);
        Cells_total(n_all,5) = cell_temp(n,4);
        Cells_total(n_all,6) = cell_temp(n,5);
        Cells_total(n_all,7) = 1/sqrt(1-cell_temp(n,5)^2);
        Cells_total(n_all,8) = cell_temp(n,7);
    end
    
    %% Borders data
    k=0;
    for n=1:numel(Junction_cad)
        % Finding borders between two cells that are completely in frame
        boundary = bwtraceboundary(I_borders, [Junction_cad(n).Extrema(3,2) + 0.5  Junction_cad(n).Extrema(3,1) - 0.5], 'NW');
        boundary_value = 0;
        Size_im = size(Cad_im);
        for i=1:size(boundary*4)
            if boundary(i,1)-1>0 && boundary(i,2)-1>0
                boundary_value(i*4-3) = L_cells(boundary(i,1)-1,boundary(i,2)-1);
            else boundary_value(i*4-3) = 0;
            end
            if boundary(i,1)-1>0 && boundary(i,2)<Size_im(2)
                boundary_value(i*4-2) = L_cells(boundary(i,1)-1,boundary(i,2)+1);
            else boundary_value(i*4-2) = 0;
            end
            if boundary(i,1)+1<Size_im(1)+1 && boundary(i,2)-1>0
                boundary_value(i*4-1) = L_cells(boundary(i,1)+1,boundary(i,2)-1);
            else boundary_value(i*4-2) = 0;
            end
            if boundary(i,1)+1<Size_im(1)+1 && boundary(i,2)+1<Size_im(2)+1
                boundary_value(i*4) = L_cells(boundary(i,1)+1,boundary(i,2)+1);
            else boundary_value(i*4-2) = 0;
            end
        end
        boundary_value = unique(boundary_value);
        boundary_value(boundary_value==0) = [];
        if length(boundary_value)==2
            q = boundary_value(1,1);
            q2 = boundary_value(1,2);
            if Junction_cad(n).Perimeter/2>5 
                if Junction_cad(n).MeanIntensity>BG2
                    % Writing down parametes of individual junctions
                    % that are longer than 5px and are surrounded by
                    % cells with area>1000px and are completely within
                    % the image
                    k=k+1;
                    intensity(k,1) = Junction_cad(n).Orientation;
                    intensity(k,2) = Junction_cad(n).MeanIntensity-BG;
                    intensity(k,3) = Junction_cad(n).MeanIntensity-BG2;
                    intensity(k,4) = Junction_cad(n).MeanIntensity;
                    intensity(k,5) = Junction_cad(n).Perimeter/2;
                    %Summarised data for all wings
                    n_0_90 = n_0_90 + 1;
                    Intensity_total(n_0_90,1) = intensity(k,1);
                    Intensity_total(n_0_90,2) = intensity(k,2);
                    Intensity_total(n_0_90,3) = intensity(k,3);
                    Intensity_total(n_0_90,4) = intensity(k,4);
                    Intensity_total(n_0_90,5) = intensity(k,5);
                end
                
            end
        end
    end
    % Sorting borders by angle
    intensity=sortrows(intensity,1);
    
    %% Summary of cell borders
     Intensity_average(g,1) = g;
     Intensity_average(g,2) =  sum(intensity(:,2))/length(intensity(:,2));
     Intensity_average(g,3) =  sum(intensity(:,3))/length(intensity(:,2));
     Intensity_average(g,4) =  sum(intensity(:,4))/length(intensity(:,2));
     Intensity_average(g,5) = BG;
     Intensity_average(g,6) = BG2;
    
    %% Fitting and plotting
    equation1=polyfit(intensity(:,1),intensity(:,2),1); % With background subtraction
    equation2=polyfit(intensity(:,1),intensity(:,3),1); % With background subtraction BG2
    equation3=polyfit(intensity(:,1),intensity(:,4),1); % Without background subtraction
    x = 0:0.5:90;
    y1 = equation1(1,1)*x + equation1(1,2);
    y2 = equation2(1,1)*x + equation2(1,2);
    y3 = equation3(1,1)*x + equation3(1,2);
    
    % Plot image of fitted intensities after background subtraction
    image1=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,2), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle after background subtraction', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y1 , '-b', 'LineWidth',3);
    precision = 2;
    text_control = ['wing', num2str(g), ': ', num2str(equation1(1,2), precision) ' + ' num2str(equation1(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([0 100 0 1000]);
    
    % Plot image of fitted intensities after background subtraction BG2
    image2=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,3), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle after background subtraction BG2', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y2 , '-b', 'LineWidth',3);
    precision = 2;
    text_control = ['wing', num2str(g), ': ', num2str(equation2(1,2), precision) ' + ' num2str(equation2(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([0 100 0 1000]);
    
    % Plot image of fitted intensities without background subtraction
    image3=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,4), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle without background subtraction', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y3 , '-b', 'LineWidth',3);
    precision = 2;
    text_control = ['wing', num2str(g), ': ', num2str(equation3(1,2), precision) ' + ' num2str(equation3(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([0 100 0 1000]);
    
    %% Writing data for individual images
    % Borders data: N_distribution_with_BG - after BG subtraction;
    cd('../');
    cd('data/Distribution with BG');
    Otput_Graph = [num2str(g),'_distribution_with_BG.tif'];
    hold off
    print(image1, '-dtiff', '-r300', Otput_Graph);
    cd('../../');
    cd('embryo-intensity/');
    
    % Borders data: N_distribution_with_BG2 - after BG subtraction, BG2;
    cd('../');
    cd('data/Distribution with BG2');
    Otput_Graph = [num2str(g),'_distribution_with_BG2.tif'];
    hold off
    print(image2, '-dtiff', '-r300', Otput_Graph);
    cd('../../');
    cd('embryo-intensity/');
    
    % N_distribution_no_BG - without BG subtraction;
    cd('../');
    cd('data/Distribution no BG');
    Otput_Graph = [num2str(g),'_distribution_no_BG.tif'];
    hold off
    print(image3, '-dtiff', '-r300', Otput_Graph);
    cd('../../');
    cd('embryo-intensity/');
    
    % N_intensity - information about individual borders;
    cd('../');
    cd('data/Intensity');
    Otput_Intensity = [num2str(g),'_intensity.csv'];
    headers = {'Angle', 'Intensity-BG', 'Intensity-BG2', 'Intensity', 'Length'};
    csvwrite_with_headers(Otput_Intensity,intensity, headers);
    close all;
    cd('../../');
    cd('embryo-intensity/');
    
    % N_area - information about individual cells
    cd('../');
    cd('data/Cells');
    Otput_Celldata = [num2str(g),'_area.csv'];
    headers = {'Cell', 'Area', 'Intensity', 'Orientation', 'Eccentricity', 'Elongation', 'Neighbours'};
    csvwrite_with_headers(Otput_Celldata,celldata,headers);
    cd('../../');
    cd('embryo-intensity/');
    
end


%% Plotting summarized data
%Fitting summarized data with lines
equation1=polyfit(Intensity_total(:,1),Intensity_total(:,2),1); %with background subtraction
equation2=polyfit(Intensity_total(:,1),Intensity_total(:,3),1); % with background subtraction, BG2
equation3=polyfit(Intensity_total(:,1),Intensity_total(:,4),1); % without background subtraction
x = 0:0.5:90;
y1 = equation1(1,1)*x + equation1(1,2);
y2 = equation2(1,1)*x + equation2(1,2);
y3 = equation3(1,1)*x + equation3(1,2);

% Plot image of fitted intensities after background subtraction
image4=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,2), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle with background subtraction', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y1 , '-b', 'LineWidth',3);
precision = 2;
text_control = ['all wings: ', num2str(equation1(1,2), precision) ' + ' num2str(equation1(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([0 100 0 1000]);

% Plot image of fitted intensities without background subtraction
image5=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,3), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle with background subtraction, BG2', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y2 , '-b', 'LineWidth',3);
precision = 2;
text_control = ['all wings: ', num2str(equation2(1,2), precision) ' + ' num2str(equation2(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([0 100 0 1000]);

% Plot image of fitted intensities without background subtraction
image6=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,4), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle without background subtraction', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y3 , '-b', 'LineWidth',3);
precision = 2;
text_control = ['all wings: ', num2str(equation3(1,2), precision) ' + ' num2str(equation3(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([0 100 0 1000]);


%% Writing combined data
%Properties of individual cells
cd('../');
cd('data/Summary');
Otput_All_Celldata = ['vertices_all.csv'];
headers = {'Wing', 'Cell', 'Area', 'Intensity', 'Orientation', 'Eccentricity', 'Elongation', 'Neighbours'};
csvwrite_with_headers(Otput_All_Celldata,Cells_total, headers);

%Number of cells in each wing
Otput_cells_number = ['cell_number.csv'];
headers2 = {'Disc', 'Cells'};
csvwrite_with_headers(Otput_cells_number,Cells_number, headers2);


%Combined all borders together
Intensity_total=sortrows(Intensity_total,1);
headers = {'Angle', 'Intensity-BG', 'Intensity-BG2', 'Intensity', 'Length'};
csvwrite_with_headers('Intensity_total.csv',Intensity_total, headers);

% Distribution_all_with_BG - after BG subtraction;
Otput_Graph = ['Distribution_all_with_BG.tif'];
hold off
print(image4, '-dtiff', '-r300', Otput_Graph);

% Distribution_all_with_BG2 - after BG subtraction;
Otput_Graph = ['Distribution_all_with_BG2.tif'];
hold off
print(image5, '-dtiff', '-r300', Otput_Graph);

% Distribution_all_no_BG - without BG subtraction;
Otput_Graph = ['Distribution_all_no_BG.tif'];
hold off
print(image6, '-dtiff', '-r300', Otput_Graph);

%Averaged intensity by angle in individual embryos
headers = {'Wing', 'Intensity-BG', 'Intensity-BG2', 'Intensity', 'BG', 'BG2'};
csvwrite_with_headers('Intensity_wing.csv',Intensity_average, headers);
cd('../../');
cd('embryo-intensity/');


close all;