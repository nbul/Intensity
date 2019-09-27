%% Clear all and initial parameters
clc
clear variables
close all


smallA = 300;
largeA = 25000;

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8 =[filedir, '/bordersnoOE'];
tifAB = [filedir, '/AB'];
tifGFP = [filedir, '/GFP'];
tifCad = [filedir, '/Cad'];


%Folder to save summarised information
if exist([filedir, '/SummaryRatio'],'dir') == 0
    mkdir(filedir,'/SummaryRatio');
end
sum_dir = [filedir, '/SummaryRatio'];

%Reading 16-bit average intensity projection files
cd(tifAB);
files_tif = dir('*.tif');

%% Assigning memory
Cells_total = zeros(1,6); % Information about all cells from all images
Cells_number = zeros(1,2); % Number of cells in all images
n_all = 0; % Counter for all cells

Intensity_total = zeros(1,4); % Intensity of all individual borders in all embryos
Intensity_average = zeros(1,9);

n_0_90 = 0; % Number of all individual borders in all embryos
n10=0; % Number of individual borders 0-10° in all embryos
n40=0; % Number of individual borders 10-40° in all embryos
n90=0; % Number of individual borders 40-90° in all embryos

precision = 2;
x = -90:0.5:90;
%dilation of vertices
se90V = strel('line', 8, 90);
se0V = strel('line', 7, 0);
%dilation of borders
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);
% dilation cell
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);

for g=1:numel(files_tif)
    %% Resetting counters
    cells_pos = 0; %Number of positively oriented cells
    cells_neg = 0; %Number of negatively oriented cells
    cell_orientation_pos = 0; %Avarage orientation of positively oriented cells
    cell_orientation_neg = 0; %Average orientation of negatively oriented cells
    cell_orientation = 0; %Final average orientation of all cells
    BG = 0; %Average intensity of background
    q=0; % First complete cell adjacent to the border
    q2=0; %Second complete cell adjacent to the border
    celldata = zeros(2,7); %Array with data of individual cells
    intensity = zeros(2,2); %Array with data of individual borders
    
    
    %% Open images and modify
    Name = [num2str(g),'.tif'];
    cd(tifCad);
    Cad = imread(Name);
    Cad2 = uint16(Cad);
    cd(tifAB);
    AB = imread(Name);
    AB2 = uint16(AB);
%     cd(tifRab);
%     Rab = imread(Name);
%     Rab2 = uint16(Rab);
    bd_dir = [tif8,'/', num2str(g)];
    cd(bd_dir);
    %Open image with vertices
    V = imread('vertices.tif');
    V = imbinarize(rgb2gray(V),0);
    % Vdil - dilated image of all vertices
    Vdil = imdilate(V, [se90V se0V]);
    % Individual vertices as objects
    cc_Vdil = bwconncomp(Vdil);
    L_Vdil = labelmatrix(cc_Vdil);
    
    I=imread('handCorrection.tif');
    I2=imbinarize(rgb2gray(I),0);
    I2(:,1) = 0;
    I2(:,end) = 0;
    I2(1,:) = 0;
    I2(end,:) = 0;
    I3 = imdilate(I2, [se90I se0I]);
    
    % I_cells - inverted image of all cells that are completely in frame;
    % s_cells - individual cells as objects
    
    I_cells=imcomplement(I3);
    I_cells=imclearborder(I_cells);
    cc_cells=bwconncomp(I_cells);
    L_cells=labelmatrix(cc_cells);
    s_cells=regionprops(cc_cells, Cad, 'MeanIntensity', 'Orientation', 'PixelIdxList', 'Eccentricity', 'Area');
    
    %Cyto_Rab = regionprops(cc_cells, Rab, 'PixelValues');
    Cyto_AB = regionprops(cc_cells, AB, 'PixelValues');
    
    for i=1:numel(s_cells)
        
    end
    
    % I_borders - all cell-cell borders of the cells that are completely in
    % frame; Junction_cad - individual borders as objects
    I_borders = imsubtract(I3, Vdil);
    I_borders = imbinarize(I_borders, 0);
    I_borders=imclearborder(I_borders);
    cc_borders=bwconncomp(I_borders);
    Junction_Cad=regionprops(cc_borders,Cad,'MeanIntensity','Orientation','Perimeter', 'centroid', 'Extrema', 'PixelList');
    Junction_AB=regionprops(cc_borders,AB,'MeanIntensity','Orientation','Perimeter', 'centroid', 'Extrema', 'PixelList');
    %% Image data: BG and orientation
    % Background and image orientation: measuring summarized background and
    % summarized orietnations of positively and negatively oriented cells
    for i=1:numel(s_cells)
        if s_cells(i).Area>smallA && s_cells(i).Area<largeA
            if s_cells(i).Orientation>0
                cell_orientation_pos=cell_orientation_pos+s_cells(i).Orientation;
                cells_pos=cells_pos+1;
            else
                cell_orientation_neg=cell_orientation_neg+s_cells(i).Orientation;
                cells_neg=cells_neg+1;
            end
            BG=BG+s_cells(i).MeanIntensity;
        end
    end
    
    % Average orientation of positively and negatively oriented cells
    if cells_pos>0
        cell_orientation_pos=cell_orientation_pos / cells_pos;
    end
    
    if cells_neg>0
        cell_orientation_neg=cell_orientation_neg / cells_neg;
    end
    
    %Getting summarized orientation of the all cells
    if abs(cell_orientation_pos-cell_orientation_neg)>90
        for i=1:numel(s_cells)
            if s_cells(i).Area>smallA && s_cells(i).Area<largeA
                if s_cells(i).Orientation<0
                    cell_orientation=cell_orientation+180+s_cells(i).Orientation;
                else
                    cell_orientation=cell_orientation+s_cells(i).Orientation;
                end
            end
            
        end
    else
        for i=1:numel(s_cells)
            cell_orientation=cell_orientation+s_cells(i).Orientation;
        end
    end
    
    %Final average orientation of all cells and background
    cell_orientation=cell_orientation/(cells_pos + cells_neg);
    BG=BG/(cells_pos + cells_neg);
    
    %% Borders data
    k=0;
    Intensity10 = 0;
    Intensity40 = 0;
    Intensity90 = 0;
    Intensity10_BG = 0;
    Intensity40_BG = 0;
    Intensity90_BG = 0;
    k10=0;
    k40=0;
    k90=0;
    cutoff = 0;
    
    for n=1:numel(Junction_Cad)
        % Finding borders between two cells that are completely in frame
        boundary = bwtraceboundary(I_borders, [Junction_Cad(n).Extrema(3,2) + 0.5  Junction_Cad(n).Extrema(3,1) - 0.5], 'NW');
        boundary_value = zeros(1,length(boundary*4));
        for i=1:length(boundary*4)
            boundary_value(i*4-3) = L_cells(boundary(i,1)-1,boundary(i,2)-1);
            boundary_value(i*4-2) = L_cells(boundary(i,1)-1,boundary(i,2)+1);
            boundary_value(i*4-1) = L_cells(boundary(i,1)+1,boundary(i,2)-1);
            boundary_value(i*4) = L_cells(boundary(i,1)+1,boundary(i,2)+1);
        end
        boundary_value = unique(boundary_value);
        boundary_value(boundary_value==0) = [];
        
        if length(boundary_value)==2
            q = boundary_value(1,1);
            q2 = boundary_value(1,2);
            if s_cells(q).Area>smallA && s_cells(q2).Area>smallA && s_cells(q).Area<largeA && s_cells(q2).Area<largeA...
                    && Junction_Cad(n).Perimeter/2>5
                
                % Normalizing orientation relative to image
                % (average cell) orientation)
                if abs(Junction_Cad(n).Orientation - cell_orientation)>90
                    if cell_orientation > 0
                        Junction_Cad(n).Orientation = abs(cell_orientation - 180 - Junction_Cad(n).Orientation);
                        Junction_AB(n).Orientation = abs(cell_orientation - 180 - Junction_AB(n).Orientation);
                    else
                        Junction_Cad(n).Orientation = abs(cell_orientation + 180 - Junction_Cad(n).Orientation);
                        Junction_AB(n).Orientation = abs(cell_orientation + 180 - Junction_AB(n).Orientation);
                    end
                else
                    Junction_Cad(n).Orientation = abs(cell_orientation - Junction_Cad(n).Orientation);
                    Junction_AB(n).Orientation = abs(cell_orientation - Junction_AB(n).Orientation);
                end
                % Writing down parametes of individual junctions
                % that are longer than 5px and are surrounded by
                % cells with area>1000px and are completely within
                % the image
                k=k+1;
                intensity(k,1) = Junction_Cad(n).Orientation;
                intensity(k,2) = Junction_Cad(n).MeanIntensity;
                intensity(k,3) = Junction_AB(n).MeanIntensity;
                intensity(k,4) = Junction_AB(n).MeanIntensity/Junction_Cad(n).MeanIntensity;
                               
                n_0_90 = n_0_90 + 1;
                Intensity_total(n_0_90,1:4) = intensity(k,1:4);
            end
            
        end
    end
    
     % Sorting borders by angle
    intensity=sortrows(intensity,1);
    % information about individual borders
    cd(sum_dir);
    Otput_Celldata = [num2str(g),'.csv'];
    headers = {'angle', 'GFP', 'AB', 'AB/GFP'};
    csvwrite_with_headers(Otput_Celldata,intensity,headers);
    
    Intensity_average(g,1) = g;
    Intensity_average(g,2:5) = mean(intensity(intensity(:,1)<=10,:),1);
    Intensity_average(g,6:9) = mean(intensity(intensity(:,1)>40,:),1);
    
    
end

Intensity_total=sortrows(Intensity_total,1);
csvwrite_with_headers('All.csv',Intensity_total,headers);

headers = {'embryo','angle', 'GFP 0-10', 'AB 0-10', 'AB/GFP 0-10', 'angle' ,'GFP 40-90','AB 40-90', 'AB/GFP 40-90'};
csvwrite_with_headers('Embryo.csv',Intensity_average,headers);

cd(currdir);
close all;
clear variables;
clc

