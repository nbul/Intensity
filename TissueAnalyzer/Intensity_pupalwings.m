%% Assigning memory
Cells_total = zeros(1,6); % Information about all cells from all images
Cells_number = zeros(1,2); % Number of cells in all images
n_all = 0; % Counter for all cells

Intensity_total = zeros(1,4); % Intensity of all individual borders in all embryos
Intensity_average = zeros(1,8);
Intensity_30 = zeros(1,4); % Intensity of individual borders 0-10° in all embryos
Intensity_60 = zeros(1,4); % Intensity of individual borders 10-40° in all embryos
Intensity_90 = zeros(1,4); % Intensity of individual borders 40-90° in all embryos
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
se90I = strel('line', 4, 90);
se0I = strel('line', 4, 0);
% dilation cell
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);

%% Analysis of individual images
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
    Cad = [num2str(g),'.tif'];
    cd(tif16_dir);
    Cad_im = imread(Cad);
    Cad_im2 = uint16(Cad_im);
    bd_dir = [tif8_dir,'/', num2str(g)];
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
    s_cells=regionprops(cc_cells, Cad_im, 'MeanIntensity', 'Orientation', 'PixelIdxList', 'Eccentricity', 'Area');
    
    % I_borders - all cell-cell borders of the cells that are completely in
    % frame; Junction_cad - individual borders as objects
    I_borders = imsubtract(I3, Vdil);
    I_borders = imbinarize(I_borders, 0);
    I_borders=imclearborder(I_borders);
    cc_borders=bwconncomp(I_borders);
    Junction_cad=regionprops(cc_borders,Cad_im,'MeanIntensity','Orientation','Perimeter', 'centroid', 'Extrema', 'PixelList');
    
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
    
    %% Cell data
    k=0;
    for n=1:numel(s_cells)
        if s_cells(n).Area>smallA && s_cells(n).Area<largeA
            newImg = imbinarize(rgb2gray(I),255);
            newImg( vertcat( s_cells(n).PixelIdxList ) ) = 1; % getting cell n only
            se90 = strel('line', 4, 90);
            se0 = strel('line', 4, 0);
            newImg = imdilate(newImg, [se90 se0]);
            cc_newImg=bwconncomp(newImg);
            s_newImg=regionprops(cc_newImg, newImg, 'Area', 'PixelList');
            pixel_value = zeros(numel(s_newImg(1).Area),1);
            for i=1:s_newImg(1).Area
                pixel_value(i) = L_Vdil(s_newImg(1).PixelList(i,2),s_newImg(1).PixelList(i,1));
            end
            pixel_value = unique(pixel_value);
            pixel_value(pixel_value==0) = [];
            % Writing data for cell n if it is bigger than 1000px in
            % celldata matrix
            k=k+1;
            celldata(k,1)=n;
            celldata(k,2)=s_cells(n).Area;
            celldata(k,3)=s_cells(n).MeanIntensity;
            celldata(k,4)=s_cells(n).Orientation;
            celldata(k,5) = s_cells(n).Eccentricity;
            celldata(k,6) = 1/sqrt(1-celldata(k,5)^2);
            celldata(k,7)= length(pixel_value);
        end
    end
    
    %% Summary of cells in the images
    cell_temp = celldata;
    Vertices3 = celldata(:,7);
    k=0;
    for n=1:length(Vertices3)
        if Vertices3(n)>8
            cell_temp(n-k,:) = [];
            k=k+1;
        end
        
        if Vertices3(n)<3
            cell_temp(n-k,:) = [];
            k=k+1;
        end
    end
    Cells_number(g,1) = g;
    Cells_number(g,2) = size(cell_temp,1);
    for n=1:size(cell_temp,1)
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

    
    %setting cutoff values for border intensity
    if choice == 1
        cutoff = BG;
    else
        cutoff = 0;
    end
    
    for n=1:numel(Junction_cad)
        % Finding borders between two cells that are completely in frame
        boundary = bwtraceboundary(I_borders, [Junction_cad(n).Extrema(3,2) + 0.5  Junction_cad(n).Extrema(3,1) - 0.5], 'NW');
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
            if s_cells(q).Area>smallA && s_cells(q2).Area>smallA && s_cells(q).Area<largeA && s_cells(q2).Area<largeA
                if Junction_cad(n).Perimeter/2>5
                    if Junction_cad(n).MeanIntensity>cutoff
                        % Normalizing orientation relative to image
                        % (average cell) orientation)
%                         if abs(Junction_cad(n).Orientation - cell_orientation)>90
%                             if cell_orientation > 0
%                                 Junction_cad(n).Orientation = abs(cell_orientation - 180 - Junction_cad(n).Orientation);
%                             else
%                                 Junction_cad(n).Orientation = abs(cell_orientation + 180 - Junction_cad(n).Orientation);
%                             end
%                         else
%                             Junction_cad(n).Orientation = abs(cell_orientation - Junction_cad(n).Orientation);
%                         end
                        % Writing down parametes of individual junctions
                        % that are longer than 5px and are surrounded by
                        % cells with area>1000px and are completely within
                        % the image
                        k=k+1;
                        intensity(k,1) = Junction_cad(n).Orientation;
                        intensity(k,2) = Junction_cad(n).MeanIntensity-BG;
                        intensity(k,3) = Junction_cad(n).MeanIntensity;
                        intensity(k,4) = Junction_cad(n).Perimeter/2;
                        
                    end
                    
                end
            end
        end
    end
    % Sorting borders by angle
   
    
    
    
    %% Summary of cell borders
    Intensity30 = mean(intensity(abs(intensity(:,1))<30,3));
    Intensity60 = mean(intensity(abs(intensity(:,1))>= 30 & abs(intensity(:,1))< 60,3));
    Intensity90 = mean(intensity(abs(intensity(:,1)) >= 60,3));
    Intensity30_BG = mean(intensity(abs(intensity(:,1))<30,2));
    Intensity60_BG = mean(intensity(abs(intensity(:,1))>= 30 & abs(intensity(:,1))< 60,2));
    Intensity90_BG = mean(intensity(abs(intensity(:,1)) >= 60,2));   
    
    Intensity_average(g,:) = [g, Intensity30, Intensity60, Intensity90, Intensity30_BG, Intensity60_BG, Intensity90_BG, BG];
    
    Intensity_30 = [Intensity_30; intensity(abs(intensity(:,1))<30,:)]; %#ok<AGROW>
    Intensity_60 = [Intensity_60; intensity(abs(intensity(:,1))>= 30 & abs(intensity(:,1))< 60,:)]; %#ok<AGROW>
    Intensity_90 = [Intensity_90; intensity(abs(intensity(:,1))>= 60,:)]; %#ok<AGROW>
    
    intensity(intensity(:,1)<0,1) = intensity(intensity(:,1)<0,1) + 180;
    Intensity_total = [Intensity_total; intensity]; %#ok<AGROW>
    
    %% Writing data for individual images
    % N_intensity - information about individual borders;
    cd(int_dir);
    Otput_Intensity = [num2str(g),'_intensity.csv'];
    headers = {'Angle', 'Intensity-BG', 'Intensity', 'Length'};
    csvwrite_with_headers(Otput_Intensity,intensity, headers);
    close all;
    
    % N_area - information about individual cells
    cd(cells_dir);
    Otput_Celldata = [num2str(g),'_area.csv'];
    headers = {'Cell', 'Area', 'Intensity', 'Orientation', 'Eccentricity', 'Elongation', 'Neighbours'};
    csvwrite_with_headers(Otput_Celldata,celldata,headers);

    
end

%% All intensities by angles
%Sorting arrays by angle and extending to 10000
Intensity_30(1,:) = [];
Intensity_30=sortrows(Intensity_30,1);
Intensity_30=wextend('ar','zpd', Intensity_30, (10000-length(Intensity_30(:,1))),'d');

Intensity_60(1,:) = [];
Intensity_60=sortrows(Intensity_60,1);
Intensity_60=wextend('ar','zpd', Intensity_60, (10000-length(Intensity_60(:,1))),'d');

Intensity_90(1,:) = [];
Intensity_90=sortrows(Intensity_90,1);
Intensity_90=wextend('ar','zpd', Intensity_90, (10000-length(Intensity_90(:,1))),'d');

Intensity_total(1,:) = [];
Intensity_total = sortrows(Intensity_total,1);

%Combining arrays into a single array
Intensity_angle2 = [Intensity_30 Intensity_60 Intensity_90];
%Removing extra 0s
Intensity_angle2(sum(Intensity_angle2, 2) == 0,:) = [];


%% Writing combined data
%Properties of individual cells
cd(sum_dir);
Otput_All_Celldata = 'vertices_all.csv';
headers = {'Wing', 'Cell', 'Area', 'Intensity', 'Orientation', 'Eccentricity', 'Elongation', 'Neighbours'};
csvwrite_with_headers(Otput_All_Celldata,Cells_total, headers);

%Number of cells in each embryo
Otput_cells_number = 'cell_number.csv';
headers2 = {'Wing', 'Cells'};
csvwrite_with_headers(Otput_cells_number,Cells_number, headers2);

%Combined intensities at cell borders by angle in all embryos
headers2 = {'0-30angels', 'Intensity', 'Intensity+BG', 'Length', '30-60angels', 'Intensity', 'Intensity+BG', 'Length', '60-90angels','Intensity', 'Intensity+BG', 'Length',};
csvwrite_with_headers('Intensity_angles.csv',Intensity_angle2, headers2);

%Averaged intensity by angle in individual embryos
headers = {'Wing', '0-30-BG', '30-60-BG', '60-90-BG', '0-30+BG', '30-60+BG', '60-90+BG', 'BG'};
csvwrite_with_headers('Intensity_embryo.csv',Intensity_average, headers);

%Combined all borders together
Intensity_total=sortrows(Intensity_total,1);
headers = {'Angle', 'Intensity-BG', 'Intensity', 'Length'};
csvwrite_with_headers('Intensity_total.csv',Intensity_total, headers);
