%% Assigning memory
Cells_total = zeros(1,6); % Information about all cells from all images
Cells_number = zeros(1,2); % Number of cells in all images
n_all = 0; % Counter for all cells

Intensity_total = zeros(1,4); % Intensity of all individual borders in all embryos
Intensity_average = zeros(1,8);
Intensity_10 = zeros(1,4); % Intensity of individual borders 0-10° in all embryos
Intensity_40 = zeros(1,4); % Intensity of individual borders 10-40° in all embryos
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
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);
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
    V = imread('vertices.png');
    V = im2bw(V,1/255);
    % Vdil - dilated image of all vertices
    Vdil = imdilate(V, [se90V se0V]);
    % Individual vertices as objects
    cc_Vdil = bwconncomp(Vdil);
    L_Vdil = labelmatrix(cc_Vdil);
    
    I=imread('tracked_bd.png');
    I2=im2bw(I,1/255);
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
    I_borders = im2bw(I_borders, 1/255);
    cc_borders=bwconncomp(I_borders);
    Junction_cad=regionprops(cc_borders,Cad_im,'MeanIntensity','Orientation','Perimeter', 'centroid', 'Extrema', 'PixelList');
    
    %% Image data: BG and orientation
    % Background and image orientation: measuring summarized background and
    % summarized orietnations of positively and negatively oriented cells
    for i=1:numel(s_cells);
        if s_cells(i).Area>300 && s_cells(i).Area<15000
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
        for i=1:numel(s_cells);
            if s_cells(i).Area>300 && s_cells(i).Area<15000
                if s_cells(i).Orientation<0
                    cell_orientation=cell_orientation+180+s_cells(i).Orientation;
                else cell_orientation=cell_orientation+s_cells(i).Orientation;
                end
            end
            
        end
    else
        for i=1:numel(s_cells);
            cell_orientation=cell_orientation+s_cells(i).Orientation;
        end
    end
    
    %Final average orientation of all cells and background
    cell_orientation=cell_orientation/(cells_pos + cells_neg);
    BG=BG/(cells_pos + cells_neg);
    
    %% Cell data
    k=0;
    for n=1:numel(s_cells)
        if s_cells(n).Area>300 && s_cells(n).Area<15000
            newImg = im2bw(I,255/255);
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
    Intensity10 = 0;
    Intensity40 = 0;
    Intensity90 = 0;
    Intensity10_BG = 0;
    Intensity40_BG = 0;
    Intensity90_BG = 0;
    k10=0;
    k40=0;
    k90=0;
    
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
            if s_cells(q).Area>300 && s_cells(q2).Area>300 && s_cells(q).Area<15000 && s_cells(q2).Area<15000
                if Junction_cad(n).Perimeter/2>5
                    if Junction_cad(n).MeanIntensity>cutoff
                        % Normalizing orientation relative to image
                        % (average cell) orientation)
                        if abs(Junction_cad(n).Orientation - cell_orientation)>90
                            if cell_orientation > 0
                                Junction_cad(n).Orientation = abs(cell_orientation - 180 - Junction_cad(n).Orientation);
                            else Junction_cad(n).Orientation = abs(cell_orientation + 180 - Junction_cad(n).Orientation);
                            end
                        else Junction_cad(n).Orientation = abs(cell_orientation - Junction_cad(n).Orientation);
                        end
                        % Writing down parametes of individual junctions
                        % that are longer than 5px and are surrounded by
                        % cells with area>1000px and are completely within
                        % the image
                        k=k+1;
                        intensity(k,1) = Junction_cad(n).Orientation;
                        intensity(k,2) = Junction_cad(n).MeanIntensity-BG;
                        intensity(k,3) = Junction_cad(n).MeanIntensity;
                        intensity(k,4) = Junction_cad(n).Perimeter/2;
                        
                        %Writing intensity into a summarizes files
                        if intensity(k,1)<=10
                            k10=k10+1;
                            Intensity10 = Intensity10 + intensity(k,2);
                            Intensity10_BG = Intensity10_BG + intensity(k,3);
                            n10=n10+1;
                            Intensity_10(n10,1) = intensity(k,1);
                            Intensity_10(n10,2) = intensity(k,2);
                            Intensity_10(n10,3) = intensity(k,3);
                            Intensity_10(n10,4) = intensity(k,4);
                        end
                        if (intensity(k,1)>10 && intensity(k,1)<=40)
                            k40=k40+1;
                            Intensity40 = Intensity40 + intensity(k,2);
                            Intensity40_BG = Intensity40_BG + intensity(k,3);
                            n40=n40+1;
                            Intensity_40(n40,1) = intensity(k,1);
                            Intensity_40(n40,2) = intensity(k,2);
                            Intensity_40(n40,3) = intensity(k,3);
                            Intensity_40(n40,4) = intensity(k,4);
                        end
                        if intensity(k,1)>40
                            k90=k90+1;
                            Intensity90 = Intensity90 + intensity(k,2);
                            Intensity90_BG = Intensity90_BG + intensity(k,3);
                            n90=n90+1;
                            Intensity_90(n90,1) = intensity(k,1);
                            Intensity_90(n90,2) = intensity(k,2);
                            Intensity_90(n90,3) = intensity(k,3);
                            Intensity_90(n90,4) = intensity(k,4);
                        end
                        n_0_90 = n_0_90 + 1;
                        Intensity_total(n_0_90,1) = intensity(k,1);
                        Intensity_total(n_0_90,2) = intensity(k,2);
                        Intensity_total(n_0_90,3) = intensity(k,3);
                        Intensity_total(n_0_90,4) = intensity(k,4);
                    end
                    
                end
            end
        end
    end
    % Sorting borders by angle
    intensity=sortrows(intensity,1);
    
    %% Summary of cell borders
    Intensity10 = Intensity10/k10;
    Intensity40 = Intensity40/k40;
    Intensity90 = Intensity90/k90;
    Intensity10_BG = Intensity10_BG/k10;
    Intensity40_BG = Intensity40_BG/k40;
    Intensity90_BG = Intensity90_BG/k90;
    Intensity_average(g,1) = g;
    Intensity_average(g,2) = Intensity10;
    Intensity_average(g,3) = Intensity40;
    Intensity_average(g,4) = Intensity90;
    Intensity_average(g,5) = Intensity10_BG;
    Intensity_average(g,6) = Intensity40_BG;
    Intensity_average(g,7) = Intensity90_BG;
    Intensity_average(g,8) = BG;
    
    %% Fitting and plotting
    equation1=polyfit(intensity(:,1),intensity(:,2),1); % With background subtraction
    equation2=polyfit(intensity(:,1),intensity(:,3),1); % Without background subtraction
    y1 = equation1(1,1)*x + equation1(1,2);
    y2 = equation2(1,1)*x + equation2(1,2);
    
    % Plot image of fitted intensities after background subtraction
    image1=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,2), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle after background subtraction', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y1 , '-b', 'LineWidth',3);
    text_control = ['embryo', num2str(g), ': ', num2str(equation1(1,2), precision) ' + ' num2str(equation1(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([0 90 0 max(intensity(:,2))]);
    
    % Plot image of fitted intensities without background subtraction
    image2=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,3), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle without background subtraction', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y2 , '-b', 'LineWidth',3);
    text_control = ['embryo', num2str(g), ': ', num2str(equation2(1,2), precision) ' + ' num2str(equation2(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([0 90 0 max(intensity(:,3))]);
    
    %% Writing data for individual images
    % Borders data: N_distribution_with_BG - after BG subtraction;
    % Borders data: N_distribution_with_BG1 - after BG subtraction, BG;
    cd(dist2_dir);
    Otput_Graph = [num2str(g),'_distribution_with_BG.tif'];
    hold off
    print(image1, '-dtiff', '-r300', Otput_Graph);
    
    % N_distribution_no_BG - without BG subtraction;
    cd(dist1_dir);                     
    Otput_Graph = [num2str(g),'_distribution_no_BG.tif'];
    hold off
    print(image2, '-dtiff', '-r300', Otput_Graph); 
    
    
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
Intensity_10=sortrows(Intensity_10,1);
Intensity_10=wextend('ar','zpd', Intensity_10, (10000-length(Intensity_10(:,1))),'d');
Intensity_40=sortrows(Intensity_40,1);
Intensity_40=wextend('ar','zpd', Intensity_40, (10000-length(Intensity_40(:,1))),'d');
Intensity_90=sortrows(Intensity_90,1);
Intensity_90=wextend('ar','zpd', Intensity_90, (10000-length(Intensity_90(:,1))),'d');

%Combining arrays into a single array
Intensity_angle = [Intensity_10 Intensity_40 Intensity_90];
%Removing extra 0s
Intensity_angle2 = Intensity_angle;
t=0;
for m=1:size(Intensity_angle, 1)
    if (Intensity_angle(m,1) == 0) && (Intensity_angle(m,2) == 0)
        Intensity_angle2(m-t, :) = [];
        t=t+1;
    end
end


%% Plotting summarized data
%Fitting summarized data with lines
equation1=polyfit(Intensity_total(:,1),Intensity_total(:,2),1); %with background subtraction
equation2=polyfit(Intensity_total(:,1),Intensity_total(:,3),1); % without background subtraction
y1 = equation1(1,1)*x + equation1(1,2);
y2 = equation2(1,1)*x + equation2(1,2);

% Plot image of fitted intensities after background subtraction
image3=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,2), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle after background subtraction', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y1 , '-b', 'LineWidth',3);
text_control = ['embryo', num2str(g), ': ', num2str(equation1(1,2), precision) ' + ' num2str(equation1(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([0 90 0 max(Intensity_total(:,2))]);

% Plot image of fitted intensities without background subtraction
image4=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,3), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle without background subtraction', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y2 , '-b', 'LineWidth',3);
text_control = ['embryo', num2str(g), ': ', num2str(equation2(1,2), precision) ' + ' num2str(equation2(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([0 90 0 max(Intensity_total(:,3))]);

%% Writing combined data
%Properties of individual cells
cd(sum_dir);
Otput_All_Celldata = 'vertices_all.csv';
headers = {'Embryo', 'Cell', 'Area', 'Intensity', 'Orientation', 'Eccentricity', 'Elongation', 'Neighbours'};
csvwrite_with_headers(Otput_All_Celldata,Cells_total, headers);

%Number of cells in each embryo
Otput_cells_number = 'cell_number.csv';
headers2 = {'Embryo', 'Cells'};
csvwrite_with_headers(Otput_cells_number,Cells_number, headers2);

%Combined intensities at cell borders by angle in all embryos
headers2 = {'0-10angels', 'Intensity', 'Intensity+BG', 'Length', '10-40angels', 'Intensity', 'Intensity+BG', 'Length', '40-90angels','Intensity', 'Intensity+BG', 'Length',};
csvwrite_with_headers('Intensity_angles.csv',Intensity_angle2, headers2);

%Averaged intensity by angle in individual embryos
headers = {'Embryo', '0-10-BG', '10-40-BG', '40-90-BG', '0-10+BG', '10-40+BG', '40-90+BG', 'BG'};
csvwrite_with_headers('Intensity_embryo.csv',Intensity_average, headers);

%Combined all borders together
Intensity_total=sortrows(Intensity_total,1);
headers = {'Angle', 'Intensity-BG', 'Intensity', 'Length'};
csvwrite_with_headers('Intensity_total.csv',Intensity_total, headers);

% Distribution_all_with_BG - after BG subtraction;
Otput_Graph = 'Distribution_all_with_BG.tif';
hold off
print(image3, '-dtiff', '-r300', Otput_Graph);

% Distribution_all_no_BG - without BG subtraction;
Otput_Graph = 'Distribution_all_no_BG.tif';
hold off
print(image4, '-dtiff', '-r300', Otput_Graph);
