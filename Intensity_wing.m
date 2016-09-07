%% Assigning memory
Cells_total = zeros(1,6); % Information about all cells from all images
Cells_number = zeros(1,2); % Number of cells in all images
n_all = 0; % Counter for all cells

Intensity_total = zeros(1,4); % Intensity of all individual borders in all wings
Intensity_average = zeros(1,4);
n_0_90 = 0; % Number of all individual borders in all wings

precision = 2;
x = -90:0.5:90;
%dilation of vertices
se90V = strel('line', 10, 90);
se0V = strel('line', 9, 0);
%dilation of borders
se90I = strel('line', 5, 90);
se0I = strel('line', 5, 0);
% dilation cell
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);

%% Analysis of individual images
for g=1:numel(files_tif)
    
    %% Resetting counters
    cell_orientation = 0; %Final average orientation of all cells
    BG=0;
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
    
    %Open image with borders
    I=imread('tracked_bd.png');
    I2=im2bw(I,1/255);
    I3 = imdilate(I2, [se90I se0I]);
    cd(currdir);
    
    % I_cells - inverted image of all cells that are completely in frame;
    % s_cells - individual cells as objects

    I_cells = imcomplement(I3);
    I_cells = imclearborder(I_cells);
    cc_cells = bwconncomp(I_cells);
    L_cells = labelmatrix(cc_cells);
    s_cells = regionprops(cc_cells, Cad_im, 'MeanIntensity', 'Orientation', 'PixelIdxList', 'Eccentricity', 'Area');
    
    % I_borders - all cell-cell borders of the cells that are completely in
    % frame; Junction_cad - individual borders as objects
    I_borders = imsubtract(I3, Vdil);
    I_borders = im2bw(I_borders, 1/255);
    cc_borders = bwconncomp(I_borders);
    Junction_cad = regionprops(cc_borders,Cad_im,'MeanIntensity','Orientation','Perimeter', 'centroid', 'Extrema', 'PixelList');
    
    %% Image data: BG and orientation
    % Background and image orientation: measuring summarized background and
    % summarized orietnations of positively and negatively oriented cells
    for i=1:numel(s_cells);
        BG=BG+s_cells(i).MeanIntensity;
    end
    
    BG=BG/numel(s_cells);
    
    %% Information about all cells that are completely within image
    % celldata will comtain information about only this image and recorded
    % at the end of cycle step
    k=0;
    for n=1:numel(s_cells)
        newImg = im2bw(I,255/255);
        newImg( vertcat( s_cells(n).PixelIdxList ) ) = 1; % getting cell n only
        newImg = imdilate(newImg, [se90I se0I]);
        cc_newImg=bwconncomp(newImg);
        s_newImg=regionprops(cc_newImg, newImg, 'Area', 'PixelList');
        pixel_value = zeros(numel(s_newImg(1).Area),1);
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
    % Cells_total will comtain information about all this image and
    % recorded at the end of complete cycle
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
        Size_im = size(Cad_im);
        for i=1:length(boundary*4)
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
        % Getting data about borders that are between two cells that are
        % completely in frame 
        if length(boundary_value)==2
            q = boundary_value(1,1);
            q2 = boundary_value(1,2);
            if Junction_cad(n).Perimeter/2>5 
                if Junction_cad(n).MeanIntensity>cutoff
                    % Writing down parametes of individual junctions
                    % that are longer than 5px and are surrounded by
                    % cells with area>1000px and are completely within
                    % the image
                    k=k+1;
                    intensity(k,1) = Junction_cad(n).Orientation;
                    intensity(k,2) = Junction_cad(n).MeanIntensity-BG;
                    intensity(k,3) = Junction_cad(n).MeanIntensity;
                    intensity(k,4) = Junction_cad(n).Perimeter/2;
                    %Summarised data for all wings
                    n_0_90 = n_0_90 + 1;
                    Intensity_total(n_0_90,1) = intensity(k,1);
                    Intensity_total(n_0_90,2) = intensity(k,2);
                    Intensity_total(n_0_90,3) = intensity(k,3);
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
     Intensity_average(g,4) = BG;
    
    %% Fitting and plotting
    equation1=polyfit(intensity(:,1),intensity(:,2),1); % With background subtraction BG
    equation2=polyfit(intensity(:,1),intensity(:,3),1); % Without background subtraction
    y1 = equation1(1,1)*x + equation1(1,2);
    y2 = equation2(1,1)*x + equation2(1,2);
    
    
    % Plot image of fitted intensities after background subtraction BG2
    image1=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,2), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle after background subtraction', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y1 , '-b', 'LineWidth',3);
    text_control = ['wing', num2str(g), ': ', num2str(equation1(1,2), precision) ' + ' num2str(equation1(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([-90 90 0 max(intensity(:,2))]);
    
    % Plot image of fitted intensities without background subtraction
    image2=figure;
    set(axes,'FontSize',16);
    plot(intensity(:,1),intensity(:,3), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
    title('Fluorescence vs relative angle without background subtraction', 'fontsize',18,'fontweight','b')
    xlabel('Relative angle', 'fontsize',16,'fontweight','b');
    ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
    hold on
    plot(x, y2 , '-b', 'LineWidth',3);
    text_control = ['wing', num2str(g), ': ', num2str(equation2(1,2), precision) ' + ' num2str(equation2(1,1), precision) '*x' ];
    text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
    axis([-90 90 0 max(intensity(:,3))]);
    
    %% Writing data for individual images
    
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


%% Plotting summarized data
%Fitting summarized data with lines
equation1=polyfit(Intensity_total(:,1),Intensity_total(:,2),1); % with background subtraction
equation2=polyfit(Intensity_total(:,1),Intensity_total(:,3),1); % without background subtraction
y1 = equation1(1,1)*x + equation1(1,2);
y2 = equation2(1,1)*x + equation2(1,2);

% Plot image of fitted intensities after background subtraction
image4=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,2), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle with background subtraction', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y1 , '-b', 'LineWidth',3);
text_control = ['all wings: ', num2str(equation1(1,2), precision) ' + ' num2str(equation1(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([-90 90 0 max(Intensity_total(:,2))]);

% Plot image of fitted intensities without background subtraction
image5=figure;
set(axes,'FontSize',16);
plot(Intensity_total(:,1),Intensity_total(:,3), 'o', 'Color','b', 'MarkerSize',4, 'MarkerFaceColor', 'b');
title('Fluorescence vs relative angle without background subtraction', 'fontsize',18,'fontweight','b')
xlabel('Relative angle', 'fontsize',16,'fontweight','b');
ylabel('Fluorescence intensity', 'fontsize',16,'fontweight','b');
hold on
plot(x, y2 , '-b', 'LineWidth',3);
text_control = ['all wings: ', num2str(equation2(1,2), precision) ' + ' num2str(equation2(1,1), precision) '*x' ];
text(98, 50, text_control, 'HorizontalAlignment','right', 'fontsize',14, 'fontweight','b');
axis([-90 90 0 max(Intensity_total(:,3))]);



%% Writing combined data
%Properties of individual cells
cd(sum_dir);
Otput_All_Celldata = 'vertices_all.csv';
headers = {'Wing', 'Cell', 'Area', 'Intensity', 'Orientation', 'Eccentricity', 'Elongation', 'Neighbours'};
csvwrite_with_headers(Otput_All_Celldata,Cells_total, headers);

%Number of cells in each wing
Otput_cells_number = 'cell_number.csv';
headers2 = {'Disc', 'Cells'};
csvwrite_with_headers(Otput_cells_number,Cells_number, headers2);


%Combined all borders together
Intensity_total=sortrows(Intensity_total,1);
headers = {'Angle', 'Intensity-BG', 'Intensity', 'Length'};
csvwrite_with_headers('Intensity_total.csv',Intensity_total, headers);

% Distribution_all_with_BG - after BG subtraction;
Otput_Graph = 'Distribution_all_with_BG.tif';
hold off
print(image4, '-dtiff', '-r300', Otput_Graph);

% Distribution_all_no_BG - without BG subtraction;
Otput_Graph = 'Distribution_all_no_BG.tif';
hold off
print(image5, '-dtiff', '-r300', Otput_Graph);

%Averaged intensity by angle in individual embryos
headers = {'Wing', 'Intensity-BG', 'Intensity', 'BG'};
csvwrite_with_headers('Intensity_wing.csv',Intensity_average, headers);
cd(currdir);

