%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
b_dir =[filedir, '/borders']; 
tif_dir = [filedir, '/tifs_original'];

%Reading 16-bit average intensity projection files
cd(tif_dir);
files_tif = dir('*.tif');

summary = zeros(numel(files_tif),16);

for k=1:numel(files_tif)
    cd(tif_dir);
    Signal=imread([num2str(k),'.tif']);
    bd_dir = [b_dir,'/', num2str(k)];
    cd(bd_dir);
    I=imread('handCorrection.tif');
    I=imbinarize(rgb2gray(I),0);
    I(:,1) = 0;
    I(:,end) = 0;
    I(1,:) = 0;
    I(end,:) = 0;
    [im_x, im_y] = size(I);
    
    [B,L,N,A] = bwboundaries(I,'holes');

    %preallocation
    counter = 0;
    MICyto = 0;
    Area = 0;
    TICyto = 0;
    Perimeter = 0;
    MIBorder = 0;
    TIBorder = 0;
    
    for l=1:numel(B)
        if sum(A(:,l)) > 0
            Itemp = imdilate(poly2mask(B{l}(:,2),B{l}(:,1),im_x,im_y), strel('diamond', 1));
            I2 = I .* Itemp;
            I2 = bwareaopen(I2,100);  
            [B2,L2, N2, A2] = bwboundaries(I2,'holes');
            im_cells_data=regionprops(L2,Signal,'Eccentricity', 'Area','Centroid', 'MeanIntensity');
            for n=1:numel(im_cells_data)
                if sum(A2(:,n)) == 0
                    counter = counter + 1;
                    MICyto(counter) = im_cells_data(n).MeanIntensity; % Mean intensity of cytoplasm
                    Area(counter) = im_cells_data(n).Area; % Area
                    TICyto(counter) = MICyto(counter) .* Area(counter); % Total intensity cytoplasm
                    Perimeter(counter) = length(B2{n});
                    
                    Itemp2 = zeros(im_x,im_y);
                    for m = 1:Perimeter(counter)
                        Itemp2(B2{n}(m,1),B2{n}(m,2)) = 1;
                    end
                    Itemp2 = imdilate(Itemp2, strel('diamond', 3));
                    Itemp2 = double(Itemp2) .* double(Signal);
                    Itemp2 = Itemp2(:);
                    Itemp2(Itemp2 == 0) = [];
                    MIBorder(counter) = mean(Itemp2);
                    TIBorder(counter) = sum(Itemp2);
                end
            end

        end
    end
   
    summary(k,:) = [k, mean(MICyto), std(MICyto), mean(TICyto), std(TICyto),mean(Area), std(Area),...
        mean(MIBorder), std(MIBorder), mean(TIBorder), std(TIBorder),mean(Perimeter), std(Perimeter),...
        mean(TICyto./TIBorder), std(TICyto./TIBorder), counter];
    
end

cd(filedir);
all = array2table(summary);
all.Properties.VariableNames = {'Wing', 'MeanCyto', 'stdMeanCyto', 'TotalCyto', 'sdtTotalCyto', 'Area', 'stdArea',...
    'MeanBorders', 'sdtMeanBorders', 'TotalBorders', 'stdTotalBorders', 'Perimeter', 'stdPerimeter',...
    'RatioCytoAJ','stdRatio','Cells'};
writetable(all,'Intensity.csv');
cd(currdir);

close all;
clear variables;
clc;