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
tif8_dir =[filedir, '/borders_bg']; 
tif16_dir = [filedir, '/tifs_original'];
sum_dir = [filedir, '/SummaryIntensity'];
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);


%Reading 16-bit average intensity projection files
cd(tif16_dir);
files_tif = dir('*.tif');
background = zeros(1,numel(files_tif));

for g=1:numel(files_tif)
    Cad = [num2str(g),'.tif'];
    cd(tif16_dir);
    Cad_im = imread(Cad);
    Cad_im2 = uint16(Cad_im);
    bd_dir = [tif8_dir,'/', num2str(g)];
    cd(bd_dir);
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
    I_cells=bwareaopen(I_cells, 10000);
    cc_cells=bwconncomp(I_cells);
    L_cells=labelmatrix(cc_cells);
    s_cells=regionprops(cc_cells, Cad_im, 'MeanIntensity');
    background(g)=s_cells.MeanIntensity;
end

cd(sum_dir);
headers = {'embryo','background'};
Otput_cells_number = 'background_stripe.csv';
csvwrite_with_headers(Otput_cells_number,[(1:1:numel(files_tif))',background'], headers);

cd(currdir);
close all;
clear variables;
clc


