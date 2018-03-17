clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

se90I = strel('line', 4, 90);
se0I = strel('line', 4, 0);
im_dir = [filedir, '/tifs_original'];
border_dir = [filedir, '/borders'];
sum_dir = [filedir, '/SummaryIntensity'];

cd(im_dir);
files = dir('*.tif');
background = zeros(numel(files),1);
for i = 1:numel(files)
    cd(im_dir);
    name = [num2str(i), '.tif'];
    Image = imread(name);
    Filter1 = imcomplement(imbinarize(Image, 'adaptive'));
    cd([border_dir,'/', num2str(i)]);
    Filter2 = imread('handCorrection.tif');
    Filter2 = imbinarize(rgb2gray(Filter2),0);   
    Filter2(:,1) = 0;
    Filter2(:,end) = 0;
    Filter2(1,:) = 0;
    Filter2(end,:) = 0;   
    Filter2 = imdilate(Filter2, [se90I se0I]);
    Filter2 = imclearborder(imcomplement(Filter2));
    Image = double(Image) .* Filter1;
    Image = Image .* Filter2;
    Image = Image(:);
    background(i) = mean(Image(Image>0));
end

background = [(1:numel(files))', background];
cd(sum_dir);
csvwrite_with_headers('Background.csv',background,{'embryo','background'});
cd(currdir)
clear variables
clc

