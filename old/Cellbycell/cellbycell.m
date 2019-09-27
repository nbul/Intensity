clc
clear variables
close all
%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);


GFP_dir =[filedir, '/GFP'];
Cherry_dir =[filedir, '/Cherry'];
border_dir =[filedir, '/borders'];
cd(GFP_dir);
files = dir('*.tif');
cd(currdir);


if exist([filedir,'/summary'],'dir') == 0
    mkdir(filedir,'/summary');
end
sum_dir = [filedir, '/summary'];

if exist([filedir,'/summary/embryos'],'dir') == 0
    mkdir(filedir,'/summary/embryos');
end
sum_embryo_dir = [filedir, '/summary/embryos'];

if exist([filedir,'/images'],'dir') == 0
    mkdir(filedir,'/images');
end
im_dir = [filedir, '/images'];


headers1 = {'GFP','Cherry'};
headers2 = {'GFP','Cherry', 'area'};
alldata = [0,0];
alldata2 = [0,0,0];
for loop=1:length(files)
    cd(GFP_dir);
    GFP_file = [num2str(loop),'.tif'];
    ImageGFP = imread(GFP_file);
    cd(Cherry_dir);
    Cherry_file = [num2str(loop),'.tif'];
    ImageCherry = imread(Cherry_file);
    ImageCherryGaus = imgaussfilt(ImageCherry,2);
    cd(border_dir);
    Path = [border_dir, '/', num2str(loop),'/'];
    cd(Path);
    Image_borders = imread('handCorrection.tif');
    [im_x, im_y] = size(ImageCherry);
    
    borders;
    
    images;
    Cherry = zeros(numel(b_valid),1);
    GFP = zeros(numel(b_valid),1);
    for k = 1:numel(b_valid)
        object_Cherry = im2double(ImageCherry) .*...
            poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        object_Cherry = object_Cherry(:);
        object_GFP = im2double(ImageGFP) .*...
            poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        object_GFP = object_GFP(:);
        Cherry(k) = mean(object_Cherry(object_Cherry>0))*255*255;
        GFP(k) = mean(object_GFP(object_GFP>0))*255*255;
    end
    image1=figure;
    plot(GFP',Cherry','o');
    cd(im_dir);
    image_filename = [num2str(loop),'_plot.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    cd(sum_embryo_dir);
    filename = [num2str(loop),'.csv'];
    csvwrite_with_headers(filename,[GFP,Cherry],headers1);
    alldata = [alldata;[GFP,Cherry]];
    
    BW = bwareaopen(imbinarize(ImageCherryGaus,'adaptive','Sensitivity', 0.4),30);
    Cherry2 = zeros(numel(b_valid),1);
    GFP2 = zeros(numel(b_valid),1);
    area = zeros(numel(b_valid),1);
    for k = 1:numel(b_valid)
        object_Cherry = im2double(ImageCherry) .* im2double(BW) .*...
            poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        object_Cherry2 = im2double(ImageCherry) .*...
            poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        object_Cherry = object_Cherry(:);
        object_GFP = im2double(ImageGFP) .* im2double(BW) .*...
            poly2mask(b_valid{k}(:,2),b_valid{k}(:,1),im_x,im_y);
        object_GFP = object_GFP(:);
        Cherry2(k) = mean(object_Cherry(object_Cherry>0))*255*255;
        GFP2(k) = mean(object_GFP(object_GFP>0))*255*255;  
        area(k) = length(object_Cherry(object_Cherry>0))/length(object_Cherry2(object_Cherry2>0))*100;
    end
    
    image1=figure;
    plot(GFP2',Cherry2','o');
    cd(im_dir);
    image_filename = [num2str(loop),'_plot2.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    image1=figure;
    plot(area',GFP','o');
    cd(im_dir);
    image_filename = [num2str(loop),'_plot3.tif'];
    print(image1, '-dtiff', '-r150', image_filename);
    close all
    cd(sum_embryo_dir);
    filename = [num2str(loop),'_2.csv'];
    csvwrite_with_headers(filename,[GFP2,Cherry2,area],headers2);
    alldata2 = [alldata2;[GFP2,Cherry2,area]];
    
    
end

alldata(alldata(:,1)==0,:) = [];
alldata2(alldata2(:,1)==0,:) = [];
cd(sum_dir);
csvwrite_with_headers('pulled.csv',alldata,headers1);
csvwrite_with_headers('pulled_thresholded.csv',alldata2,headers2);
image1=figure;
plot(alldata(:,1),alldata(:,2),'o');
print(image1, '-dtiff', '-r150', 'pulled_plot.tif');
close all

image1=figure;
plot(alldata2(:,1),alldata2(:,2),'o');
print(image1, '-dtiff', '-r150', 'pulled_plot2.tif');
close all

image1=figure;
plot(alldata2(:,3),alldata2(:,2),'o');
print(image1, '-dtiff', '-r150', 'pulled_plot3.tif');
close all

cd(currdir);
clear variables
clc