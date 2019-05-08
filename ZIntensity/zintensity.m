%% Clear all and initial parameters
clc
clear variables
close all


%Folder to save information about cells
if exist([filedir, '/Zdata'],'dir') == 0
    mkdir(filedir,'/Zdata');
end
result_dir = [filedir, '/Zdata'];

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8_dir =[filedir, '/borders'];
oib_dir = [filedir, '/originals'];

step = 0.38;
x = step*(0:1:5);

%Reading 16-bit average intensity projection files
cd(oib_dir);
files_oib = dir('*.oib');

%dilation of borders
se90I = strel('line', 2, 90);
se0I = strel('line', 2, 0);
MeanBorder = zeros(numel(files_oib),6);
SDBorder = zeros(numel(files_oib),6);
Width = zeros(numel(files_oib),1);
Centre = zeros(numel(files_oib),1);
gof = zeros(numel(files_oib),1);
for g=1:numel(files_oib)
    
     %% Open images and modify
    cd(oib_dir);
    Number1 = [num2str(g),'.oib'];
    OIB=bfopen(Number1);

    Series = OIB{1,1};
    seriesCount = size(Series, 1)/2; %display size to check type of file
    Series_plane1= Series{1,1};
     
     
    bd_dir = [tif8_dir,'/', num2str(g)];
    cd(bd_dir);
    
    I=imread('handCorrection.tif');
    I2=imbinarize(rgb2gray(I),0);
    I2(:,1) = 0;
    I2(:,end) = 0;
    I2(1,:) = 0;
    I2(end,:) = 0;
    I3 = imdilate(I2, [se90I se0I]);
    Border = zeros(length(I3(I3>0)),seriesCount);
    for q=1:seriesCount
        
        Series_plane = Series{q*2-1,1};
        Series_border = im2double(Series_plane) .* I3;
        Series_border = Series_border(:);
        Series_border(Series_border == 0) = [];
        Border(:,q) = Series_border;
    end
    Border = Border ./ max(Border,[],2);
    MeanBorder(g,:) = mean(Border,1);
    SDBorder(g,:) = std(Border,1);
    [fittemp, goftemp] = fit(x',MeanBorder(g,:)','gauss1');
    gof(g) = goftemp.rsquare;
    Centre(g) = fittemp.b1;
    Width(g) = 2*sqrt(2*log(2))*fittemp.c1;
end

cd(result_dir);

csvwrite('Meanbysection.csv',[x',MeanBorder']);
csvwrite('SDbysection.csv',[x',SDBorder']);
headers = {'Width', 'Centre', 'gof'};
csvwrite_with_headers('Fit.csv',[Width,Centre,gof], headers);

cd(currdir);
clc;
clear variables;