%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(currdir);
javaaddpath(currdir);
filedir = uigetdir();
cd(filedir);

genotype = 'CD8';

parameters = inputdlg({'Genotype'},'Parameters',1,{num2str(genotype)});
% Redefine extension
genotype = parameters{1};

%Folder with CD8
if exist([filedir, '/', genotype],'dir') == 0
    mkdir(filedir, ['/', genotype]);
end
exp_dir = [filedir, '/', genotype];

%Folder with noCD8
if exist([filedir, '/no', genotype],'dir') == 0
    mkdir(filedir, ['/no', genotype]);
end
contr_dir = [filedir, '/no', genotype];

%Folder with originals
if exist([filedir, '/originals'],'dir') == 0
    mkdir(filedir, '/originals');
end
or_dir = [filedir, '/originals'];

%% folders for projections
cd(exp_dir);
if exist([exp_dir, '/borders'],'dir') == 0
    mkdir(exp_dir, '/borders');
end
expb_dir = [exp_dir, '/borders'];

if exist([exp_dir, '/tifs_original'],'dir') == 0
    mkdir(exp_dir, '/tifs_original');
end
expi_dir = [exp_dir, '/tifs_original'];

cd(contr_dir);
if exist([contr_dir, '/borders'],'dir') == 0
    mkdir(contr_dir, '/borders');
end
contrb_dir = [contr_dir, '/borders'];

if exist([contr_dir, '/tifs_original'],'dir') == 0
    mkdir(contr_dir, '/tifs_original');
end
contri_dir = [contr_dir, '/tifs_original'];

cd(filedir);
files = dir('*.oib');
%% Projections
for i=1:numel(files)
    name.original = [num2str(i), '.oib'];
    original = bfopen(name.original);
    
    Series = original{1,1};
    seriesCount = size(Series, 1); %display size to check type of file
    Cad.average = zeros(size(Series{1,1},1),size(Series{1,1},2));
    Cad.max = zeros(size(Series{1,1},1),size(Series{1,1},2));
    CD8.average = zeros(size(Series{1,1},1),size(Series{1,1},2));
    for plane = 1:(seriesCount/2)
        Cad.average = plus(Cad.average,double(Series{plane*2-1,1}));
        CD8.average = plus(CD8.average,double(Series{plane*2,1}));
        Cad.max = max(Cad.max, double(Series{plane*2-1,1}));
    end   
    Cad.average = uint16(Cad.average*2/seriesCount);
    CD8.average = imadjust(uint16(CD8.average*2/seriesCount));
    Cad.background = imopen(Cad.max,strel('disk',50));
    Cad.max = Cad.max - Cad.background;
    Cad.max = imadjust(uint16(Cad.max));
    CD8.average = imgaussfilt(CD8.average,2);
    thresh = graythresh(im2double(CD8.average));
    CD8.threshold = imbinarize(im2double(CD8.average),thresh); 
    CD8.threshold = imfill(CD8.threshold, 'holes');
    CD8.threshold = bwareaopen(CD8.threshold,10);
    CD8.threshold = imdilate(CD8.threshold,strel('disk',18));
    CD8.threshold = imfill(CD8.threshold, 'holes');
    Cad.CD8 = Cad.max;
    Cad.CD8(CD8.threshold == 0) = 0;
    CD8.threshold = imerode(CD8.threshold,strel('disk',24));
    Cad.noCD8 = Cad.max;
    Cad.noCD8(CD8.threshold == 1) = 0;
    cd(expb_dir);
    imwrite(Cad.CD8, [num2str(i), '.tif']);
    cd(contrb_dir);
    imwrite(Cad.noCD8, [num2str(i), '.tif']);
    cd(expi_dir);
    imwrite(Cad.average, [num2str(i), '.tif']);
    cd(contri_dir);
    imwrite(Cad.average, [num2str(i), '.tif']);
    cd(filedir);
    movefile([num2str(i), '.oib'], or_dir);
end

cd(currdir);
close all
clear variables
clc