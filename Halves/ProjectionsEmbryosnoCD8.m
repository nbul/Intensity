%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(currdir);
filedir = uigetdir();
cd(filedir);

color0 = questdlg(strcat('Which color is E-cad'),'Settings','Green', 'Red', 'Red');
if strcmp(color0, 'Green')
    Cadcolor = 0;
else 
    Cadcolor = 1;
end

color0 = questdlg(strcat('Which color to quntify'),'Settings','Green', 'Red', 'Green');
if strcmp(color0, 'Green')
    Bcolor = 0;
else
    Bcolor = 1;
end

%Folder with originals
if exist([filedir, '/originals'],'dir') == 0
    mkdir(filedir, '/originals');
end
or_dir = [filedir, '/originals'];

%% folders for projections
if exist([filedir, '/borders'],'dir') == 0
    mkdir(filedir, '/borders');
end
b_dir = [filedir, '/borders'];

if exist([filedir, '/tifs_original'],'dir') == 0
    mkdir(filedir, '/tifs_original');
end
i_dir = [filedir, '/tifs_original'];

files = dir('*.oib');
%% Projections
for i=1:numel(files)
    name.original = [num2str(i), '.oib'];
    original = bfopen(name.original);
    
    Series = original{1,1};
    seriesCount = size(Series, 1); %display size to check type of file
    Intens = zeros(size(Series{1,1},1),size(Series{1,1},2));
    Border = zeros(size(Series{1,1},1),size(Series{1,1},2));
    for plane = 1:(seriesCount/2)
        Intens = plus(Intens,double(Series{plane*2-Cadcolor,1}));
        Border = max(Border, double(Series{plane*2-Bcolor,1}));
    end   
    Intens = uint16(Intens*2/seriesCount);
    Border = imadjust(uint16(Border*2/seriesCount));
    BG = imopen(Border,strel('disk',50));
    Border = Border - BG;
    Border = imadjust(uint16(Border));
    cd(b_dir);
    imwrite(Border, [num2str(i), '.tif']);
    cd(i_dir);
    imwrite(Intens, [num2str(i), '.tif']);
    cd(filedir);
    movefile([num2str(i), '.oib'], or_dir);
end

cd(currdir);
close all
clear variables
clc