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
oib_dir = [filedir, '/originals'];
%dilation of borders
se90I = strel('line', 5, 90);
se0I = strel('line', 5, 0);
ch= 1;
bin_size = 0.38;
fileformat = '.oib';
parameters = inputdlg({'The channel to analyse:',...
    'z-section distance','Format of original files'},...
        'Parameters',1,{num2str(ch), num2str(bin_size), fileformat});
    % Redefine extension
    ch = str2double(parameters{1});
    bin_size = str2double(parameters{2});
    fileformat = parameters{3};
    
channel = ['C=',num2str(ch),'/'];

%Folder to save information about cells
if exist([filedir, '/Cells'],'dir') == 0
    mkdir(filedir,'/Cells');
end
cells_dir = [filedir, '/Cells'];

%Folder to save information about cells
if exist([filedir, '/Summary'],'dir') == 0
    mkdir(filedir,'/Summary');
end
sum_dir = [filedir, '/Summary'];



%Reading 16-bit average intensity projection files
cd(oib_dir);
files_tif = dir(['*', fileformat]);

averagedata = zeros(numel(files_tif),15);
pulled = zeros(1,8);

for k=1:numel(files_tif)
    readoriginal;
    borders;
    imageoutput;
    planebyplane;
    widthandmean;
end

cd(sum_dir);

T2 = array2table(averagedata);

T2.Properties.VariableNames = {'wing', 'Mean','SDmean',...
    'Total','SDtotal', 'MeanBG','SDmeanBG', 'TotalBG','SDtotalBG',...
    'Area','SDarea','Eccentricity','SDecc', 'BG', 'cells'};
writetable(T2, 'avaragesbywing.csv');

pulled = pulled(2:end,:);
%pulled2 = pulled(

T3 = array2table(pulled);

T3.Properties.VariableNames = {'wing','cell', 'Mean', 'Total',...
    'MeanBG', 'TotalBG', 'Area','Eccentricity'};
writetable(T3, 'pulleddata.csv');


cd(currdir);
clear variables;
close all;
clc;