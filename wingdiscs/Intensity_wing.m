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
se90I = strel('line', 20, 90);
se0I = strel('line', 20, 0);
ch= 1;
cutoff = 0.9;
bin_size = 0.15;
fileformat = '.oib';
parameters = inputdlg({'The channel to analyse:','Cutoff for the goodness of fit:',...
    'z-section distance','Format of original files'},...
        'Parameters',1,{num2str(ch), num2str(cutoff), num2str(bin_size), fileformat});
    % Redefine extension
    ch = str2double(parameters{1});
    cutoff = str2double(parameters{2});
    bin_size = str2double(parameters{3});
    fileformat = parameters{4};
    
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

if exist([filedir, '/Distributions'],'dir') == 0
    mkdir(filedir,'/Distributions');
end
Dist_dir = [filedir, '/Distributions'];

if exist([filedir, '/Distributions_fit'],'dir') == 0
    mkdir(filedir,'/Distributions_fit');
end
Fit_dir = [filedir, '/Distributions_fit'];



%Reading 16-bit average intensity projection files
cd(oib_dir);
files_tif = dir(['*', fileformat]);

averagedata = zeros(numel(files_tif),15);
averagedata2 = zeros(numel(files_tif),15);
pulled = zeros(1,9);
pulled2 = zeros(1,9);

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
    'Total','SDtotal', 'Zposition', 'SDZ',...
    'Width','SDwidth', 'pvalue','SDp', 'Area','SDarea','Eccentricity','SDecc'};
writetable(T2, 'avaragesbywing.csv');

pulled = pulled(2:end,:);
%pulled2 = pulled(

T3 = array2table(pulled);

T3.Properties.VariableNames = {'wing','cell', 'Mean', 'Total', 'Zposition', ...
    'Width', 'pvalue', 'Area','Eccentricity'};
writetable(T3, 'pulleddata.csv');

T4 = array2table(averagedata2);

T4.Properties.VariableNames = {'wing', 'Mean','SDmean',...
    'Total','SDtotal', 'Zposition', 'SDZ',...
    'Width','SDwidth', 'pvalue','SDp', 'Area','SDarea','Eccentricity','SDecc'};
writetable(T4, 'avaragesbywing_cutoff.csv');

pulled2 = pulled2(2:end,:);

T5 = array2table(pulled2);

T5.Properties.VariableNames = {'wing','cell', 'Mean', 'Total', 'Zposition', ...
    'Width', 'pvalue', 'Area','Eccentricity'};
writetable(T5, 'pulleddata_cutoff.csv');

cd(currdir);
clear variables;
close all;
clc;