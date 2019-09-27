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

ch= 1;
parameters = inputdlg({'The channel to analyse'},...
        'Parameters',1,{num2str(ch)});
    % Redefine extension
    cutoff = str2double(parameters{1});
    
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
files_tif = dir('*.oib');

for k=1:numel(files_tif)
    readoriginal;
    
end

