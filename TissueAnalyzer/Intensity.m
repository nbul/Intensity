%% Clear all and initial parameters
clc
clear variables
close all

%% Default settings and script choice
usedefault2 = questdlg(strcat('Which kind of data is analysed'),'Settings','Embryo','Wing','Wing');

usedefault = questdlg(strcat('Do you want to consider borders with intensity below background'),'Settings','Yes','No','No');
if strcmp(usedefault, 'No')
    choice = 1;
else
    choice = 0;
end

usedefault3 = questdlg(strcat('Which kind of data is analysed'),'Settings','Airyscan','Olympus','Olympus');

if strcmp(usedefault3, 'Olympus')
    smallA = 300;
    largeA = 25000;
else
    smallA = 1500;
    largeA = 35000;
end

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8_dir =[filedir, '/borders']; 
tif16_dir = [filedir, '/tifs_original'];

%Folder to save information about cells
if exist([filedir, '/Cells'],'dir') == 0
    mkdir(filedir,'/Cells');
end
cells_dir = [filedir, '/Cells'];

%Folder to save images with distributions after BG subtraction 
if exist([filedir, '/Distribution no BG'],'dir') == 0
    mkdir(filedir,'/Distribution no BG');
end
dist1_dir = [filedir, '/Distribution no BG'];

%Folder to save images with distributions without BG subtraction
if exist([filedir, '/Distribution with BG'],'dir') == 0
    mkdir(filedir,'/Distribution with BG');
end
dist2_dir = [filedir, '/Distribution with BG'];

%Folder to save intensity information from individual images
if exist([filedir, '/Intensity'],'dir') == 0
    mkdir(filedir,'/Intensity');
end
int_dir = [filedir, '/Intensity'];

%Folder to save summarised information
if exist([filedir, '/SummaryIntensity'],'dir') == 0
    mkdir(filedir,'/SummaryIntensity');
end
sum_dir = [filedir, '/SummaryIntensity'];

%Reading 16-bit average intensity projection files
cd(tif16_dir);
files_tif = dir('*.tif');

if strcmp(usedefault2, 'Wing')
    Intensity_wing;
else
    Intensity_embryo;
end
cd(currdir);
close all;
clear variables;
clc

