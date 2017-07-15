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

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
%Folders with images
tif8_dir =[filedir, '/borders']; 
tif16_dir = [filedir, '/tifs_original'];

%Folder to save information about cells
mkdir(filedir,'/Cells');
cells_dir = [filedir, '/Cells'];

%Folder to save images with distributions after BG subtraction 
mkdir(filedir,'/Distribution no BG');
dist1_dir = [filedir, '/Distribution no BG'];

%Folder to save images with distributions without BG subtraction
mkdir(filedir,'/Distribution with BG');
dist2_dir = [filedir, '/Distribution with BG'];

%Folder to save intensity information from individual images
mkdir(filedir,'/Intensity');
int_dir = [filedir, '/Intensity'];

%Folder to save summarised information
mkdir(filedir,'/SummaryIntensity');
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

