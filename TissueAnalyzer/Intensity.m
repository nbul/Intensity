%% Clear all and initial parameters
clc
clear variables
close all

%% Default settings and script choice
usedefault2 = questdlg(strcat('Which kind of data is analysed'),'Settings','Embryo','Wing','Pupal wing', 'Wing');

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



%Reading 16-bit average intensity projection files
cd(tif16_dir);
files_tif = dir('*.tif');

if strcmp(usedefault2, 'Wing')
    Intensity_wing;
elseif strcmp(usedefault2, 'Embryo')
    Intensity_embryo;
else
    Intensity_pupalwings;
end
cd(currdir);
close all;
clear variables;
clc

