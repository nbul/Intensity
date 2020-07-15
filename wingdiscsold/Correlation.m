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

Cad = bfopen('4_Cad.tif');
Cad = Cad{1,1};

AP1 = bfopen('4_AP1.tif');
AP1 = AP1{1,1};

SeriesCount = size(Cad, 1);

Cad_signal = [Cad{1, 1}(:)];
AP1_signal = [AP1{1, 1}(:)];

for i= 2:SeriesCount
    Cad_signal = [Cad_signal; Cad{i, 1}(:)];
    AP1_signal = [AP1_signal; AP1{i, 1}(:)];
end

[R,P] = corrcoef(double(Cad_signal), double(AP1_signal))