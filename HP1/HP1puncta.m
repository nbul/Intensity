%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);


HP1 = '7DAPI.tif';
HP1_im = imread(HP1);
HP1_im = imadjust(imgaussfilt(HP1_im,2));
level = graythresh(HP1_im);
%level = level * 1.7;

BW = double(imbinarize(HP1_im,level));
imshow(BW);

Mask = '7mask.tif';
Mask_im = imread(Mask);
Mask_im = double(imbinarize(Mask_im));

Ecadpos = double(immultiply(Mask_im, BW));
Ecadneg = double(immultiply(imcomplement(Mask_im), BW));

C = imfuse(double(HP1_im), imclearborder(BW), 'ColorChannels', 'green-magenta');


% 19.6078 pixels is 1 µm
% 1px is 0.051 µm in length
% 1px is 0.0026 µm2 in area

%% E-cad-positive HP1 puncta
%Ecadpos=imclearborder(Ecadpos);
cc_Ecadpos=bwconncomp(Ecadpos);
EcadposPuncta=regionprops(cc_Ecadpos,HP1_im,'Area','MeanIntensity');

Pos = [EcadposPuncta.Area];
SizePos = mean(Pos) * 0.0026;
NumPos = length(Pos)/(sum(Mask_im(:)) * 0.0026);
ConcPos = sum(Ecadpos(:))/sum(Mask_im(:)) * 100;
IntPos = mean([EcadposPuncta.MeanIntensity]);

%% E-cad-negative HP1 puncta
%Ecadneg=imclearborder(Ecadneg);
cc_Ecadneg=bwconncomp(Ecadneg);
EcadnegPuncta=regionprops(cc_Ecadneg,HP1_im,'Area', 'MeanIntensity');

Neg = [EcadnegPuncta.Area];
SizeNeg = mean(Neg) * 0.0026;
NumNeg = length(Neg)/(sum(imcomplement(Mask_im(:))) * 0.0026);
ConcNeg = sum(Ecadneg(:))/sum(imcomplement(Mask_im(:))) * 100;
IntNeg = mean([EcadnegPuncta.MeanIntensity]);