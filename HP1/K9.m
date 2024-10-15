%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
files_tif = dir('*.tif');
data_pulled3 = zeros(1,4);
data_pulled4 = zeros(1,4);

for loop = 1:length(files_tif)
    Signal = imread(files_tif(loop).name);    
    x = sscanf(files_tif(loop).name,'%d.tif');
    Signal2 = imgaussfilt(Signal,2);
    T = graythresh(rescale(Signal2));
    Signal3 = imbinarize(rescale(Signal2), T*1.6);
    Signal3 = bwareaopen(Signal3,150);
    Signal3 = imclearborder(Signal3);
    Signal4 = imbinarize(rescale(Signal2), 'adaptive', 'Sensitivity', 0.2);
    Signal4 = bwareaopen(Signal4,200);
    Signal4 = imclearborder(Signal4);
    cc3 = bwconncomp(Signal3);
    cc4 = bwconncomp(Signal4);
    data3= struct2table(regionprops(cc3,Signal, 'Area','MeanIntensity'));
    data4= struct2table(regionprops(cc4,Signal, 'Area','MeanIntensity'));
    data_pulled3 = [data_pulled3; ones(size(data3,1),1)*x, (1:size(data3,1))', table2array(data3)];
    data_pulled4 = [data_pulled4; ones(size(data4,1),1)*x, (1:size(data4,1))', table2array(data4)];
end
data_pulled3(1,:) = [];
data_pulled4(1,:) = [];
data_table3 = array2table(data_pulled3);
data_table4 = array2table(data_pulled4);
data_table3.Properties.VariableNames = {'Image', 'Object','Area','MeanIntensity'};
data_table4.Properties.VariableNames = {'Image', 'Object','Area','MeanIntensity'};
writetable(data_table3, 'summaryOtsu.csv');
writetable(data_table4, 'summaryAdaptive.csv');
cd(currdir);
clc
clear variables;