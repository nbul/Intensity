%% Clear all and initial parameters
clc
clear variables
close all

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

ecc_dir =[filedir, '/Cells'];
int_dir = [filedir, '/Intensity'];

cd(ecc_dir);
files = dir('*.csv');

average = zeros(numel(files),16);

for g = 1:numel(files)
    name.ecc = [num2str(g), '_area.csv'];
    name.int = [num2str(g), '_intensity.csv'];
    cd(ecc_dir);
    Cells = csvread(name.ecc,1,0);
    cd(int_dir);
    Int = csvread(name.int,1,0);
    average(g,1) = g;
    % eccentricity
    average(g,2) = mean(Cells(:,5));
    average(g,3) = std(Cells(:,5))/ sqrt(length(Cells(:,5))-1);
    % after BG subtraction
    average(g,4) = mean(Int(Int(:,1)<=10,2));
    average(g,5) = std(Int(Int(:,1)<=10,2))/ sqrt(length(Int(Int(:,1)<=10,2))-1);
    average(g,6) = mean(Int(Int(:,1)>=40,2));
    average(g,7) = std(Int(Int(:,1)>=40,2))/ sqrt(length(Int(Int(:,1)>=40,2))-1);
    
    % before BG subtraction
    average(g,10) = mean(Int(Int(:,1)<10,3));
    average(g,11) = std(Int(Int(:,1)<10,3))/ sqrt(length(Int(Int(:,1)<10,3))-1);
    average(g,12) = mean(Int(Int(:,1)>=40,3));
    average(g,13) = std(Int(Int(:,1)>=40,3))/ sqrt(length(Int(Int(:,1)>=40,3))-1);
    
    %number cells
    average(g,16) = length(Cells(:,5));
end

average(:,8) = average(:,6)./average(:,4);
average(:,9) = average(:,8) .* sqrt(average(:,5).*average(:,5)./average(:,4)./average(:,4)...
    + average(:,7).*average(:,7)./average(:,6)./average(:,6));

average(:,14) = average(:,12)./average(:,10);
average(:,15) = average(:,14) .* sqrt(average(:,11).*average(:,11)./average(:,10)./average(:,10)...
    + average(:,13).*average(:,13)./average(:,12)./average(:,12));
headers = {'Image', 'Eccentricity','sem', '0-10afterBG','sem', '40-90afterBG','sem', 'AsymmetryafterBG', 'sem'...
    '0-10beforeBG','sem', '40-90beforeBG','sem', 'AsymmetrybeforeBG', 'sem', 'Cells'};

cd([filedir, '/SummaryIntensity']);
csvwrite_with_headers('Ecc_vs_Asymmetry.csv',average,headers);
cd(currdir);

clear variables
clc
