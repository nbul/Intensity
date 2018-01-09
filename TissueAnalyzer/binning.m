clc
clear variables
close all

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

files = dir('*.csv');

data = [0,0,0,0];
%% Combining data
for i=1:numel(files)
    data_temp = csvread(files(i).name, 1,0);
    data = [data; data_temp];
end

data(1,:) = [];
data = sortrows(data,1);

%% Binning
step = 10;
data_binned = zeros(90/step,7);
for i=step:step:90
    data_temp2 = zeros(1,4);
    counter = 0;
    for k=1:length(data)
        if data(k,1)<i && data(k,1)>=i-step
            counter = counter + 1;
            data_temp2(counter,:) = data(k,:);           
        end
    end   
    data_binned(i/step,1) = mean(data_temp2(:,1));
    data_binned(i/step,2) = std(data_temp2(:,1)) / size(data_temp2(:,1),1);
    data_binned(i/step,3) = mean(data_temp2(:,2));
    data_binned(i/step,4) = std(data_temp2(:,2)) / size(data_temp2(:,1),1);
    data_binned(i/step,5) = mean(data_temp2(:,3));
    data_binned(i/step,6) = std(data_temp2(:,3)) / size(data_temp2(:,1),1);
    data_binned(i/step,7) = size(data_temp2(:,1),1);
end

csvwrite('DataBinned.csv',data_binned);
cd(currdir);
clc
clear variables
close all
