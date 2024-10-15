%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);
files_tif = dir('*_cp_masks.png');
data_pulled = zeros(1,7);

for loop = 1:length(files_tif)
    Mask = imread(files_tif(loop).name);
    Mask2 = Mask;
    Mask2(ismember(Mask, union(Mask(:, [1 end]), Mask([1 end], :)))) = 0;
    s = regionprops(Mask2, {'Centroid', 'PixelIdxList'});
    Mask3 = Mask2;
    A = unique(Mask2);
    for i = 2:length(unique(Mask2))
        number = A(i);
        Mask3(Mask2(:,:)==number) = i-1;
    end
    x = sscanf(files_tif(loop).name,'%d_cp_masks.png');
    HP1 = imread([num2str(x), '.tif']);
    data_image2 = regionprops('table', Mask3, HP1, 'MeanIntensity','Area','Eccentricity');
    HP1_3 = imgaussfilt(HP1,2);
    T = graythresh(HP1_3);
    HP1_2 = imbinarize(HP1_3, T*1.7);
    %HP1_2 = imbinarize(HP1_3, 'adaptive', 'Sensitivity', 0.4);
    %montage({HP1_2,imadjust(double(HP1)/256/256)})
    foci_all = zeros(length(unique(Mask3))-1,3);
    for i = 2:length(unique(Mask3))
        temp = Mask3 * 0;
        temp(Mask3(:,:) == i-1) = 1;
        temp2 = imbinarize(temp) .* HP1_2;
        cc = bwconncomp(temp2);
        foci = struct2table(regionprops(cc,HP1,'Area','MeanIntensity'));
        foci_all(i-1,:) = [size(foci,1),table2array(mean(foci,1))];
    end
    data_image3 = array2table(foci_all);
    data_image3.Properties.VariableNames = {'FociNumber','FociArea','FociIntensity'};
    data_all = [data_image2,data_image3];
    data_pulled = [data_pulled; ones(size((data_all),1),1) * x, table2array(data_all)];
end
data_pulled(1,:) = [];
data_table = array2table(data_pulled);
data_table.Properties.VariableNames = {'Image', 'Area',...
    'Eccentricity', 'MeanIntensity',...
'FociNumber','FociArea','FociIntensity'};
writetable(data_table, 'summary.csv');
cd(currdir);
clc
clear variables;