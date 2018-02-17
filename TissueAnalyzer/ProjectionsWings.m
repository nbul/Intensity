%% Clear all and initial parameters
clc
clear variables
close all

%% Determening paths and setting folders
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
cd(filedir);

genotype = 'CD8';

parameters = inputdlg({'Genotype'},'Parameters',1,{num2str(genotype)});
% Redefine extension
genotype = parameters{1};

%Folder with CD8
if exist([filedir, '/', genotype],'dir') == 0
    mkdir(filedir, ['/', genotype]);
end
exp_dir = [filedir, '/', genotype];

%Folder with noCD8
if exist([filedir, '/no', genotype],'dir') == 0
    mkdir(filedir, ['/no', genotype]);
end
contr_dir = [filedir, '/no', genotype];

%% folders for projections
cd(exp_dir);
if exist([exp_dir, '/borders'],'dir') == 0
    mkdir(exp_dir, '/borders');
end
expb_dir = [exp_dir, '/borders'];

if exist([exp_dir, '/tifs_original'],'dir') == 0
    mkdir(exp_dir, '/tifs_original');
end
expi_dir = [exp_dir, '/tifs_original'];

%Folder with originals
if exist([exp_dir, '/originals'],'dir') == 0
    mkdir(exp_dir, '/originals');
end
expo_dir = [exp_dir, '/originals'];

cd(contr_dir);
if exist([contr_dir, '/borders'],'dir') == 0
    mkdir(contr_dir, '/borders');
end
contrb_dir = [contr_dir, '/borders'];

if exist([contr_dir, '/tifs_original'],'dir') == 0
    mkdir(contr_dir, '/tifs_original');
end
contri_dir = [contr_dir, '/tifs_original'];

if exist([contr_dir, '/originals'],'dir') == 0
    mkdir(contr_dir, '/originals');
end
contro_dir = [contr_dir, '/originals'];

cd(filedir);
files = dir('*.oib');
%% Projections
for i=1:(numel(files)/2)
    name.original = [num2str(i), '.oib'];
    name.control = [num2str(i), '_c.oib'];
    original = bfopen(name.original);
    control = bfopen(name.control);
    
    Series.original = original{1,1};
    Series.control = control{1,1};
    seriesCount = size(Series.original, 1); %display size to check type of file
    CadCD8.average = zeros(size(Series.original{1,1},1),size(Series.original{1,1},2));
    CadCD8.max = zeros(size(Series.original{1,1},1),size(Series.original{1,1},2));
    CadnoCD8.average = zeros(size(Series.control{1,1},1),size(Series.control{1,1},2));
    CadnoCD8.max = zeros(size(Series.control{1,1},1),size(Series.control{1,1},2));
    
    for plane = 1:(seriesCount/2)
        CadCD8.average = plus(CadCD8.average,double(Series.original{plane*2-1,1}));
   
        CadnoCD8.average = plus(CadnoCD8.average,double(Series.control{plane*2-1,1}));
        
    end 
    for plane = 2:(seriesCount/2-1)
         CadCD8.max = max(CadCD8.max, double(Series.original{plane*2-1,1}));
         CadnoCD8.max = max(CadnoCD8.max, double(Series.control{plane*2-1,1}));
    end
    
    CadCD8.average = uint16(CadCD8.average/seriesCount);
    CadnoCD8.average = uint16(CadnoCD8.average/seriesCount);
    
    CadCD8.background = imopen(CadCD8.max,strel('disk',50));
    CadCD8.max = CadCD8.max - CadCD8.background;
    CadCD8.max = imadjust(uint16(CadCD8.max));
    
    CadnoCD8.background = imopen(CadnoCD8.max,strel('disk',50));
    CadnoCD8.max = CadnoCD8.max - CadnoCD8.background;
    CadnoCD8.max = imadjust(uint16(CadnoCD8.max));
    
   
    cd(expb_dir);
    imwrite(CadCD8.max, [num2str(i), '.tif']);
    cd(contrb_dir);
    imwrite(CadnoCD8.max, [num2str(i), '.tif']);
    cd(expi_dir);
    imwrite(CadCD8.average, [num2str(i), '.tif']);
    cd(contri_dir);
    imwrite(CadnoCD8.average, [num2str(i), '.tif']);
    cd(filedir);
    movefile([num2str(i), '.oib'], expo_dir);
    movefile([num2str(i), '_c.oib'], contro_dir);
end

cd(currdir);
close all
clear variables
clc