curve = struct([]);
gof = struct([]);
zposition = zeros(numel(b_valid),1);
width1 = zeros(numel(b_valid),1);
pvalue = zeros(numel(b_valid),1);
Intensity = zeros(numel(b_valid),1);
Total = zeros(numel(b_valid),1);


dist_length = bin_size * (numel(zsection)-1);
binrange = 0 : bin_size : dist_length;
bincenter=binrange(1:(end-1)) + bin_size/2;
mI_bg = (mI -  mean(zsection{1}(:)));

image3 = figure;

for i=1:numel(b_valid)
    [curve{i},gof{i}] = fit(binrange',mI_bg(:,i),'gauss1');
    zposition(i) = curve{i}.b1';
    width1(i) = sqrt(curve{i}.c1*curve{i}.c1/2)';
    pvalue(i) = gof{i}.rsquare';
    Intensity(i) = mean(mI_bg(:,i));
    Total(i) = trapz(mI_bg(:,i))*msize(i);
    subplot(4, ceil(numel(b_valid)/4),i);
    plot(binrange', mI_bg(:,i), 'o', binrange', curve{i}(binrange'));
    title(num2str(gof{i}.rsquare));
    hold on;
end

cd(Fit_dir);
image_filename = [num2str(k),'_distribution_fit.fig'];
savefig(image3,image_filename);
close all



T = table((1:numel(b_valid))', Intensity, Total, zposition,...
    width1, pvalue, cell_data(:,3), cell_data(:,4));
T.Properties.VariableNames = {'cell', 'Mean', 'Total', 'Zposition', ...
    'Width', 'pvalue', 'Area','Eccentricity'};
Tclean = T(T.pvalue>cutoff,:);


writetable(T, [num2str(k),'_fit.csv']);
writetable(Tclean, [num2str(k),'_fit_cutoff.csv']);

averagedata(k,1) = k;
averagedata(k,2:2:end) = mean(T{:,2:end},1);
averagedata(k,3:2:end) = std(T{:,2:end},1);

averagedata2(k,1) = k;
averagedata2(k,2:2:end) = mean(Tclean{:,2:end},1);
averagedata2(k,3:2:end) = std(Tclean{:,2:end},1);

pulled = [pulled; [k*ones(numel(b_valid),1), (1:numel(b_valid))', ...
    Intensity, Total, zposition,...
    width1, pvalue, cell_data(:,3), cell_data(:,4)]];

pulled2 = [pulled2; [k*ones(height(Tclean),1),table2array(Tclean)]];