Total_bg = zeros(numel(b_valid),1);
Total = zeros(numel(b_valid),1);

mI_bg = (mI -  BG);


for i=1:numel(b_valid)
    Total(i) = mI(:,i)*length(b_valid{i});
    Total_bg(i) = mI_bg(:,i)*length(b_valid{i});
end



T = table((1:numel(b_valid))', mI', Total,...
     mI_bg', Total_bg,cell_data(:,3), cell_data(:,4));
T.Properties.VariableNames = {'cell', 'Mean', 'Total',...
    'MeanBG', 'TotalBG', 'Area','Eccentricity'};

writetable(T, [num2str(k),'_result.csv']);

averagedata(k,1) = k;
averagedata(k,2:2:end-2) = mean(T{:,2:end},1);
averagedata(k,3:2:end-2) = std(T{:,2:end},1);
averagedata(k,end-1) = BG;
averagedata(k,end) = numel(b_valid);

pulled = [pulled; [k*ones(numel(b_valid),1), (1:numel(b_valid))', ...
    mI', Total, mI_bg', Total_bg cell_data(:,3), cell_data(:,4)]];
