cd(oib_dir);

%% read image
Number1 = [num2str(k),'.oib'];

I=bfopen(Number1);
counter = 0;
Series = I{1,1};
seriesCount = size(Series, 1); %display size to check type of file
%% collecting write frames only
for i=1:seriesCount
    metadata = Series{i, 2};
    if contains(metadata,channel) == true
      counter = counter + 1;
      
    end
end