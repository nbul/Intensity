msize = zeros(numel(b_valid),1);

mI = zeros(length(zsection), numel(b_valid));

for t=1:numel(b_valid)
    mask = zeros(im_x,im_y);
    for i=1:length(b_valid{t})
        mask(b_valid{t}(i,1),b_valid{t}(i,2)) = 1;
    end
    mask = imdilate(mask, [se90I se0I]);
    msize(t) = sum(mask(mask>0));
    for j=1:length(zsection)
        mask2 = double(mask) .* double(zsection{j});
        mI(j,t) = mean(mask2(mask2>0));
    end
end

image2 = figure;
for t=1:numel(b_valid)
    plot(mI(:,t),'-', 'LineWidth', 2);
    hold on;
end

cd(Dist_dir);
image_filename = [num2str(k),'_distributions.tif'];
print(image2, '-dtiff', '-r150', image_filename);
close all

