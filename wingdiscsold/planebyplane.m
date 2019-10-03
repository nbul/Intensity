msize = zeros(numel(b_valid),1);

mI = zeros(1, numel(b_valid));

for t=1:numel(b_valid)
    mask = zeros(im_x,im_y);
    for i=1:length(b_valid{t})
        mask(b_valid{t}(i,1),b_valid{t}(i,2)) = 1;
    end
    mask = imdilate(mask, [se90I se0I]);
    mask2 = double(mask) .* double(Projection);
    mI(t) = mean(mask2(mask2>0));
end

%% Background
maskBG = zeros(im_x,im_y);
for t=1:numel(b_valid)   
    for i=1:length(b_valid{t})
        maskBG(b_valid{t}(i,1),b_valid{t}(i,2)) = 1;
    end  
end

maskBG = imdilate(maskBG, [se90I se0I]);
maskBG = imclearborder(imcomplement(maskBG));
mask3 = double(maskBG) .* double(Projection);
BG = mean(mask3(mask3>0));